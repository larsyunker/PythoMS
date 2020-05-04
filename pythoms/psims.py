"""Classes for interacting with HUPO-PSI-MS controlled variable names and descriptions"""
import os
import obonet
import textwrap
import urllib
from xml.etree import ElementTree
cv_param_def = None  # default state for loaded CV parameters


def find_obos(paths: list):
    """
    locates any *.obo files in the specified directory

    :param paths: paths to search
    :return: list of paths
    :rtype: list
    """
    locations = []
    for eachpath in paths:
        for root, dirs, files in os.walk(eachpath):
            for ind in files:
                if ind.endswith('.obo'):
                    locations.append(os.path.join(root, ind))
    return locations


def newest_version(paths: list):
    """
    Select the newest possible obo version in the supplied paths

    :param paths: list of paths
    :return: selected path
    :rtype: str
    """
    ver = 0.
    out = ''
    for loc in paths:
        hndl = open(loc, 'rt')
        lines = hndl.readlines()
        hndl.close()
        for line in lines:
            if line.startswith('format-version:'):
                if float(line.split()[1]) > ver:
                    ver = float(line.split()[1])
                    out = loc
                break
    return out


def verify_obo(path):
    """
    Verifies that the provided obo is a HUPO-PSI obo

    :param path: path to obo file
    :return: whether the file schema is HUPO-PSI
    :rtype: bool
    """
    loaded = obonet.read_obo(path)
    return all([
        loaded.graph['ontology'] == 'ms',  # MS ontology
        'publisher: HUPO Proteomics Standards Initiative Mass Spectrometry Standards Working Group and HUPO Proteomics '
        'Standards Initiative Proteomics Informatics Working Group' in loaded.graph['remark'],  # correct publisher
    ])


def find_local_obo(custom_paths: list = None):
    """
    Finds local paths to obo files

    :param custom_paths: additional paths besides default
    :return: path to newest, local obo or None
    :rtype: str
    """
    paths = [
        os.getcwd(),  # by default, check CWD
    ]
    if custom_paths is not None:
        paths.extend(paths)
    possible = find_obos(paths)  # retrieve possible locations
    verified = [path for path in possible if verify_obo(path) is True]  # verify correct paths
    if len(verified) == 0:
        return None  # if there are no valid possibilities, return None
    return newest_version(verified)


# url for accessing the most up-to-date version of the psi-ms obo file
OBOURL = 'https://raw.githubusercontent.com/HUPO-PSI/psi-ms-CV/master/psi-ms.obo'
local_obo = find_local_obo()  # try to find a local obo and use instead
if local_obo is not None:
    OBOURL = local_obo


def interpret_term_string(string: str):
    """
    Interprets a term string and returns a dictionary of separated values

    :param string: string to parse
    :return: parsed string
    :rtype: dict
    """
    lst = string.split('\n')  # split into a line-separated list
    out = {}
    for line in lst:
        if line == '':  # ignore empty lines
            continue
        name, description = line.split(': ')  # split into name, description
        out[name] = description
    return out  # return the dictionary


def stringtodigit(string: str):
    """
    attempts to convert a unicode string to float or integer

    :param str string: String to attempt conversion on
    :return:
    """
    try:
        value = float(string)  # try converting to float
    except ValueError:
        try:
            value = int(string)  # try converting to int
        except ValueError:
            value = string  # otherwise keep as unicode
    return value


class CVParam(object):
    _value = None
    unitacc = None
    _definition = None
    _xref = None
    _comment = None
    _is_a = []
    _synonym = None
    _obsolete = False
    _relationship = []
    _property_value = None
    _name = None

    def __init__(self,
                 id: str = None,
                 name: str = None,
                 definition: str = None,
                 xref: str = None,
                 comment: str = None,
                 is_a: str = None,
                 obsolete: bool = False,
                 synonym: str = None,
                 relationship: str = None,
                 property_value: str = None,
                 value=None,
                 **kwargs,  # catch in case they implement anything new
                 ):
        """
        Class to define a controlled vocabulary object. If values are not provided, they will be retrieved automatically
        from the definitions in the psi-ms.obo file.

        :param id: MS id of the value
        :param name: name of the value
        :param definition: definition of the value
        :param xref: cross reference
        :param comment: general comment
        :param is_a: association key to other accession values
        :param obsolete: whether the accession key is obsolete
        :param synonym: synonym for the value
        :param value: any value associated with the CV parameter
        """
        # store values
        self.id = id
        self._name = name
        self._definition = definition
        self._xref = xref
        self._comment = comment
        self.is_a = is_a
        self._synonym = synonym
        self.obsolete = obsolete
        self.relationship = relationship
        self._property_value = property_value
        self.value = value

        if 'unitAccession' in kwargs:
            self.unitacc = kwargs['unitAccession']

    def _get_value_or_default(self, item):
        """
        Attempt to retrieve local value, and if not, return defined value

        :param item: item to retrieve
        :return: value
        """
        if getattr(self, item) is not None:
            return getattr(self, item)
        elif cv_param_def is not None:
            return getattr(cv_param_def[self.id], item)
        return None

    @property
    def name(self):
        return self._get_value_or_default('_name')

    @property
    def value(self):
        """value associated with the parameter (if any)"""
        return self._value

    @value.setter
    def value(self, value):
        if value == '':
            value = None
        if value is not None:
            self._value = stringtodigit(value)

    @property
    def unit(self):
        if self.unitacc is None:
            return None
        return cv_units[self.unitacc].name

    @property
    def definition(self):
        return self._get_value_or_default('_definition')

    @property
    def xref(self):
        return self._get_value_or_default('_xref')

    @property
    def comment(self):
        return self._get_value_or_default('_comment')

    @property
    def is_a(self):
        return self._get_value_or_default('_is_a')

    @is_a.setter
    def is_a(self, value):
        if value is not None:
            if type(value) != list:
                raise TypeError(f'is_a must be of type: list, provided: {type(value)}')
            self._is_a = list(value)

    @property
    def synonym(self):
        return self._get_value_or_default('_synonym')

    @property
    def obsolete(self):
        return self._get_value_or_default('_obsolete')

    @obsolete.setter
    def obsolete(self, value):
        if value is not None:
            self._obsolete = bool(value)

    @property
    def relationship(self):
        return self._get_value_or_default('_relationship')

    @relationship.setter
    def relationship(self, value):
        if value is not None:
            if type(value) != list:
                raise TypeError(f'is_a must be of type: list, provided: {type(value)}')
            self._relationship = list(value)

    @property
    def property_value(self):
        return self._get_value_or_default('_property_value')

    def __repr__(self):
        return f'{self.__class__.__name__}({self.id})'

    def __str__(self):
        return f'{self.id} {self.name}'


class CVParameterSet(object):
    def __init__(self,
                 **incoming
                 ):
        """
        Creates an object which allows access to a set of CV parameters

        :param incoming: key, kwargs sets for instantiation of CVParam objects
        """
        self.cv_values = {}
        for acc in incoming:
            self.cv_values[acc] = CVParam(
                id=acc,
                **incoming[acc],
            )

    def __repr__(self):
        return f'{self.__class__.__name__}({len(self.cv_values)})'

    def __str__(self):
        return f'{self.__class__.__name__} with {len(self.cv_values)} keys'

    def __getitem__(self, item):
        """:rtype: CVParam"""
        if item in self.cv_values.keys():
            return self.cv_values[item]
        else:
            try:
                item = self.acc_from_name(item)
                return self.cv_values[item]
            except KeyError:
                raise KeyError(f'No accession number or name could be found matching the provided item {item}')

    def __iter__(self):
        """:rtype: CVParam"""
        for key in self.cv_values.keys():
            yield self.cv_values[key]

    def __len__(self):
        return len(self.cv_values)

    def __contains__(self, item):
        if item in self.ids:  # defined accession key
            return True
        elif item in self.names:  # defined names
            return True
        elif item in self.__dict__:
            return True
        return False

    @property
    def names(self):
        """Defined names"""
        return {self.cv_values[key].name for key in self.cv_values}

    @property
    def ids(self):
        return {key for key in self.cv_values}

    def keys(self):
        """
        :return: defined accession keys
        :rtype: set
        """
        return self.ids

    def acc_from_name(self, name):
        """
        Finds the accession key from the name provided

        :param name: name of the value
        :return: accession key
        :rtype: str
        """
        name = name.lower()  # avoids case sensitivity
        for acc in self.cv_values:
            if self.cv_values[acc].name.lower() == name:
                return acc
        raise KeyError(f'No accession key was found matching the supplied name {name}.')

    def print_properties(self, key):
        """Prints the properties of the provided key"""
        cv_param_def.print_properties(key)  # print the detailed information
        if self.cv_values[key].value is not None:
            print(f'value: {self.cv_values[key].value}')  # followed by the value

    @classmethod
    def branch_to_dict(cls,
                       branch: ElementTree.Element,
                       search_children: bool = True,
                       ) -> dict:
        """
        Searches through children of a branch and retrieves cvparameter values and details.

        :param branch: branch to search
        :param search_children: whether to search child branches for additional cvparameters
        :return: dictionary of cv parameters
        """
        _cvparam_name = '{http://psi.hupo.org/ms/mzml}cvParam'
        out = {}
        for param_element in branch.findall(_cvparam_name):
            acc = param_element.attrib.get('accession')  # accession key
            out[acc] = {}
            for attribute, value in param_element.attrib.items():  # pull all the attributes
                if attribute != 'accession':
                    # attempt to convert to integer or float, keep as string otherwise
                    out[acc][attribute] = stringtodigit(value)
        if search_children is True:
            for child in branch:
                if child.tag != _cvparam_name:
                    out.update(
                        cls.branch_to_dict(
                            child,
                            search_children=True,
                        ),
                    )
        return out

    @classmethod
    def create_from_branch(cls,
                           branch: ElementTree.Element,
                           search_children: bool = True,
                           ) -> "CVParameterSet":
        """
        Creates a class instance from an XML branch

        :param branch: XML branch to interpret
        :param search_children: whether to search child branches for additional cvparameters
        :return: cv parameter set associated with that branch
        """
        return cls(
            **cls.branch_to_dict(
                branch=branch,
                search_children=search_children,
            )
        )


class CVParameterDefinitions(CVParameterSet):
    def __init__(self,
                 obo_path=OBOURL,
                 ):
        """
        Loads and interprets a PSI-MS obo file into a python-interpretable format.

        :param obo_path: file path or url to an obo file
        """
        CVParameterSet.__init__(self)
        try:
            self.obo_file = obonet.read_obo(obo_path)  # read the obo file
        except (FileNotFoundError, urllib.error.HTTPError):
            raise FileNotFoundError(f'An obo file could not be found at the provided path or URL: {obo_path}')
        if obo_path == OBOURL:  # remind the user to cite the pulication
            print('Data was read from PSI-MS, please cite DOI: 10.1093/database/bat009')
        self.format_version = self.obo_file.graph['format-version']
        self.data_version = self.obo_file.graph['data-version']
        for acc in self.obo_file:
            dct = self.obo_file.nodes[acc]
            if 'def' in dct:  # if the invalid key def is in the dictionary, convert and remove
                dct['definition'] = dct['def']
                del dct['def']
            self.cv_values[acc] = CVParam(
                id=acc,
                **dct,
            )

    def __repr__(self):
        return f'{self.__class__.__name__}(v{self.format_version})'

    def __str__(self):
        return f'{self.__class__.__name__} format version {self.format_version} with {len(self.cv_values)} keys'

    def print_properties(self, key):
        cvparam = self.__getitem__(key)
        print(f'ID: {cvparam.id}')
        print(f'name: {cvparam.name}')
        if cvparam.obsolete is True:
            print('This parameter is obsolete.')
        if cvparam.definition is not None:
            definition = "\n    ".join(textwrap.wrap(cvparam.definition))
            print(f'definition: {definition}')
        if cvparam.is_a is not None:
            print('is a:')
            for key in cvparam.is_a:
                print(f'    {str(self.__getitem__(key))}')
        if cvparam.relationship is not None:
            print('relationships: ')
            for rel in cvparam.relationship:
                context, key = rel.split(' ')
                print(f'    {context} {str(self.__getitem__(key))}')
        if cvparam.synonym is not None:
            print('synonyms: ')
            for syn in cvparam.synonym:
                print(f'    {syn}')


cv_param_def = CVParameterDefinitions()

unit_definitions = {
    # x units
    'MS:1000040': 'm/z',
    'UO:0000010': 'second',
    'UO:0000017': 'micrometer',
    'UO:0000018': 'nanometer',
    'UO:0000028': 'millisecond',
    'UO:0000031': 'minute',
    'UO:0000221': 'dalton',
    'UO:0000222': 'kilodalton',
    # y units
    'MS:1000131': 'number of detector counts',
    'MS:1000132': 'percent of base peak',
    'MS:1000814': 'counts per second',
    'MS:1000905': 'percent of base peak times 100',
    'UO:0000187': 'percent',
    'UO:0000269': 'absorbance unit',
    # other
    'MS:1000807': 'Th/s',
    'UO:0000002': 'mass unit',
    'UO:0000008': 'meter',
    'UO:0000012': 'kelvin',
    'UO:0000021': 'gram',
    'UO:0000027': 'degree Celsius',
    'UO:0000098': 'milliliter',
    'UO:0000106': 'hertz',
    'UO:0000110': 'pascal',
    'UO:0000112': 'joule',
    'UO:0000166': 'parts per notation unit',
    'UO:0000169': 'parts per million',
    'UO:0000175': 'gram per liter',
    'UO:0000185': 'degree',
    'UO:0000218': 'volt',
    'UO:0000228': 'tesla',
    'UO:0000266': 'electronvolt',
    'UO:0000268': 'volt per meter',
}

# object for automatic retrieval of unit names from accession values
cv_units = CVParameterSet(
    **{acc: {'name': unit_definitions[acc]} for acc in unit_definitions}
)
