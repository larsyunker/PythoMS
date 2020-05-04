"""
IGNORE:
CHANGELOG:
-
---2.7 building

to add:
    try to extract timepoints and tic from chromatogramList (x values are sorted, so this probably won't work)
IGNORE
"""
import sys
import os
import re
import zlib
import gzip
import base64
import struct
import subprocess
import logging
import pathlib
import warnings
from xml.etree import ElementTree
from tqdm import tqdm
import numpy as np
import scipy as sci
from random import random
from typing import Generator, Iterable, Union
from .spectrum import Spectrum
from .psims import CVParameterSet, stringtodigit
from .tome import resolution, locate_in_list, trimspectrum

# decoding formats for decoding mzML binary data array strings
decode_formats = {
    'MS:1000519': ['<', 'i'],  # signed 32-bit little-endian integer
    # 'MS:1000520':['',''], # [OBSOLETE] Signed 16-bit float
    'MS:1000521': ['<', 'f'],  # 32-bit precision little-endian floating point conforming to IEEE-754
    'MS:1000522': ['<', 'l'],  # Signed 64-bit little-endian integer
    'MS:1000523': ['<', 'd'],  # 64-bit precision little-endian floating point conforming to IEEE-754.
}

# prefix for xml elements
_xml_element_prefix = '{http://psi.hupo.org/ms/mzml}'

logger = logging.getLogger(__name__)


class BoundsError(Warning):
    """A warning class to handle bounds errors when integrating (used only by PyRSIR)"""

    def __init__(self):
        self.warned = {}

    def printwarns(self):
        """prints the number of warnings if merited"""
        if len(self.warned) > 0:
            logger.warning(
                'The following peaks exceeded the bounds of the spectrum n number of times:\n'
                f'{", ".join([f"{name}: {self.warned[name]}" for name in self.warned])}'
            )

    def warn(self, name, intstart, intend, mzstart, mzend):
        """warns the user if there was a mismatch"""
        if name not in self.warned:
            logger.warning(
                f'The peak "{name}" ({intstart}-{intend}) is outside of the bounds of the spectrum being summed '
                f'm/z {mzstart:.1f}-{mzend:.1f}\n'
            )
            self.warned[name] = 1
        else:
            self.warned[name] += 1


def branch_attributes(branch: ElementTree.Element):
    """
    Pulls all the attributes of an xml.dom.minidom xml branch.
    These are generally things like index, id, etc.

    :param xml.dom.minidom branch: An xml.dom.minidom object.
    :return: A dictionary of attributes with each key being the attribute name and its value being the value of that
        attribute.
    :rtype: dict

    **Notes**

    The script will attempt to convert any values to float or
    integer in order to reduce TypeErrors when trying to use
    the extracted values.
    """
    return {key: stringtodigit(val) for key, val in branch.attrib.items()}


def file_present(filepath):
    """checks for the presence of the specified file or directory in the current working directory"""
    tf = os.path.isfile(filepath)  # look for file first
    if tf is False:  # if file cannot be found, look for directory
        tf = os.path.isdir(filepath)
    return tf


def decodeformat(p: CVParameterSet, speclen: int):
    """
    Determines the decode format from the accession parameter

    :param p: extracted CVParamterSet of the data array
    :param speclen: length of the spectrum (retrievable from the XML file)
    :return: decode format
    :rtype: str
    """
    for key in set(decode_formats) & p.keys():  # find the set combination of the possibilities
        return f'{decode_formats[key][0]}{speclen}{decode_formats[key][1]}'  # create the decode format


def associate_data_type(cv_params: CVParameterSet):
    """
    Associates a parameter set with a numpy data type

    :param cv_params: set of cv parameters associated with a list of values
    :return: numpy data type
    """
    param_map = {
        '523': np.float64,
        '521': np.float32,
        '522': np.int64,
        '519': np.int32,
        '520': np.float16,
    }
    for key, dtype in param_map.items():
        if f'MS:1000{key}' in cv_params:
            return dtype
    raise ValueError('a data type could not be determined from this parameter set')


def gettext(nodelist):
    """gets text from a simple XML object"""
    rc = []
    for node in nodelist:
        if node.nodeType == node.TEXT_NODE:
            rc.append(node.data)
    return ''.join(rc)


def spectrum_array(spectrum: ElementTree.Element) -> np.ndarray:
    """
    Extracts and converts binary data to a numpy ndarray.

    :param spectrum: A spectrum branch element. This element is expected to have two child nodes containing
        binaryDataArrays.
    """
    # spectrum length (defined in the spectrum attributes)
    speclen = int(spectrum.attrib.get('defaultArrayLength'))
    out = np.ndarray((2, speclen))
    # iterate over the binary data arrays
    for ind, bda in enumerate(_findall_ele('binaryDataArrayList', 'binaryDataArray', parent=spectrum)):
        p = CVParameterSet.create_from_branch(bda)  # grab cvparameters

        # pull the binary string
        binary_string = _find_ele('binary', parent=bda).text

        # decode the string
        decoded = base64.standard_b64decode(binary_string)

        # if the string is compressed, decompress
        if 'MS:1000574' in p:
            decoded = zlib.decompress(decoded)

        # retrieve the data from buffer with the associated format
        out[ind] = np.frombuffer(
            decoded,
            associate_data_type(p)
        )

    return out


def get_element_units(spectrum: ElementTree.Element):
    """
    Retrieves the units of the provided element. Assumes that the element has binaryDataArrayList and binaryDataArray
    children which each have a unit-type controlled variable.

    :param spectrum: XML spectrum Element
    :return: x and y units corresponding to the spectrum
    """
    units = []
    # iterate over the binary data arrays
    for bda in _findall_ele('binaryDataArrayList', 'binaryDataArray', parent=spectrum):
        p = CVParameterSet.create_from_branch(bda)  # grab cvparameters
        for cv in p:
            if cv.unit is not None:
                units.append(cv.unit)
    return units


def extract_spectrum(spectrum: ElementTree.Element, units: bool = False) -> list:
    """
    Extracts and converts binary data to two lists.

    :param spectrum: A spectrum branch element. This element is expected to have two child nodes containing
        binaryDataArrays.
    :param units: whether to extract the units from the spectrum
    :return:
    """
    # todo go through all usages and refactor to array handling
    out = spectrum_array(spectrum).tolist()
    if units is not False:  # extends the units onto out
        out.extend(get_element_units(spectrum))
    return out


def pw_convert(filename,
               bit=64,
               compression=True,
               gzip=True,
               verbose=True,
               out_directory=None
               ):
    """
    Runs msconvert.exe from ProteoWizard to convert Waters .RAW format to .mzXML
    which can then be parsed by python.

    module requirements: os, subprocess, sys

    ProteoWizard must be installed for this script to function.
    go to
    http://proteowizard.sourceforge.net/downloads.shtml
    to download

    This script assumes that the ProteoWizard is installed under either
    c:\\program files\\proteowizard
    or
    c:\\program files (x86)\\proteowizard

    If you use this python script to convert to mzML, you should cite the paper of the folks who wrote the program
    Chambers, M.C. Nature Biotechnology 2012, 30, 918-920
    doi 10.1038/nbt.2377

    :param filename: file name to convert
    :param bit: floating point bit precision (32 or 64)
    :param compression: enable zlib compression of data in output file (saves space but increases processing time)
    :param gzip: enable gzip compression of output file (saves disk space but increases processing time)
    :param verbose: verbose for subprocess call
    :param out_directory: optional output directory (if not specified, the file will be saved to the same directory as
        the data)
    :return: the file path for the generated file
    """

    def find_all(fname, path):
        """
        Finds all files of a given name within a specified directory.
        Adapted from http://stackoverflow.com/questions/1724693/find-a-file-in-python

        Module dependancies: os
        """
        locations = []
        for root, dirs, files in os.walk(path):
            if fname in files:
                locations.append(os.path.join(root, fname))
        return locations

    if sys.platform != 'win32':
        raise OSError(
            'The function that converts to mzML is limited to Windows operating systems.\n'
            'You can manually convert to *.mzML using the proteowizard standalone package '
            'and supply that mzML file to this script')
    locs = []
    for val in ['c:\\program files\\proteowizard',
                'c:\\program files (x86)\\proteowizard']:  # searches for msconvert.exe in expected folders
        locs.extend(find_all('msconvert.exe', val))

    if len(locs) == 0:  # if script cannot find msconvert.exe
        raise IOError(
            'The python script could not find msconvert.exe\n'
            'Please ensure that ProteoWizard is installed in either:\n'
            'c:\\program files\\proteowizard\nor\nc:\\program files (x86)\\proteowizard')

    filename = pathlib.Path(filename)
    if out_directory is None:
        out_directory = filename.parent

    if bit not in [32, 64]:
        raise ValueError(
            f'an invalid floating point precision was specified"{bit}".')

    callstring = " ".join([
        f'{locs[-1]} "{filename}"',  # main call
        f'-o "{out_directory}"',  # output directory
        '--mzML',
        '--gzip' if gzip else '',  # gzip compression
        f'--{bit}',  # floating point precision
        '--verbose' if verbose else '',  # verbose mode
    ])

    out_exten = f'.mzML{".gz" if gzip else ""}'

    logger.info(f'Generating mzML file from {filename}')
    subprocess.call(callstring)
    logger.info('conversion DONE')
    filename.with_suffix(f'.mzML{".gz" if gzip else ""}')


def fix_extension(fn):
    """tries to fix invalid file extensions"""
    oopsx = {'.mzm': 'l', '.mz': 'ml', '.m': 'zml', '.': 'mzml'}  # incomplete mzml extensions
    oopsr = {'.ra': 'w', '.r': 'aw', '.': 'raw'}  # incomplete raw extionsions
    oopsg = {'.mzml.g': 'z', '.mzml.': 'gz', '.mzml': '.gz', '.mzm': 'l.gz', '.mz': 'ml.gz', '.m': 'zml.gz',
             '.': 'mzml.gz'}  # incomplete gz extensions
    # looks for missing extensions first
    if file_present(fn + '.mzml.gz') is True:
        return fn + '.mzml.gz'
    if file_present(fn + '.mzml') is True:
        return fn + '.mzml'
    for key in oopsg:  # tries to complete mzml.gz shortenings
        if fn.lower().endswith(key) is True:
            if file_present(fn + oopsg[key]) is True:
                return fn + oopsg[key]
    for key in oopsx:  # tries to complete mzml shortenings
        if fn.lower().endswith(key) is True:
            if file_present(fn + oopsx[key]) is True:
                return fn + oopsx[key]
    for key in oopsr:  # tries to complete raw shortenings
        if fn.lower().endswith(key) is True:
            if file_present(fn + oopsr[key]) is True:
                return fn + oopsr[key]
    if file_present(fn + '.raw') is True:  # finally looks for raw file
        return fn + '.raw'
    raise FileNotFoundError(f'The file {fn} could not be located in the current working directory')


_waters_id_re = re.compile('function=(?P<fn>\d+)\sprocess=(?P<proc>\d+)\sscan=(?P<scan>\d+)')
_agilent_id_re = re.compile('scanId=(?P<scan>\d+)')


def fps(element: ElementTree.Element):
    """
    Extracts function #, process #, and scan # from the idstring of a spectrum branch

    :param element: XML element to retrieve from
    :return:
    """
    # pull id string from scan attribute
    idstring = element.attrib['id']
    match = _waters_id_re.match(idstring)
    if match is not None:
        return (
            int(match.group('fn')),
            int(match.group('proc')),
            int(match.group('scan')),
        )
    else:
        # todo create generalized catch
        match = _agilent_id_re.match(idstring)
        return (
            1,
            None,
            int(match.group('scan'))
        )


def scan_properties(parameters: Union[CVParameterSet, ElementTree.Element]):
    """
    Determines the scan properties of the provided spectrum.

    :param parameters: CVParam parameters
    :return:
    """
    """determines the scan properties of the provided spectrum"""
    mstypes = {  # ms accession keys and their respective names (for spectrum identification)
        'MS:1000928': 'calibration spectrum',
        'MS:1000294': 'mass spectrum',
        'MS:1000322': 'charge inversion mass spectrum',
        'MS:1000325': 'constant neutral gain spectrum',
        'MS:1000326': 'constant neutral loss spectrum',
        'MS:1000328': 'e/2 mass spectrum',
        'MS:1000341': 'precursor ion spectrum',
        'MS:1000343': 'product ion spectrum',
        'MS:1000579': 'MS1 spectrum',
        'MS:1000580': 'MSn spectrum',
        'MS:1000581': 'CRM spectrum',
        'MS:1000582': 'SIM spectrum',
        'MS:1000583': 'SRM spectrum',
    }
    othertypes = {  # other accession keys (non-MS)
        'MS:1000620': 'PDA spectrum',
        'MS:1000804': 'electromagnetic radiation spectrum',
        'MS:1000805': 'emission spectrum',
        'MS:1000806': 'absorption spectrum',
    }
    out = {}
    if isinstance(parameters, CVParameterSet):  # handed a cvparam class object (expected)
        p = parameters
    else:  # handed a tree or branch (generate the cvparam class object)
        p = CVParameterSet.create_from_branch(parameters)
    for acc in p.keys() & mstypes.keys():  # check for ms spectrum
        out['acc'] = acc  # accession code
        out['name'] = mstypes[acc]  # name of spectrum
        out['type'] = 'MS'  # it is a mass spectrum
        out['level'] = p['MS:1000511'].value  # ms level
        out['window'] = [p['MS:1000501'].value, p['MS:1000500'].value]  # scan window
        if 'MS:1000129' in p:  # negative scan
            out['mode'] = '-'
        elif 'MS:1000130' in p:  # positive scan
            out['mode'] = '+'
        if 'MS:1000827' in p:  # if there is an isolation window target m/z
            out['target'] = p['MS:1000827'].value
        # if MSn > 2, not sure how to handle this (will have to be hard coded later as I have no examples)
        elif out['level'] > 2:
            raise ValueError(
                'This script has not been coded to handle MSn > 2, please contact the author of the class')
        return out

    for acc in p.keys() & othertypes.keys():  # if the scan is something else
        out['acc'] = acc  # accession code
        out['name'] = othertypes[acc]  # name of spectrum
        if 'MS:1000804' in p:  # if it is a UV-Vis
            out['type'] = 'UV'
        else:  # other other type (not handled by script)
            raise KeyError(
                'The script has not been coded to handle spectra types other than MS and UV-Vis. '
                'Please contact the authors to get this functionality included.')
        return out


def _find_ele(*paths, parent: ElementTree.Element = None) -> ElementTree.Element:
    """
    Performs a find on the element tree.

    Each element name must be prefixed by "{http://psi.hupo.org/ms/mzml}", so this function saves that step.

    :param paths: path names
    :param parent: parent to search (if not provided, root will be used)
    :return: element tree at location
    """
    return parent.find(
        "/".join([
            f'{_xml_element_prefix}{path}' for path in paths
        ])
    )


def _findall_ele(*paths, parent: ElementTree.Element = None) -> Iterable[ElementTree.Element]:
    """
    Performs a findall on the element tree.

    Each element name must be prefixed by "{http://psi.hupo.org/ms/mzml}", so this function saves that step.

    :param paths: path names
    :param parent: parent to search (if not provided, root will be used)
    :return: element tree at location
    """
    return parent.findall(
        "/".join([
            f'{_xml_element_prefix}{path}' for path in paths
        ])
    )


class mzML(object):
    def __init__(self,
                 filename: str,
                 verbose: bool = True,
                 precision: int = 64,
                 compression: bool = True,
                 gzip_file: bool = True,
                 obo: str = None,
                 ftt: bool = False,
                 **kwargs
                 ):
        """
        A class for loading and extracting data from an mzML file.

        :param str filename: The name of the mzML or mass spectrometric data file. Accepted file types are listed below,
            and this script can automatically convert some proprietary file types to mzML by calling ProteoWizard
            (see notes).
        :param bool verbose: Chatty enable or disable. It can be useful to enable this when processing large files or long
            acquisitions, as many of the methods have progress reporters.
        :param int precision: The floating point precision to use if converting to mzML. Default 64 (although this
            appears to have a minimal effect in the experience of the author). This can be set to 32 to decrease mzML
            file sizes.
        :param bool compression: Whether or not to compress the mzML files when converting. This can decrease file
            sizes at a slight cost in processing time.
        :param bool gzip: Whether or not to gzip the mzML files when converting. This substantially decreases file
            sizes (mass spectrometric data compresses very well when gzipped). This will slightly increase processing
            time.
        :param str obo: A specific path or URL to an *.obo file defining the accession keys used in mzML files. If this
            is not specified, the default accession URL will be used to download the required obo file. This should not
            be necessary normally, as most of the commonly encountered accession keys are hard-coded into this
            script. The script will raise an error if it encounters an undefined accession key.
        :param bool ftt:  Whether to run the function_timetic() method on initialization. This is useful if you require
            access to the total ion current and time lists for each function in the mzML file. This does increase file
            load times quite significantly (~6x slower).

        **Notes**

        An mzML file is a data format for mass spectrometric data which can be parsed by python (avoiding the pitfalls
        associated with the proprietary files usually generated by the mass spectrometers themselves). The mzML file
        structures are expected to conform to those outlined in the HUPO Proteomics Standards Working Group. More
        information can be found at https://raw.githubusercontent.com/HUPO-PSI/psi-ms-CV/master/psi-ms.obo

        If you wish to use the format conversion functionality of this script, you will need to download and install
        ProteoWizard, which can be found at http://proteowizard.sourceforge.net/

        """
        # store keyword settings
        self.verbose = verbose
        self.precision = precision
        self.compression = compression
        self.gzip_file = gzip_file
        self.obo = obo

        self.filename = self.check_for_file(filename)

        logger.info(f'Loading {self.filename} into memory')
        if self.filename.lower().endswith('.mzml.gz'):  # if mzml is gzipped
            handle = gzip.open(self.filename)  # unzip the file
        else:
            handle = self.filename
        try:
            self.tree = ElementTree.parse(handle)  # full mzML file
            self.root = self.tree.getroot()
        except:
            raise IOError(
                'The mzML file "%s" could not be loaded. The file is either unsupported, corrupt, or incomplete.' % self.filename)

        try:  # number of spectra
            nscans = _find_ele('mzML', 'run', 'spectrumList', parent=self.root)
            self.nscans = int(nscans.attrib['count'])
        except IndexError:  # no spectra
            self.nscans = 0
        try:
            ncrhoms = _find_ele('mzML', 'run', 'chromatogramList', parent=self.root)
            self.nchroms = int(ncrhoms.attrib['count'])
        except IndexError:
            self.nchroms = 0

        self.functions = {}
        for spectrum in self._spectra_elements:
            try:
                # try to retrieve function, process, and scan from attributes
                func, proc, scan = fps(spectrum)  # extract each value and convert to integer
            except ValueError:
                # if undefined, assume only one function
                func = 1
            if func not in self.functions:  # if function is not defined yet
                p = CVParameterSet.create_from_branch(spectrum)  # pull spectrum's cvparameters
                self.functions[func] = {
                    'sr': [int(spectrum.attrib.get('index')), None],  # the scan index range that the function spans
                    'nscans': 1,  # number of scans
                }
                self.functions[func].update(scan_properties(p))  # update with scan properties
            else:
                self.functions[func]['sr'][1] = int(
                    spectrum.attrib.get('index'))  # otherwise set the scan index range to the current index
                self.functions[func]['nscans'] += 1
        try:
            p = CVParameterSet.create_from_branch(spectrum)  # pull properties of final spectrum
            self.duration = p['MS:1000016'].value  # final start scan time
        except (KeyError, UnboundLocalError):  # if there are no spectra, set to None
            # todo figure out a catch to retrieve time from other sources (e.g. TIC)
            self.duration = None

        self._BE = BoundsError()  # load warning instance for integration
        self.ftt = False
        if ftt is True:
            self.function_timetic()

    def __str__(self):
        """The string that is returned when printed"""
        return f'{self.__class__.__name__} {self.nscans} spectra, {self.nchroms} chromatograms'

    def __repr__(self):
        """The representation that is returned"""
        return "%s('%s')" % (self.__class__.__name__, self.filename)

    def __len__(self):
        return self.nscans

    def __getitem__(self, ind):
        """retrieves a scan or summed scans"""
        if isinstance(ind, slice):  # if getitem is trying to slice
            """
            returns the summed scans with the supplied indicies
            slice will assume that the intended function is 1
            """
            if ind.start is None:  # no start
                start = 0
            else:
                start = ind.start
            if ind.stop is None:  # no stop
                stop = self.functions[1]['sr'][1]
            else:
                stop = ind.stop
            return self.sum_scans(start, stop, mute=True)

        elif type(ind) is int:  # scan index number
            """will return the spectrum of the scan index provided"""
            if ind < 0 or ind > self.nscans:
                raise IndexError("The scan index number #%d is outside of the mzML's scan index range (0-%d)" % (
                ind, self.nscans - 1))
            for spectrum in self._spectra_elements:
                attr = branch_attributes(spectrum)
                if attr['index'] == ind:
                    return extract_spectrum(spectrum)

        elif type(ind) is float:  # timepoint in function 1
            """float will assume the intended function was 1"""
            if ind < 0 or ind > self.duration:
                raise ValueError(
                    "The supplied time %.3f is outside of this file's time range (0 - %.3f)" % (ind, self.duration))
            ind = self.scan_index(ind)
            for spectrum in self._spectra_elements:
                attr = branch_attributes(spectrum)
                if attr['index'] == ind:
                    return extract_spectrum(spectrum)

    @property
    def _spectra_elements(self) -> Generator[ElementTree.Element, None, None]:
        """generator of spectra elements"""
        for element in _findall_ele('mzML', 'run', 'spectrumList', 'spectrum', parent=self.root):
            yield element

    @property
    def _chromatogram_elements(self) -> Generator[ElementTree.Element, None, None]:
        """generator of chromatogram elements"""
        for element in _findall_ele('mzML', 'run', 'chromatogramList', 'chromatogram', parent=self.root):
            yield element

    @property
    def _mute_tqdm(self) -> bool:
        """whether to mute tqdm"""
        return not self.verbose

    def foreachchrom(self, fn):
        """
        a decorator function that will apply the supplied function to every chromatogram in the mzml file
        the supplied function will be handed the chromatogram XML object as the first argument
        the decorated function will return a list of outputs of the supplied function where each index corresponds to a scan

        e.g.::
            loaded = mzML(filename)

            @loaded.foreachchrom
            def do_this(chrom):
                # extract the attributes using the mzML.attributes() method
                attr = loaded.attributes(chrom)
                return attr['id'] # return the name of the chromatogram

            do_this()

        """

        def foreachchrom(*args, **kwargs):
            """decorates the supplied function to run for every scan"""
            out = []
            msg = f'applying function to chromatogram'
            for chromatogram in tqdm(self._chromatogram_elements, desc=msg, disable=self._mute_tqdm):
                out.append(fn(chromatogram, *args, **kwargs))
            return out

        return foreachchrom

    def foreachscan(self, fn):
        """
        a decorator function that will apply the supplied function to every spectrum in the mzml file
        the supplied function will be handed the spectrum XML object as the first argument
        the decorated function will return a list of outputs of the supplied function where each index corresponds to a scan

        e.g.::

            loaded = mzML(filename)

            @loaded.foreachscan
            def do_this(scan):
                p = loaded.cvparam(scan) # pull spectrum's cvparameters
                sst = p['MS:1000016'] # start scan time
                x,y = loaded.extract_spectrum(scan,False) # extract the x,y spectrum
                # return the start scan time, x list, and y list
                return sst,x,y

            do_this() # do it
        """

        def foreachscan(*args, **kwargs):
            """decorates the supplied function to run for every scan"""
            out = []
            msg = 'applying function to spectrum'
            for spectrum in tqdm(self._spectra_elements, desc=msg, disable=self._mute_tqdm):
                out.append(fn(spectrum, *args, **kwargs))
            return out

        return foreachscan

    def associate_to_function(self, affin=None, level=None, dct=None):
        """
        Associates a given species to the appropriate function number
        in the mzML data file.

        **Parameters**

        affin: '+', '-', or 'UV'
            The affinity of the species. i.e. to positive mode,
            negative mode, or UV-Vis spectra respectively.

        level: *integer* or None
            If the species is found in an MS/MS function,
            the MS^n level can be specified here.

        dct: *dictionary*
            If details are known about the species' affinity,
            they can be provided in dictionary format.
            Specifically, this function looks for the keys:
            'function', 'affin', and 'level'.


        **Returns**

        function number: *integer*
            Returns the appropriate function number in which
            the given species should be found.


        **Notes**

        If nothing is provided to this method, it will return
        the integer 1 (assuming that the species will be found
        in the first function).

        """
        if dct is not None:  # if function was handed a dictionary
            if 'function' in dct:
                return dct['function']
            if 'affin' in dct:
                affin = dct['affin']
            if 'level' in dct:
                level = dct['level']

        if affin is None and level is None:
            return min(self.functions.keys())  # assume first function

        elif affin == 'UV':  # if UV-Vis affinity
            for fn in self.functions:  # determine which function is UV-Vis
                if self.functions[fn]['acc'] == 'MS:1000804':
                    return fn
            raise ValueError('There is no electromagnetic radiation spectrum function in this mzML file')

        elif affin in ['+', '-']:  # if affinity to mass spectrum
            levelcount = 0  # counter for number of matches to this affinity and level
            for fn in self.functions:
                if self.functions[fn]['type'] == 'MS':  # if fn is ms
                    if self.functions[fn]['mode'] == affin:  # if mode mathes
                        # if there is no level specified, assume 1
                        if level is None and self.functions[fn]['level'] == 1:
                            fnout = fn
                            levelcount += 1
                        elif self.functions[fn]['level'] == level:  # if level matches
                            fnout = fn
                            levelcount += 1
            if levelcount > 1:
                raise ValueError(
                    f"There affinity specification of mode: {affin}, level: '{level}' matches more than one function "
                    f"in the mzML file. \nTo process this species, be more specific in your level specification or "
                    f"assign it to a specific function number by adding a 'function' key to its dictionary.")
            return fnout
        else:  # if some other affinity
            raise ValueError('The specified affinity "%s" is not supported.' % affin)

    def auto_resolution(self, n=10, function=None, npeaks=4):
        """
        Attempts to automatically determine the resolution of the spectrometer
        that the provided mzML data file was recorded on.
        The method will find n random samples of the entire spectrum and
        calculate the resolution of each of those samples and return the
        average resolution.

        :param int n: The number of psuedo-random samples of the spectrum to determine
            the resolution of. Default 10.
        :param int function: The mzML function number to calculate the resolution of. Default 1.
        :param int npeaks: number of peaks to to try to find
        :return: Estimated resolution of the spectrum
        :rtype: float
        """
        def findsomepeaks(y):
            """roughly locates 4 peaks by maximum values in the spectrum and returns their index"""
            split = int(len(y) / npeaks)
            start = 0
            end = start + split
            splity = []
            for i in range(npeaks):
                splity.append(sci.asarray(y[start:end]))
                start += split
                end += split
            out = []
            for ind, section in enumerate(splity):
                maxy = max(section)
                if maxy == max(section[1:-1]):  # if max is not at the edge of the spectrum
                    out.append(sci.where(section == maxy)[0][0] + split * ind)
            return out

        if function is None:  # if no function is provided, use first
            function = self.associate_to_function()
        if self.functions[function]['type'] != 'MS':
            raise ValueError(
                'The auto_resolution function only operates on mass spectrum functions. '
                'Type of specified function %d: %s' % (function, self.functions[function]['type']))
        ranges = []  # list of scan intervals

        if self.functions[function]['nscans'] <= 20:  # if the number of scans is less than 20
            ranges = [[1, self.functions[function]['nscans']]]
        else:
            while len(ranges) < n:  # generate 10 pseudo-random intervals to sample
                ran = int(random() * self.functions[function]['nscans']) + self.functions[function]['sr'][0]
                if ran - 10 >= self.functions[function]['sr'][0] and ran + 10 <= self.functions[function]['sr'][1]:
                    ranges.append([ran - 10, ran + 10])
        summed = []
        msg = 'Estimating resolution of the instrument'
        for ind, rng in enumerate(tqdm(ranges, desc=msg, disable=self._mute_tqdm)):
            summed.append(  # sum those scans and append output
                self.sum_scans(
                    rng[0], rng[1],
                    function,
                    2,
                    mute=True
                )
            )
        res = []
        for spec in summed:  # calculate resolution for each scan range
            inds = findsomepeaks(spec[1])  # find some peaks
            for ind in inds:  # for each of those peaks
                res.append(resolution(spec[0], spec[1], ind, threshold=10))
        res = [y for y in res if y is not None]  # removes None values (below S/N)
        return sum(res) / len(res)  # return average

    def check_for_file(self, fn):
        """checks for the mzML file in the working directory and converts it if necessary"""
        # cast path-like to string to enable extension check
        if type(fn) is not str:
            fn = str(fn)

        valid = [  # supported extensions
            '.raw',
            '.mzml.gz',
            '.mzml',
        ]
        if fn.lower().endswith('.raw') is True:  # extension is raw
            if file_present(fn[:-4] + '.mzML.gz') is True:  # if corresponding gzipped mzml is present
                return fn[:-4] + '.mzML.gz'
            if file_present(fn[:-4] + '.mzML') is True:  # if corresponding mzml is present
                return fn[:-4] + '.mzML'
            # otherwise convert and return mzml
            return pw_convert(fn, self.precision, self.compression, self.gzip_file, verbose=self.verbose)
        elif file_present(fn) is True:  # if the specified file is present
            for exten in valid:  # checks for supported extensions
                if fn.lower().endswith(exten) is True:
                    return fn
            # otherwise asks user whether to continue
            if input(
                    'The extension of the supplied filename "%s" is unexpected and may not be supported.\n'
                    'Do you wish to proceed with file loading? [Y/N] ' % fn).lower() in ['y', 'yes']:
                return fn
            else:
                sys.exit('The user cancelled mzML loading.')
        else:
            fn = fix_extension(fn)  # try to fix extension
            if fn.lower().endswith('.raw') is True:  # convert if only raw file is found
                return pw_convert(fn, self.precision, self.compression, self.gzip_file, verbose=self.verbose)
            return fn

    def function_timetic(self):
        """
        extracts timepoints and tic lists for each function
        this function is separate from mzml contents because it would increase load times significantly (~6x)
        """
        msg = 'extracting timepoints and total ion current values'
        for function in self.functions:  # add timepoint and tic lists
            self.functions[function]['timepoints'] = []  # list for timepoints
            self.functions[function]['tic'] = []  # list for total ion current values
            if 'level' in self.functions[function] and self.functions[function]['level'] > 1:
                self.functions[function]['ce'] = []  # list for collision energies
        for spectrum in tqdm(self._spectra_elements, desc=msg, disable=self._mute_tqdm):
            attr = branch_attributes(spectrum)
            function, proc, scan = fps(spectrum)  # determine function, process, and scan numbers
            p = CVParameterSet.create_from_branch(spectrum)  # pull spectrum's cvparameters
            self.functions[function]['timepoints'].append(p['MS:1000016'].value)  # start scan time
            self.functions[function]['tic'].append(p['MS:1000285'].value)  # total ion current
            if 'MS:1000045' in p:
                self.functions[function]['ce'].append(p['MS:1000045'].value)  # collision energy
        self.ftt = True

    def integrate(self, name, start, end, x, y):
        """
        Integrates y values given x bounds in a paired set of lists (e.g. a m/z list and an intensity list)

        name: name of the peak being integrated (only used for warning purposes)
        start: float
            start x value
        end: float or None
            end x value
            None will return the nearest value to the provided start value
        x: list of x values
        y: list of y values (paired with x)

        returns: integral
        """
        if start > max(x) or start < min(x):  # check that start is within the m/z bounds
            self._BE.warn(name, start, end, min(x), max(x))
        if end is None:  # if only a start value is supplied, return closest to that value
            try:  # try to find the value in the list
                return y[locate_in_list(x, start)]
            except TypeError:  # if the value is not in the list, return 0
                return 0
        if end > max(x):  # check that end is within the m/z bounds
            self._BE.warn(name, start, end, min(x), max(x))
        else:
            l = locate_in_list(x, start, 'greater')
            r = locate_in_list(x, end, 'lesser')
            if l <= r:
                return sum(y[l:r])
            else:  # catch for if there are no values in the bounds
                return 0

    def pull_chromatograms(self):
        """
        Pulls mzML chromatograms

        returns:
        dictionary = {'chromatogram 1 id', 'chromatogram 2 id', ...}
        dictionary['chromatogram 1 id'] = {
        'x': list of x values
        'y': list of y values (paired with x)
        'xunit': unit of the x values
        'yunit': unit of the y values
        }
        """
        msg = 'extracting chromatogram'
        chroms = {}  # dictionary of chromatograms
        for chromatogram in tqdm(self._chromatogram_elements, desc=msg, disable=self._mute_tqdm):
            attr = branch_attributes(chromatogram)  # pull attributes
            x, y, xunit, yunit = extract_spectrum(chromatogram, True)  # extract x list, y list, and units
            chroms[attr['id']] = {'x': x, 'y': y, 'xunit': xunit, 'yunit': yunit}
        return chroms

    def pull_species_data(self, sp, sumspec=False):
        """
        Extracts integrated data at every timepoint for all species specified in the sp dictionary
        This function is intended to by called by PyRSIR.py

        sp: dictionary
        sp = {species1, species2, ...} //one key for every species to track
        sp[species] = {
        'bounds':[species x start, species x end], //start and end x values to integrate between
        'affin':['+' or '-' or 'UV'}, //which spectrum to look for this species in
        'level':integer, //if applicable, the MSn level (optional, but adds specificity)
        'function':integer, //the specific function in which to find this species (optional; overrides affin and level)
        }

        sumspec: bool
            toggles summing of all spectra together (creates an additional output item)
            also sums the spectra of mass spectrum species to generate an isotope pattern used by the bounds

        output:
            filled dictionary, each subkey will have:
            'raw': list of raw integrated values dictacted by the bounds
            'function': the function that the species was associated with

            if sumspec is true, will also output a dictionary of Spectrum objects
            the keys of this dictionary are the function numbers

        explicitly interprets full scan mass spectra and UV species
        """
        warnings.warn(
            'This is a legacy method used with the older PyRSIR, please use the new RSIRTarget and and PyRSIR class.',
            DeprecationWarning,
            stacklevel=2,
        )
        if sumspec is True:
            spec = {}
            for function in self.functions:  # create spectrum objects for all MS species
                if self.functions[function]['type'] == 'MS':
                    spec[function] = Spectrum(3)
        for species in sp:  # look for and assign function affinity
            sp[species]['function'] = self.associate_to_function(
                dct=sp[species])  # associate each species in the spectrum with a function
            if 'raw' not in sp[species]:  # look for empty raw list
                sp[species]['raw'] = []
        if self.ftt is False:  # if timepoints and tic values have not been extracted yet, extract those
            self.function_timetic()

        if self.verbose is True:
            prog = self.Progress(  # generate progress instance
                string='Extracting species data from spectrum',
                last=self.nscans,
                writeevery=5
            )
        for spectrum in self._spectra_elements:
            function, proc, scan = fps(spectrum)  # pull function, process, and scan numbers
            attr = branch_attributes(spectrum)  # get attributes
            if self.verbose is True:
                prog.write(attr['index'] + 1)  # outtput progress
                # self.sys.stdout.write('\rExtracting species data from spectrum #%d/%d  %.1f%%' %(attr['index']+1,self.nscans,float(attr['index']+1)/float(self.nscans)*100.))
            x, y = extract_spectrum(spectrum)  # generate spectrum
            if sumspec is True and function == 1:
                spec[function].add_spectrum(x, y)
            for key in sp:  # integrate each peak
                if sp[key]['function'] == function:  # if species is related to this function
                    if self.functions[function]['type'] == 'MS':
                        sp[key]['raw'].append(
                            self.integrate(key, sp[key]['bounds'][0], sp[key]['bounds'][1], x, y))  # integrate
                    if self.functions[function]['type'] == 'UV':
                        sp[key]['raw'].append(self.integrate(key, sp[key]['bounds'][0], sp[key]['bounds'][1], x,
                                                             y) / 1000000.)  # integrates and divides by 1 million bring it into au
        if self.verbose is True:
            prog.fin()  # write done
            # self.sys.stdout.write(' DONE\n')
        self._BE.printwarns()  # print bounds warnings (if any)
        if sumspec is True:
            return sp, spec
        return sp, None

    @property
    def scans(self):
        """a generator for scans in """
        return 

    def retrieve_scans(self, start=None, end=None, mzstart=None, mzend=None, function=None, mute=False, outside=False):
        """
        Retrieves the specified scans or time range from the specified function

        start: integer or float
            the point to start retrieving scans
            if integer, this will be a start scan number
            if float, this will be the start time
        end: (optional) integer or float
            the end point to stop retrieving scans
            same options as start
        mzstart: (optional) integer or float
            left m/z bound
        mzend: (optional) integer or float
            right m/z bound
        fn: integer
            the function to pull scans from (default 1)
        mute: bool
            overrides the verbose setting of the mzml instance
        outside: bool
            Whether to include the next point outside of the specified m/z bounds.
            This is useful for line continuity if the spectrum is to be used for
            rendering images.

        returns a list with each index corresponding to a scan, with two sublists for x and y data
        """
        if function is None:  # if not specified, retrieve first function
            function = self.associate_to_function()
        # find spectrum indicies to extract between
        if function not in self.functions:
            raise ValueError('The function "%d" is not in this mzml file.' % function)
        start = self.scan_index(start, function, bias='greater')
        end = self.scan_index(end, function, bias='lesser')
        if self.ftt is False:  # extract the timepoints and etc from the mzml
            self.function_timetic()
        msg = 'retrieving scans'
        out = []
        for spectrum in tqdm(self._spectra_elements, desc=msg, disable=self._mute_tqdm, total=end-start):
            attr = branch_attributes(spectrum)
            # func,proc,scan = self.fps(spectrum) # determine function, process, and scan numbers
            p = CVParameterSet.create_from_branch(spectrum)
            if attr['index'] > end:
                break
            if start <= attr['index'] <= end:  # within the index bounds
                x, y = extract_spectrum(spectrum)
                if mzstart is not None or mzend is not None:
                    if mzstart is None:
                        l = min(x)
                    else:
                        l = mzstart
                    if mzend is None:
                        r = max(x)
                    else:
                        r = mzend
                    spec = trimspectrum(x, y, l, r, outside)
                    out.append(spec)
                else:
                    out.append([x, y])
        if len(out) == 1:  # if only one scan, return that scan
            return out[0]
        return out

    def scan_index(self, scan=None, function=1, bias='lesser'):
        """
        Determines the index for a scan or timepoint in a given function

        :param int, float scan: The scan number (int) or time point (float) to find.
        :param int function: The mzml function to look in
        :param str bias: Bias of index finding (options dictacted by locate_in_list() )
        :return: scan index
        :rtype: int
        """
        if function not in self.functions:
            raise KeyError('The function %d is not in this mzML file.' % function)
        if scan is None:  # if no scan number is specified
            if bias == 'greater':  # used for start point
                return self.functions[function]['sr'][0]
            if bias == 'lesser':  # used for end point
                return self.functions[function]['sr'][1]
        if type(scan) is float:  # timepoint
            if self.ftt is False:
                self.function_timetic()
            # return located index plus start of the scan range
            return locate_in_list(self.functions[function]['timepoints'], scan, bias=bias) + self.functions[function]['sr'][0]
        elif type(scan) is int:  # scan number
            if scan < 1:
                raise ValueError('The scan number must be greater or equal to 1 (specified: %d)' % scan)
            if scan > self.functions[function]['nscans']:
                raise ValueError(f'The scan number {scan} exceeds the number of scans in function {function} '
                                 f'({self.functions[function]["nscans"]})')
            # return scan minus 1 (to shift into index domain) plus the start location index
            return scan - 1 + self.functions[function]['sr'][0]
        else:
            raise ValueError(f'An unexpected scan type was handed to the scan_index function ("{scan}", '
                             f'type: {type(scan)})')

    def sum_scans(self,
                  start=None,
                  end=None,
                  function=None,
                  dec=3,
                  mute=False
                  ):
        """
        Sums the specified scans together. If the scan range moves into another function, an error is raised.
        This method has a lower memory overhead than retrieve_scans().

        :param float, int start: start point to begin summing. ``int`` is interpreted as a scan number, ``float`` is
            interpreted as a time point in the acquisition.
        :param float, int end: end point to finish summing. Parameters are the same as with start.
        :param int function: mzML function to sum. If this is not provided, the first function will be used.
        :param int dec: number of decimal places to track in the spectrum (lower values lower memory overhead).
        :param bool mute: override chatty mode of mzML object
        :return: summed spectrum in the format ``[[m/z values], [intensity values]]``
        :rtype: list
        """

        # if no function is specified, use the first function
        if function is None:
            if len(self.functions) == 0:
                raise IndexError('The sum_scans method requires functions to be associated with the mzML file. There '
                                 'are none associated with this file. ')
            function = min(self.functions.keys())
        elif function not in self.functions:  # if fn is not defined
            raise KeyError(f'The function {function} is not defined in the mzML object. Available options: '
                           f'{", ".join([str(key) for key in self.functions.keys()])}')
        if self.functions[function]['type'] != 'MS':
            raise ValueError(f'The sum_scans function does not have the functionality to sum non-mass spec scans.'
                             f'The specified function {function} is of type {self.functions[function]["type"]}')
        start = self.scan_index(start, function, 'greater')
        end = self.scan_index(end, function, 'lesser')

        spec = Spectrum(  # create Spectrum object
            dec,
            start=self.functions[function]['window'][0],
            end=self.functions[function]['window'][1]
        )

        msg = 'combining spectra'
        for spectrum in tqdm(self._spectra_elements, desc=msg, disable=self._mute_tqdm or mute, total=end-start):
            attr = branch_attributes(spectrum)  # get attributes
            if attr['index'] > end:
                break
            if start <= attr['index'] <= end:  # if within the specified bounds
                x, y = extract_spectrum(spectrum)  # pull spectrum
                spec.add_spectrum(x, y)  # add spectrum to Spectrum object
        out = spec.trim()
        return out


if __name__ == '__main__':
    filename = 'MultiTest'
    mzml = mzML(filename, verbose=True, ftt=True)
    # sp = {
    # 'pos':{'bounds':[325,327],'affin':'+','spectrum':Spectrum(3),'raw':[]},
    # 'neg':{'bounds':[348,350],'affin':'-','spectrum':Spectrum(3),'raw':[]},
    # 'uv':{'bounds':[378,None],'affin':'UV','raw':[]}
    # }
