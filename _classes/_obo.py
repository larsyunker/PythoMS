"""
this class is already incorporated into the mzml class
it is retained in a separate file for convenience when searching for an accession key or name of an accession key
"""


class obo(object):
    """
    locates *.obo files and converts the most recent version into a python dictionary for parsing
    by default the class will look in the script's directory and the current working directory
    a path to a directory or to an obo file can also be suppliedS
    """

    def __init__(self, lookin=None):
        self.paths = []
        self.os = __import__('os')
        if lookin is not None:  # if a path is supplied
            if self.os.path.isdir(lookin) or self.os.path.isfile(lookin):
                self.paths.append(lookin)
            elif lookin.lower().endswith('.obo'):
                if self.os.path.isfile(lookin):
                    self.paths.append(lookin)
            else:
                raise ValueError('The supplied path to look in "%s" does not exist' % lookin)
        self.paths.append(self.os.path.realpath(__file__))  # the script's directory
        self.paths.append(self.os.getcwd())  # current working directory
        self.obo_path = self.find_obos()  # find obo files in the directory
        self.obo_path, self.ver = self.newest_version()  # refine to newest version obo file
        self.obodict = self.obo_to_dict(self.obo_path)

    def __str__(self):
        return '*.obo file version: %s path: %s' % (str(self.ver), self.obo_path)

    def __repr__(self):
        return 'obo(%s)' % self.obo_path

    def __getitem__(self, key):
        """
        returns the name of the accession key or the accession key that goes with the supplied name that is provided
        more extensive information about the accession key can be obtained by calling obo.print_properties(key)
        """
        if isinstance(key, slice):
            raise ValueError('Slicing of the obo object is not supported')
        if type(key) != str:
            raise ValueError('Only accession keys may be retrieved from the obo object')
        try:
            return self.obodict[key]['name']  # try to find name of accession key
        except KeyError:
            try:
                return self.acc_from_name(key)  # try to find accession key from name
            except KeyError:
                raise KeyError(
                    'They item "%s" is not defined as an accession key or a name of an accession key in the obo file.' % key)

    def acc_from_name(self, name):
        """finds the accession key from the name supplied"""
        for acc in self.obodict:
            if self.obodict[acc]['name'] == name:
                return acc
        raise KeyError('No accession key was found matching the supplied name "%s"' % name)

    def find_obos(self):
        """locates any *.obo files in the specified directory"""
        locations = []
        for eachpath in self.paths:
            for root, dirs, files in self.os.walk(eachpath):
                for ind in files:
                    if ind.endswith('.obo'):
                        locations.append(self.os.path.join(root, ind))
        if len(locations) == 0:
            raise IOError(
                'An *.obo file could not be located in the directory\n'
                '%s\nThis file can be obtained from the GitHub psi-ms-CV repository:\n'
                'https://github.com/HUPO-PSI/psi-ms-CV/blob/master/psi-ms.obo')
        return locations

    def newest_version(self):
        """of the supplied obo file paths, selects the most recent version"""
        ver = 0.
        out = ''
        for loc in self.obo_path:
            hndl = open(loc, 'rt')
            lines = hndl.readlines()
            hndl.close()
            for line in lines:
                if line.startswith('format-version:'):
                    if float(line.split()[1]) > ver:
                        ver = float(line.split()[1])
                        out = loc
                    break
        return out, ver

    def obo_to_dict(self, filepath):
        """converts the *.obo file into a dictionary"""
        hndl = open(filepath, 'rt')
        lines = hndl.readlines()
        hndl.close()

        out = {}
        for ind, line in enumerate(lines):
            if line.startswith('[Term]'):
                offset = 1  # offset initial state
                key = lines[ind + offset][4:-1]  # accession ID
                out[key] = {}
                while lines[ind + offset].startswith('[Term]') is False:
                    offset += 1
                    if ind + offset == len(lines):
                        break
                    if lines[ind + offset].startswith('name:'):  # accession name
                        out[key]['name'] = lines[ind + offset][6:-1]
                    elif lines[ind + offset].startswith('def:'):  # definition
                        out[key]['def'] = lines[ind + offset].split('"')[1]
                    elif lines[ind + offset].startswith(
                            'xref:'):  # xref (the coder cannot see a particular use for this key
                        out[key]['xref'] = lines[ind + offset].split()[1]
                    elif lines[ind + offset].startswith('is_a:'):  # is a something
                        if 'is_a' not in out[key]:
                            out[key]['is_a'] = []
                        out[key]['is_a'].append(lines[ind + offset].split()[1])
                        if lines[ind + offset].split()[1] not in out:  # catches undefined things like UO:0000000 ! unit
                            out[lines[ind + offset].split()[1]] = {'name': lines[ind + offset].split('!')[1][1:-1]}
                    elif lines[ind + offset].startswith('relationship:'):  # relationship keys
                        if 'relationship' not in out[key]:
                            out[key]['relationship'] = {}
                        spl = lines[ind + offset].split()
                        if spl[1] not in out[key]['relationship']:
                            out[key]['relationship'][spl[1]] = []
                        out[key]['relationship'][spl[1]].append(spl[2])
                    elif lines[ind + offset].startswith('synonym:'):  # synonym
                        if 'synonym' not in out[key]:
                            out[key]['synonym'] = []
                        out[key]['synonym'].append(lines[ind + offset].split()[1].split('"')[1])
                    elif lines[ind + offset].startswith('is_obsolete:'):  # obsolete
                        out[key]['obsolete'] = True
                    # the above keys seemed most pertinent, but additional keys can be coded here
        return out

    def print_properties(self, key):
        """prints detailed properties of the specified accession key"""
        if key not in self.obodict:
            raise KeyError('The accession "%s" is not defined in the obo file' % key)
        import sys
        sys.stdout.write('Accession: %s\n' % key)
        sys.stdout.write('Name: %s\n' % self.obodict[key]['name'])
        if 'obsolete' in self.obodict[key]:
            sys.stdout.write('This accession is obsolete.\n')
        sys.stdout.write('Definition: %s\n' % self.obodict[key]['def'])
        if 'is_a' in self.obodict[key]:
            for item in self.obodict[key]['is_a']:
                sys.stdout.write('Is a: %s (%s)\n' % (item, self.obodict[item]['name']))
        if 'relationship' in self.obodict[key]:
            sys.stdout.write('Relationships:\n')
            for subkey in self.obodict[key]['relationship']:
                for item in self.obodict[key]['relationship'][subkey]:
                    sys.stdout.write('%s: %s\n' % (subkey, item))
        if 'synonym' in self.obodict[key]:
            sys.stdout.write('Synonyms:\n')
            for item in self.obodict[key]['synonym']:
                sys.stdout.write('%s\n' % item)


if __name__ == '__main__':
    loadedobo = obo('C:\\1Files\\GitHub Repos\\psi-ms-CV')
    print
    loadedobo
