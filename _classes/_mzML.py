"""
Class for interpreting and processing mzML files

CHANGELOG:
    
    ---2.5 building

to add:
    removed obsolete functions in a couple of months (once all calls to them have been found)
    create attributes class (similar to cvparam) which will interpret the id string to as detailed a degree as possible
    add functionality to apply a calibration
    add plotspectrum call (which would use the tome_v02 plotter)
    identify spectrometer in softwareList
    try to extract timepoints and tic from chromatogramList (x values are sorted, so this probably won't work)
    obtain example files from different manufacturers and validate the script with those
    ask proteowizard support if it is possible to extract source conditions
"""

class mzML(object):
    def __init__(self,filename,**kwargs):
        """interprets and extracts information from a mzML (mass spectrum) file"""
        # check and set kewyord arguments
        self.ks = { # default keyword arguments
        'verbose': True, # toggle verbose
        'precision': 64, # floating point precision for values (32 or 64)
        'compression':True, # compression of binary strings (can substantially reduce file sizes)
        'gzip': True, # toggle gzip compression of mzml file (reduces file sizes even further
        'obo': None, # specific path to an *.obo file or the directory containing one
        }
        if set(kwargs.keys()) - set(self.ks.keys()): # check for invalid keyword arguments
            string = ''
            for i in set(kwargs.keys()) - set(self.ks.keys()):
                string += ` i`
            raise KeyError('Unsupported keyword argument(s): %s' %string)
        self.ks.update(kwargs) # update defaules with provided keyword arguments
        
        # import required functions
        self.sys = __import__('sys')
        self.os = __import__('os')
        self.filename = self.checkforfile(filename)
        self.b64 = __import__('base64')
        self.st = __import__('struct')
        self.zlib = __import__('zlib')
        from bisect import bisect_left as bl
        self.bl = bl # for convenience of calls
        self.sys.path.append(self.os.path.dirname(self.os.path.realpath(__file__))) # required so that this class can access other classes in the same directory
        
        # load file and determine key properties
        if self.ks['verbose'] is True:
            self.sys.stdout.write('Loading %s into memory' %self.filename)
            self.sys.stdout.flush()
        if self.filename.lower().endswith('.mzml.gz'): # if mzml is gzipped
            import gzip
            handle = gzip.open(self.filename) # unzip the file
        else:
            handle = self.filename
        import xml.dom.minidom
        try:
            self.tree = xml.dom.minidom.parse(handle) # full mzML file
        except:
            raise IOError('The mzML file "%s" could not be loaded. The file is either unsupported, corrupt, or incomplete.' %self.filename)
        self.mzmlcontents() # extract the contents of the mzML
        self.ftt = False # signifier for whether or not function_timetic has been run
        if self.ks['verbose'] is True:
            self.sys.stdout.write(' DONE\n')
        
    def __str__(self):
        """The string that is returned when printed"""
        return 'mzML file: "%s"\n# spectra: %d\n# chromatograms: %d' %(self.filename, self.nscans, self.nchroms)
    
    def __repr__(self):
        """The representation that is returned"""
        return "%s('%s')" %(self.__class__.__name__,self.filename)
    
    def __getitem__(self,ind):
        """
        returns the summed scans with the supplied scan number
        if a slice is provided, the script will attempt to sum those scans
        """
        if isinstance(ind,slice): # if getitem is trying to slice
            if ind.start is None: # no start
                start = 0
            else:
                start = ind.start
            if ind.stop is None: # no stop
                stop = self.nscans
            else:
                stop = ind.stop
        elif type(ind) is int:
            if ind == 0 or ind > self.nscans:
                raise IndexError("The scan #%d is outside of the mzML's scan range (1-%d)" %(ind,self.nscans))
            start = ind
            stop = ind
            if start < 0: # if negative indicies are supplied
                start = self.nscans+start+1
            if stop < 0:
                stop = self.nscans+stop+1
            return self.sumscans(sr=[start,stop],mute=True)[0]
        elif type(ind) is float:
            raise IndexError('Indexing by floating point (time point) is currently unsupported')
    
    def __add__(self,x):
        return 'Addition to the mzML class is unsupported'
    def __sub__(self,x):
        return 'Subtraction from the mzML class is unsupported'
    def __mul__(self,x):
        return 'Multiplication of the mzML class is unsupported'
    def __div__(self,x):
        return 'Division of the mzML class is unsupported'
    
    class BoundsError(Warning):
        """A warning class to handle bounds errors when integrating"""
        def __init__(self):
            self.warned = {}
        def printwarns(self):
            """prints the number of warnings if merited"""
            if len(self.warned) > 0:
                import sys
                sys.stdout.write('The following peaks exceeded the bounds of the spectrum n number of times:\n(number of scans in the file: %d)\n'%self.nscans)
                for name in self.warned:
                    sys.stdout.write('"%s": %d\n' %(name,self.warned[name]))
        def warn(self,name,intstart,intend,mzstart,mzend):
            """warns the user if there was a mismatch"""
            if self.warned.has_key(name) is False:
                import sys
                sys.stdout.write('\nThe peak "%s" (%f-%f) is outside of the bounds of the spectrum being summed m/z %.1f-%.1f\n' %(name,intstart,intend,mzstart,mzend))
                self.warned[name] = 1
            else:
                self.warned[name] += 1

    class cvparam(object):
        """interprets all cvparameter tags within the provided branch and can provide details as required with flexible input"""
        def __init__(self,branch):
            self.parameters = self.interpret(branch)
            self.psorted = sorted(self.parameters) # sorted key list for item getting
        
        def __getitem__(self,key):
            """
            returns the value of the accesssion key supplied
            (can also be handed an item name, provided that it is defined in this cvParameter set)
            """
            if isinstance(key,slice):
                raise ValueError('Slicing of a cvparam is not supported')
            if type(key) is int:
                return self.psorted[key] # return sorted parameter index
            #if type(key) != str:
            #    raise ValueError('Only accession keys or names may be retrieved from the cvparam object (retrieval attempt: "%s")' %str(key))
            if self.parameters.has_key(key) is True: # if the key matches an accession key
                return self.returnvalue(key)
            try:
                return self.parameters[self.acc_by_name(key)]['value'] # if a name was provided and it matches an accession key
            except KeyError:
                raise KeyError('The accession key "%s" is not present in this cvParameters set')
            
        def acc_by_name(self,name):
            """trys to find an accession by a provided name"""
            for key in self.parameters:
                if self.parameters[key]['name'] == name:
                    return key
            raise KeyError('An accession key matching the name "%s" is not defined in this cvParameters set' %name)
        
        def has_key(self,key):
            """checks for the presence of an accession key in the cvparameters set"""
            return self.parameters.has_key(key)
        
        def interpret(self,branch):
            """
            retrieves all the cvParam tags for a given branch
            returns a dictionary with keys corresponding to accession number
            each key is a subdictionary with the attributes of the cvParam
            """
            def stringtodigit(string):
                """attempts to convert a unicode string to float or integer"""
                try:
                    value = int(string) # try converting to integer
                except ValueError:
                    try:
                        value = float(string) # try converting to float
                    except ValueError:
                        value = string # otherwise keep as unicode
                return value
            out = {}
            for cvParam in branch.getElementsByTagName('cvParam'):
                acc = cvParam.getAttribute('accession') # accession key
                out[acc] = {}
                for attribute,value in cvParam.attributes.items(): # pull all the attributes
                    if attribute != 'accession':
                        out[acc][attribute] = stringtodigit(value) # attempt to convert to integer or float, keep as string otherwise
            return out
        
        def returnvalue(self,key):
            """returns the value of a particular accession key in the parameters dictionary"""
            try:
                return self.parameters[key]['value']
            except KeyError:
                raise KeyError('The accession key "%s" is not present in this cvParameters set')
        
        def unitname(self):
            """
            looks for units in the cvParameters and returns the appropriate unit name
            all units as of 2016-07-12 are defined here, but all possible units can be found in
            https://raw.githubusercontent.com/HUPO-PSI/psi-ms-CV/master/psi-ms.obo
            """
            unitkeys = {
            # x units
            'MS:1000040':'m/z',
            'UO:0000010':'second',
            'UO:0000017':'micrometer',
            'UO:0000018':'nanometer',
            'UO:0000028':'millisecond',
            'UO:0000031':'minute',
            'UO:0000221':'dalton',
            'UO:0000222':'kilodalton',
            # y units
            'MS:1000131':'number of detector counts',
            'MS:1000132':'percent of base peak',
            'MS:1000814':'counts per second',
            'MS:1000905':'percent of base peak times 100',
            'UO:0000187':'percent',
            'UO:0000269':'absorbance unit',
            # other
            'MS:1000807':'Th/s',
            'UO:0000002':'mass unit',
            'UO:0000008':'meter',
            'UO:0000012':'kelvin',
            'UO:0000021':'gram',
            'UO:0000027':'degree Celsius',
            'UO:0000098':'milliliter',
            'UO:0000106':'hertz',
            'UO:0000110':'pascal',
            'UO:0000112':'joule',
            'UO:0000166':'parts per notation unit',
            'UO:0000169':'parts per million',
            'UO:0000175':'gram per liter',
            'UO:0000185':'degree',
            'UO:0000218':'volt',
            'UO:0000228':'tesla',
            'UO:0000266':'electronvolt',
            'UO:0000268':'volt per meter',
            }
            for key in self.parameters:
                if self.parameters[key].has_key('unitAccession'): # find unitAccession (the accession code defines the unit)
                    if self.parameters[key]['unitAccession'] in unitkeys: # if it is defined, return unit code
                        return unitkeys[self.parameters[key]['unitAccession']]
                    else:
                        if mzML.__dict__.has_key('loadedobo') is False:
                            mzML.loadedobo = mzML.obo(mzML.ks['obo']) # load obo file
                        return mzML.loadedobo[self.parameters[key]['unitAccession']] # return accession key name as defined in obo file
        
    class obo(object):
        """
        locates *.obo files and converts the most recent version into a python dictionary for parsing
        by default the class will look in the script's directory and the current working directory
        a path to a directory or to an obo file can also be suppliedS
        """
        def __init__(self,lookin=None):
            self.paths = []
            self.os = __import__('os')
            if lookin is not None: # if a path is supplied
                if self.os.path.isdir(lookin) or self.os.path.isfile(lookin):
                    self.paths.append(lookin)
                elif lookin.lower().endswith('.obo'):
                    if self.os.path.isfile(lookin):
                        self.paths.append(lookin)
                else:
                    raise ValueError('The supplied path to look in "%s" does not exist' %lookin)
            self.paths.append(self.os.path.realpath(__file__)) # the script's directory
            self.paths.append(self.os.getcwd()) # current working directory
            self.obo_path = self.find_obos() # find obo files in the directory
            self.obo_path,self.ver = self.newest_version() # refine to newest version obo file
            self.obodict = self.obo_to_dict(self.obo_path)
            
        def __str__(self):
            return '*.obo file version: %s path: %s' %(str(self.ver),self.obo_path)
        def __repr__(self):
            return 'obo(%s)'%self.obo_path
        
        def __getitem__(self,key):
            """
            returns the name of the accession key or the accession key that goes with the supplied name that is provided
            more extensive information about an accession key can be obtained by calling obo.print_properties(key)
            """
            if isinstance(key,slice):
                raise ValueError('Slicing of the obo object is not supported')
            if type(key) != str:
                raise ValueError('Only accession keys may be retrieved from the obo object')
            try:
                return self.obodict[key]['name'] # try to find name of accession key
            except KeyError:
                try:
                    return self.acc_from_name(key) # try to find accession key from name
                except KeyError:
                    raise KeyError('They item "%s" is not defined as an accession key or a name of an accession key in the obo file.' %key)
        
        def acc_from_name(self,name):
            """finds the accession key from the name supplied"""
            for acc in self.obodict:
                if self.obodict[acc]['name'] == name:
                    return acc
            raise KeyError('No accession key was found matching the supplied name "%s"' %name) 
                
        def find_obos(self):
            """locates any *.obo files in the specified directory"""
            locations = []
            for eachpath in self.paths:
                for root,dirs,files in self.os.walk(eachpath):
                    for ind in files:
                        if ind.endswith('.obo'):
                            locations.append(self.os.path.join(root,ind))
            if len(locations) == 0:
                raise IOError('An *.obo file could not be located in the directory\n%s\nThis file can be obtained from the GitHub psi-ms-CV repository:\nhttps://github.com/HUPO-PSI/psi-ms-CV/blob/master/psi-ms.obo')
            return locations
        
        def newest_version(self):
            """of the supplied obo file paths, selects the most recent version"""
            ver = 0.
            out = ''
            for loc in self.obo_path:
                hndl = open(loc,'rt')
                lines = hndl.readlines()
                hndl.close()
                for line in lines:
                    if line.startswith('format-version:'):
                        if float(line.split()[1]) > ver:
                            ver = float(line.split()[1])
                            out = loc
                        break
            return out,ver
            
        def obo_to_dict(self,filepath):
            """converts the *.obo file into a dictionary"""
            hndl = open(filepath,'rt')
            lines = hndl.readlines()
            hndl.close()
            
            out = {}
            for ind,line in enumerate(lines):
                if line.startswith('[Term]'):
                    offset = 1 # offset initial state
                    key = lines[ind+offset][4:-1] # accession ID
                    out[key] = {}
                    while lines[ind+offset].startswith('[Term]') is False:
                        offset += 1
                        if ind+offset == len(lines):
                            break
                        if lines[ind+offset].startswith('name:'): # accession name
                            out[key]['name'] = lines[ind+offset][6:-1]
                        elif lines[ind+offset].startswith('def:'): # definition
                            out[key]['def'] = lines[ind+offset].split('"')[1]
                        elif lines[ind+offset].startswith('xref:'): # xref (the coder cannot see a particular use for this key
                            out[key]['xref'] = lines[ind+offset].split('"')[1]
                        elif lines[ind+offset].startswith('is_a:'): # is a something
                            if out[key].has_key('is_a') is False:
                                out[key]['is_a'] = []
                            out[key]['is_a'].append(lines[ind+offset].split()[1])
                            if out.has_key(lines[ind+offset].split()[1]) is False: # catches undefined things like UO:0000000 ! unit
                                out[lines[ind+offset].split()[1]] = {'name':lines[ind+offset].split('!')[1][1:-1]}
                        elif lines[ind+offset].startswith('relationship:'): # relationship keys
                            if out[key].has_key('relationship') is False:
                                out[key]['relationship'] = {}
                            spl = lines[ind+offset].split()
                            if out[key]['relationship'].has_key(spl[1]) is False:
                                out[key]['relationship'][spl[1]] = []
                            out[key]['relationship'][spl[1]].append(spl[2])
                        elif lines[ind+offset].startswith('synonym:'): # synonym
                            if out[key].has_key('synonym') is False:
                                out[key]['synonym'] = []
                            out[key]['synonym'].append(lines[ind+offset].split()[1].split('"')[1])
                        elif lines[ind+offset].startswith('is_obsolete:'): # obsolete
                            out[key]['obsolete'] = True
                        # the above keys seemed most pertinent, but additional keys can be coded here
            return out
        
        def print_properties(self,key):
            """prints detailed properties of the specified accession key"""
            if self.obodict.has_key(key) is False:
                raise KeyError('The accession "%s" is not defined in the obo file' %key)
            import sys
            sys.stdout.write('Accession: %s\n' %key)
            sys.stdout.write('Name: %s\n' %self.obodict[key]['name'])
            if self.obodict[key].has_key('obsolete'):
                sys.stdout.write('This accession is obsolete.\n')
            sys.stdout.write('Definition: %s\n' %self.obodict[key]['def'])
            if self.obodict[key].has_key('is_a') is True:
                for item in self.obodict[key]['is_a']:
                    sys.stdout.write('Is a: %s (%s)\n' %(item,self.obodict[item]['name']))
            if self.obodict[key].has_key('relationship') is True:
                sys.stdout.write('Relationships:\n')
                for subkey in self.obodict[key]['relationship']:
                    for item in self.obodict[key]['relationship'][subkey]:
                        sys.stdout.write('%s: %s\n' %(subkey,item))
            if self.obodict[key].has_key('synonym') is True:
                sys.stdout.write('Synonyms:\n')
                for item in self.obodict[key]['synonym']:
                    sys.stdout.write('%s\n' %item)    
    
    def _foreachchrom(self,fn):
        """
        a decorator function that will apply the supplied function to every chromatogram in the mzml file
        the supplied function will be handed the chromatogram XML object as the first argument
        the decorated function will return a list of outputs of the supplied function where each index corresponds to a scan
        """
        def foreachchrom(*args,**kwargs):
            """decorates the supplied function to run for every scan"""
            out = []
            for chromatogram in self.tree.getElementsByTagName('chromatogram'):
                current = int(chromatogram.getAttribute('index'))+1
                if self.ks['verbose'] is True:
                    self.sys.stdout.write('\rApplying function "%s" to chromatogram #%d/%d %.1f%%' %(fn.__name__,current,self.chroms,float(current)/float(self.nchroms)*100.))
                out.append(fn(chromatogram,*args,**kwargs))
            if self.ks['verbose'] is True:
                self.sys.stdout.write(' DONE\n')
            return out
        return foreachchrom
    
    def _foreachscan(self,fn):
        """
        a decorator function that will apply the supplied function to every spectrum in the mzml file
        the supplied function will be handed the spectrum XML object as the first argument
        the decorated function will return a list of outputs of the supplied function where each index corresponds to a scan
        """
        def foreachscan(*args,**kwargs):
            """decorates the supplied function to run for every scan"""
            out = []
            for spectrum in self.tree.getElementsByTagName('spectrum'):
                current = int(spectrum.getAttribute('index'))+1
                if self.ks['verbose'] is True:
                    self.sys.stdout.write('\rApplying function "%s" to scan #%d/%d %.1f%%' %(fn.__name__,current,self.nscans,float(current)/float(self.nscans)*100.))
                out.append(fn(spectrum,*args,**kwargs))
            if self.ks['verbose'] is True:
                self.sys.stdout.write(' DONE\n')
            return out
        return foreachscan        
    
    def associate_to_function(self,affin=None,level=None,dct=None):
        """
        associates an affinity and/or level with a function in the mzML instance
        affin: '+' '-' or 'UV'
            MS mode or UV-Vis
        level: integer
            MSn level (only applicable if a MS spectrum)
        dct: dictionary
            can also hand a dictionary, in which the script will look for affin and level keys
        """
        if dct is not None: # if function was handed a dictionary
            if dct.has_key('function'):
                return dct['function']
            if dct.has_key('affin'):
                affin = dct['affin']
            if dct.has_key('level'):
                level = dct['level']
        
        if affin is None and level is None:
            return 1 # assume first function
        
        elif affin == 'UV': # if UV-Vis affinity
            for fn in self.functions: # determine which function is UV-Vis
                if self.functions[fn]['acc'] == 'MS:1000804':
                    return fn
            raise ValueError('There is no electromagnetic radiation spectrum function in this mzML file')
        
        elif affin in ['+','-']: # if affinity to mass spectrum
            levelcount = 0 # counter for number of matches to this affinity and level
            for fn in self.functions:
                if self.functions[fn]['type'] == 'MS': # if fn is ms
                    if self.functions[fn]['mode'] == affin: # if mode mathes
                        if level is None and self.functions[fn]['level'] == 1: # if there is no level specified, assume 1
                            fnout = fn
                            levelcount += 1
                        elif self.functions[fn]['level'] == level: # if level matches
                            fnout = fn
                            levelcount += 1
            if levelcount > 1:
                raise ValueError("There affinity specification of mode: %s, level: '%d' matches more than one function in the mzML file.\nTo process this species, be more specific in your level specification or assign it to a specific function number by adding a 'function' key to its dictionary." %(affin,level))
            return fnout
        else: # if some other affinity
            raise ValueError('The specified affinity "%s" is not supported.' %affin)
                                                                                                                                                
    def attributes(self,branch):
        """pulls all attributes of a supplied branch and creates a dictionary of them"""
        def stringtodigit(string):
            """attempts to convert a unicode string to float or integer"""
            try:
                value = int(string) # try converting to integer
            except ValueError:
                try:
                    value = float(string) # try converting to float
                except ValueError:
                    value = string # otherwise keep as unicode
            return value
        out = {}
        for pair in branch.attributes.items():
            out[pair[0]] = stringtodigit(pair[1])
        return out
    
    def autoresolution(self,n=10,fn=1):
        """
        automatically determines the resolution of the spectrometer that recorded the mzml file
        resolution is based on the average resolution of 10 pseudo-random samples
        each sample spectrum is split into 4 sections and 4 peaks are found to calculate the resolution
        """
        def findsomepeaks(y):
            """roughly locates 4 peaks by maximum values in the spectrum and returns their index"""
            split = int(len(y)/4)
            start = 0
            end = start+split
            splity = []
            for i in range(4):
                splity.append(sci.asarray(y[start:end]))
                start += split
                end += split
            out = []
            for ind,section in enumerate(splity):
                maxy = max(section)
                if maxy == max(section[1:-1]): # if max is not at the edge of the spectrum
                    out.append(sci.where(section==maxy)[0][0]+split*ind)
            return out
            
        def resolution(x,y,index=None):
            """
            Finds the resolution and full width at half max of a spectrum
            x: list of mz values
            y: corresponding list of intensity values
            index: index of maximum intensity (optional; used if the resolution of a specific peak is desired)
            
            returns resolution
            """
            y = sci.asarray(y) # convert to array for efficiency
            if index is None: # find index and value of maximum
                maxy = max(y)
                index = sci.where(y==maxy)[0][0]
            else:
                maxy = y[index]
            if maxy/(sum(y)/len(y)) < 10: # if intensity to average is below this threshold (rough estimate of signal/noise)
                return None
            halfmax = maxy/2
            indleft = int(index)-1 # generate index counters for left and right walking
            indright = int(index)+1
            while y[indleft] > halfmax: # while intensity is still above halfmax
                indleft -= 1
            while y[indright] > halfmax:
                indright += 1
            return x[index]/(x[indright]-x[indleft]) # return resolution (mz over full width at half max)
        
        from random import random
        import scipy as sci
        if self.functions[fn]['type'] != 'MS':
            raise ValueError('The autoresolution function only operates on mass spectrum functions. Type of specified function %d: %s' %(fn,self.function[fn]['type']))
        ranges = [] # list of scan intervals
        
        while len(ranges) < n: # generate 10 pseudo-random intervals to sample
            ran = int(random()*self.functions[fn]['nscans']) + self.functions[fn]['sr'][0]
            if ran-10 >= self.functions[fn]['sr'][0] and ran+10 <= self.functions[fn]['sr'][1]:
                ranges.append([ran-10,ran+10])
        summed = []
        for ind,rng in enumerate(ranges):
            if self.ks['verbose'] is True:
                self.sys.stdout.write('\rEstimating resolution of the instrument %.0f%%' %(float(ind+1)/float(n)*100.))
            summed.append(self.sum_scans(rng[0],rng[1],fn,2,True)) # sum those scans and append output
        res = []
        for spec in summed: # calculate resolution for each scan range
            inds = findsomepeaks(spec[1]) # find some peaks
            for ind in inds: # for each of those peaks
                res.append(resolution(spec[0],spec[1],ind))
        if self.ks['verbose'] is True:
            self.sys.stdout.write(' DONE\n')
        res = [y for y in res if y is not None] # removes None values (below S/N)
        return sum(res)/len(res) # return average
        
    def checkforfile(self,fn):
        """checks for file and converts if necessary"""
        def version_input(string):
            """checks the python version and uses the appropriate version of user input"""
            import sys
            if sys.version.startswith('2.7'):
                return raw_input('%s' %string)
            if sys.version.startswith('3.'):
                return input('%s' %string)
            else:
                raise EnvironmentError('The version_input method encountered an unsupported version of python.')
        
        valid = ['.raw','.mzml.gz','.mzml'] # supported extensions
        if fn.lower().endswith('.raw') is True: # extension is raw
            if self.filepresent(fn[:-4]+'.mzML.gz') is True: # if corresponding gzipped mzml is present
                return fn[:-4]+'.mzML.gz'
            if self.filepresent(fn[:-4]+'.mzML') is True: # if corresponding mzml is present
                return fn[:-4]+'.mzML'
            return self.pwconvert(fn,self.ks['precision'],self.ks['compression'],self.ks['gzip']) # otherwise convert and return mzml
        elif self.filepresent(fn) is True: # if the specified file is present
            for exten in valid: # checks for supported extensions
                if fn.lower().endswith(exten) is True:
                    return fn
            # otherwise asks user whether to continue
            if version_input('The extension of the supplied filename "%s" is unexpected and may not be supported.\nDo you wish to proceed with file loading? [Y/N] ' %fn).lower() in ['y','yes']:
                return fn
            else:
                self.sys.exit('The user cancelled mzML loading.')
        else:
            fn = self.fixextension(fn) # try to fix extension
            if fn.lower().endswith('.raw') is True: # convert if only raw file is found
                return self.pwconvert(fn,self.ks['precision'],self.ks['compression'],self.ks['gzip'])
            return fn
    
    def extractspectrum(self,spectrum,units=False):
        """pulls and converts binary data to list"""
        def gettext(nodelist):
            """gets text from a simple XML object"""
            rc = []
            for node in nodelist:
                if node.nodeType == node.TEXT_NODE:
                    rc.append(node.data)
            return ''.join(rc)
        
        def decodeformat(p,speclen):
            """determines the decode format from the accession parameter"""
            formats = {
            'MS:1000519':['<','i'], # signed 32-bit little-endian integer
            #'MS:1000520':['',''], # [OBSOLETE] Signed 16-bit float
            'MS:1000521':['<','f'], # 32-bit precision little-endian floating point conforming to IEEE-754
            'MS:1000522':['<','l'], # Signed 64-bit little-endian integer
            'MS:1000523':['<','d'], # 64-bit precision little-endian floating point conforming to IEEE-754.
            }
            for key in formats:
                if p.has_key(key): # find accession number match
                    return formats[key][0]+str(speclen)+formats[key][1]
            
        speclen = int(spectrum.getAttribute('defaultArrayLength')) # spectrum length (defined in the spectrum attricubes)
        out = []
        if units is True:
            units = []
        for binary in spectrum.getElementsByTagName('binaryDataArray'):
            p = self.cvparam(binary) # grab cvparameters
            if p.has_key('MS:1000574') is True: # determine whether the binary string is zlib compressed
                compressed = True
            else:
                compressed = False
            unpack_format = decodeformat(p,speclen) # determine unpack format
            string = gettext(binary.getElementsByTagName('binary')[0].childNodes) # pull the binary string
            decoded = self.b64.decodestring(string) # decode the string
            if compressed is True: # if the string is compressed, decompress
                decoded = self.zlib.decompress(decoded)
            out.append(list(self.st.unpack(unpack_format,decoded))) # unpack the string
            if units is not False:
                units.append(p.unitname())
        if units is not False: # extends the units onto out
            out.extend(units)
        return out
    
    def filepresent(self,filepath):
        """checks for the presence of the specified file or directory in the current working directory"""
        tf = self.os.path.isfile(filepath) # look for file first
        if tf is False: # if file cannot be found, look for directory
            tf = self.os.path.isdir(filepath)
        return tf
        
    def fixextension(self,fn):
        """tries to fix invalid file extensions"""
        oopsx = {'.mzm':'l','.mz':'ml','.m':'zml','.':'mzml'} # incomplete mzml extensions
        oopsr = {'.ra':'w','.r':'aw','.':'raw'} # incomplete raw extionsions
        oopsg = {'.mzml.g':'z','.mzml.':'gz','.mzml':'.gz','.mzm':'l.gz','.mz':'ml.gz','.m':'zml.gz','.':'mzml.gz'} #incomplete gz extensions
        # looks for missing extensions first
        if self.filepresent(fn+'.mzml.gz') is True:
            return fn+'.mzml.gz'
        if self.filepresent(fn+'.mzml') is True: 
            return fn+'.mzml'
        for key in oopsg: # tries to complete mzml.gz shortenings
            if fn.lower().endswith(key) is True:
                if self.filepresent(fn+oopsg[key]) is True:
                    return fn+oopsg[key]
        for key in oopsx: # tries to complete mzml shortenings
            if fn.lower().endswith(key) is True:
                if self.filepresent(fn+oopsx[key]) is True:
                    return fn+oopsx[key]
        for key in oopsr: # tries to complete raw shortenings
            if fn.lower().endswith(key) is True:
                if self.filepresent(fn+oopsr[key]) is True:
                    return fn+oopsr[key]
        if self.filepresent(fn+'.raw') is True: # finally looks for raw file
            return fn+'.raw'
        raise IOError('The file "%s" could not be located in the current working directory'%(fn)) # if it can't be found, raise IOError
    
    def fps(self,branch):
        """
        extracts function #, process #, and scan # from the idstring of a spectrum branch
        returns function, process, scan as integers
        """
        idstring = branch.getAttribute('id').split() # pull id string from scan attribute
        return [int(x.split('=')[1]) for x in idstring] # return each value after converting to integer
    
    def function_timetic(self):
        """
        extracts timepoints and tic lists for each function
        this function is separate from mzml contents because it would increase load times significantly (~6x)
        """
        if self.ftt is True: # if this has already been done, break out
            return None
        for func in self.functions: # add timepoint and tic lists
            self.functions[func]['timepoints'] = [] # list for timepoints
            self.functions[func]['tic'] = [] # list for total ion current values
            if self.functions[func].has_key('level') and self.functions[func]['level'] > 1:
                self.functions[func]['ce'] = [] # list for collision energies
        for spectrum in self.tree.getElementsByTagName('spectrum'):
            attr = self.attributes(spectrum)
            func,proc,scan = self.fps(spectrum) # determine function, process, and scan numbers
            if self.ks['verbose'] is True:
                self.sys.stdout.write('\rExtracting timepoints and total ion current values from mzML %.1f%%' %(float(attr['index']+1)/float(self.nscans)*100.))
            p = self.cvparam(spectrum) # pull spectrum's cvparameters
            self.functions[func]['timepoints'].append(p['MS:1000016'])
            self.functions[func]['tic'].append(p['MS:1000285'])
            if p.has_key('MS:1000045'):
                self.functions[func]['ce'].append(p['MS:1000045'])
        if self.ks['verbose'] is True:
            self.sys.stdout.write(' DONE\n')
        self.ftt = True
            
    def integrate(self,name,start,end,x,y):
        """
        Integrates y values given x bounds in a paired set of lists (e.g. a m/z list and an intensity list)
        
        name: name of the peak being integrated (only used for warning purposes)
        start: start x value
        end: end x value
        x: list of x values
        y: list of y values (paired with x)
        
        returns: integral value
        """
        if start > max(x) or start < min(x): # check that start is within the m/z bounds
            self.BE.warn(name,start,end,min(x),max(x))
        if end is None: # if only a start value is supplied, return closest to that value
            return y[self.locateinlist(x,start)]
        if end > max(x): # check that end is within the m/z bounds
            self.BE.warn(name,start,end,min(x),max(x))
        return sum(y[self.locateinlist(x,start,'greater'):self.locateinlist(x,end,'lesser')]) # integrate using the nearest values inside the bounds        
    
    def locateinlist(self,lst,value,bias='closest'):
        """
        Finds index in a sorted list of the value closest to a given value
        
        If two numbers are equally close, return the smallest number.
        based on http://stackoverflow.com/questions/12141150/from-list-of-integers-get-number-closest-to-a-given-value
        
        lst: list
            list of values to search
        value: float or int
            value number to find
        bias: 'lesser','greater', or 'closest'
            default 'left'
            'lesser' will return the position of the value just less than the provided value
            'greater' will return the position of the value just greater than the provided value
            'closest' will return the index of the nearest value to the one provided
        """                
        pos = self.bl(lst, value)
        if pos == 0: # if at start of list
            return pos
        elif pos == len(lst): # if insertion is beyond index range
            return pos -1 
        if lst[pos] == value: # if an exact match is found
            return pos
        if bias == 'greater': # return value greater than the value (bisect_left has an inherent bias to the right)
            return pos
        if bias == 'lesser': # return value lesser than the provided
            return pos -1
        if bias == 'closest': # check differences between index and index-1 and actual value, return closest
            adjval = abs(lst[pos-1] - value)
            curval = abs(lst[pos] - value)
            if adjval < curval: # if the lesser value is closer
                return pos-1
            if adjval == curval: # if values are equidistant
                return pos-1
            else:
                return pos
    
    def mzmlcontents(self):
        """finds the total number of scans, the number of chromatograms, and the scan range for each function in the mzml file"""
        self.nscans = int(self.tree.getElementsByTagName('spectrumList')[0].getAttribute('count')) # number of spectra
        self.nchroms = int(self.tree.getElementsByTagName('chromatogramList')[0].getAttribute('count')) # number of chromatograms
        self.functions = {}
        for spectrum in self.tree.getElementsByTagName('spectrum'):
            func,proc,scan = self.fps(spectrum) # extract each value and convert to integer
            if func not in self.functions: # if function is not defined yet
                p = self.cvparam(spectrum) # pull spectrum's cvparameters
                self.functions[func] = {
                'sr':[int(spectrum.getAttribute('index')),None], # the scan index range that the function spans
                'nscans':1, # number of scans
                }
                self.functions[func].update(self.scanproperties(p)) # update with scan properties
            else:
                self.functions[func]['sr'][1] = int(spectrum.getAttribute('index')) # otherwise set the scan index range to the current index
                self.functions[func]['nscans'] += 1
    
    def pull_chromatograms(self):
        """
        Pulls mzML chromatograms
        returns:
            dictionary of all chromatograms with keys denoted by the chromatogram id
            each dictionary has subkeys for 'x':xlist, 'y':ylist, 'xunit':unit of x, 'yunit':unit of y
        """
        chroms = {} #dictionary of chromatograms
        for chromatogram in self.tree.getElementsByTagName('chromatogram'):
            attr = self.attributes(chromatogram) # pull attributes
            if self.ks['verbose'] is True:
                self.sys.stdout.write('\rExtracting chromatogram #%s/%i  %.1f%%' %(attr['index']+1,self.nchroms,float(attr['index']+1)/float(self.nchroms)*100.))
                self.sys.stdout.flush()
            x,y,xunit,yunit = self.extractspectrum(chromatogram,True) # extract x list, y list, and units
            chroms[attr['id']] = {'x':x, 'y':y, 'xunit':xunit, 'yunit':yunit}
        if self.ks['verbose'] is True:
            self.sys.stdout.write(' DONE\n')
        return chroms

    def pullmsmsspectra(self):
        """
        OBSOLETE (use retrieve_scans or sum_scans)
        Extracts MSMS spectra from the mzML file
        
        Exclusively finds mass spectra with level greater or equal to 2
        groups spectra by key defined by the isolation window target
        
        each dictionary has subkeys for each time point which has subsubkeys containing the information of that scan
        
        also returns a dictionary of m/z limits using the same keys as msms
        
        output format:
        msms = {function, function,...} one entry for every msms function
        msms[function] = [scan1,scan2,...] one entry for every scan in that function (ordered)
        msms[function][index] = {'CE':value,'x':list,'y':list}
        """
        self.sys.exit('The pullmsmsspectra function is obsolete. Identify the appropriate function to sum from mzML.functions and use retrieve_scans instead.')
    
    def pullscansfromfn(self,fn,scan=None,endscan=None):
        """
        OBSOLETE (use retrieve_scans or sum_scans)
        pulls specific scan(s) from the specified function number in the mzml
        fn: int
            function number (as defined in the mzML file)
        scan: int
            scan number to return (scan number within that function)
        endscan: int
            end scan of scan range (optional)
        
        returns: [xlist,ylist]
        """
        self.sys.exit('The pullscansfromfn function is obsolete. Identify the appropriate function to sum from mzML.functions and use retrieve_scans instead.')
    
    def pull_species_data(self,sp,sumspec=False):
        """
        Extracts integrated data at every timepoint for all species specified in the sp dictionary
        
        sp: dictionary
        sp = {species1, species2, ...} //one key for every species to track
        sp[species] = {
        'bounds':[species x start, species x end], //start and end x values to integrate between
        'affin':['+' or '-' or 'UV'}, //which spectrum to look for this species in
        'level':integer, //if applicable, the MSn level (optional, but add specificity)
        'function':integer, //the specific function in which to find this species (optional; overrides affin and level)
        'raw':[], //an empty list for the raw values to be inserted into
        'spectrum':Spectrum object //a spectrum object to contain the summed isotope pattern given by the bounds (only applicable for mass spectrum species)
        }
        
        sumspec: bool
            toggles summing of all spectra together (creates an additional output item)
            also sums the spectra of mass spectrum species to generate an isotope pattern used by the bounds
        
        output:
            filled dictionary with integrated values inserted into 'raw' key
            if sumspec is true, will also output an [x,y] list
        
        explicitly interprets full scan mass spectra and UV species
        """
        if sumspec is True:
            from _Spectrum import Spectrum
            spec = Spectrum(3)
        for species in sp: # look for and assign function affinity
            sp[species]['function'] = self.associate_to_function(dct=sp[species]) # associate each species in the spectrum with a function
            if sp[species].has_key('raw') is False: # look for empty raw list
                sp[species]['raw'] = []
        if self.ftt is False: # if timepoints and tic values have not been extracted yet, extract those
            self.function_timetic()
        self.BE = self.BoundsError() # load warning instance for integration
        for spectrum in self.tree.getElementsByTagName('spectrum'):
            func,proc,scan = self.fps(spectrum) # pull function, process, and scan numbers
            attr = self.attributes(spectrum) # get attributes
            if self.ks['verbose'] is True:
                self.sys.stdout.write('\rExtracting species data from spectrum #%d/%d  %.1f%%' %(attr['index']+1,self.nscans,float(attr['index']+1)/float(self.nscans)*100.))
            x,y = self.extractspectrum(spectrum) # generate spectrum
            if sumspec is True:
                spec.addspectrum(x,y)
            for key in sp: # integrate each peak
                if sp[key]['function'] == func: # if species is related to this function
                    if self.functions[func]['type'] == 'MS':
                        sp[key]['raw'].append(self.integrate(key,sp[key]['bounds'][0],sp[key]['bounds'][1],x,y)) # integrate
                    if self.functions[func]['type'] == 'UV':
                        sp[key]['raw'].append(self.integrate(key,sp[key]['bounds'][0],sp[key]['bounds'][1],x,y)/1000000.) # integrates and divides by 1 million bring it into au
        if self.ks['verbose'] is True:
            self.sys.stdout.write(' DONE\n')
        self.BE.printwarns() # print bounds warnings (if any)
        if sumspec is True:
            return sp,spec  
        else:
            return sp
    
    def pullspectra(self,sr=None,tr=None,mzrange=None,mute=False):
        """
        OBSOLETE (use retrieve_scans or sum_scans)
        iterates through a mzML tree and pulls full scan mass spectra

        sr: [start,end] or None
            scan range to sum (default 'all')
        tr: [start,end] or None
            time range to sum within
        mzrange: [start,end] or None
            limit the m/z range to extract (default None)
            None will keep whatever is present in the file
            this is a memory saving utility
        
        output:
            dictionary of timpoints with subkeys for m/z, intensity, and scan #
            also outputs the updated scan range and mz range
        
        exclusively pulls mass spectra with MS level 1 (for MSMS spectra use pullmsmsspectra)
        """
        self.sys.exit('The pullspectra function is obsolete. Identify the appropriate function to sum from mzML.functions and use retrieve_scans instead.')
        
    def pulluvspectra(self,sr=None,tr=None):
        """
        OBSOLETE (use retrieve_scans or sum_scans)
        extracts UV spectra from the mzML file
        
        the scan range or time range can be limited by sr and trange respectively (trange overrides sr)
        sr: [start scan, end scan]
        trange: [start time, end time]
        output:
            dictionary of timpoints with subkeys for intensity (corresponding to the uvlambda list) and scan #
        """
        self.sys.exit('The pulluvspectra function is obsolete. Identify the appropriate function to sum from mzML.functions and use retrieve_scans instead.') 
        
    def pwconvert(self,filename,bit=64,compression=True,gzip=True):
        """
        Runs msconvert.exe from ProteoWizard to convert Waters .RAW format to .mzXML
        which can then be parsed by python.
        
        module requirements: os, subprocess, sys
        
        ProteoWizard must be installed for this script to function.
        go to 
        http://proteowizard.sourceforge.net/downloads.shtml
        to download
        
        This script assumes that the ProteoWizard is installed under either
        c:\program files\proteowizard
        or
        c:\program files (x86)\proteowizard
        
        If you use this python script to convert to mzML, you should cite the paper of the folks who wrote the program
        Chambers, M.C. Nature Biotechnology 2012, 30, 918-920
        doi 10.1038/nbt.2377
        """
        def find_all(fname,path):
            """
            Finds all files of a given name within a specified directory.
            Adapted from http://stackoverflow.com/questions/1724693/find-a-file-in-python
            
            Module dependancies: os
            """
            locations = []
            for root,dirs,files in self.os.walk(path):
                if fname in files:
                    locations.append(self.os.path.join(root,fname))                   
            return locations
        
        if self.sys.platform != 'win32':
            raise OSError('The function that converts to mzML is limited to Windows operating systems.\nYou can manually convert to *.mzML using the proteowizard standalone package and supply that mzML file to this script')
        locs = []
        for val in ['c:\\program files\\proteowizard','c:\\program files (x86)\\proteowizard']: #searches for msconvert.exe in expected folders
            locs.extend(find_all('msconvert.exe',val))
                    
        if len(locs)==0: # if script cannot find msconvert.exe
            raise IOError('The python script could not find msconvert.exe\nPlease ensure that ProteoWizard is installed in either:\nc:\\program files\\proteowizard\nor\nc:\\program files (x86)\\proteowizard')
        
        outname = filename[:-4]+'.mzML'
        callstring = locs[-1]+' "'+filename+'" --mzML'
        if bit == 32 or bit == 64:
            callstring += ' --'+str(bit)
        else:
            raise ValueError('ProteoWizard conversion was called with an invalid floating point precision "%s".' %str(bit))
        
        if compression is True: # call for compression
            callstring += ' --zlib'
        
        exten ='*.mzML'
        if gzip is True: # call to gzip entire mzml
            callstring += ' --gzip'
            outname += '.gz'
            exten += '.gz'
        
        import subprocess
        if self.ks['verbose'] is True:
            callstring += ' --verbose'
            self.sys.stdout.write('Generating %s file from %s' %(exten,filename))
            self.sys.stdout.flush()
            subprocess.call(callstring)
            self.sys.stdout.write(' DONE\n')
            self.sys.stdout.flush()  
        else:
            subprocess.call(callstring)
        return outname
    
    def retrieve_scans(self,start=None,end=None,fn=1,mute=False):
        """
        retrieves the specified scans or time range from the specified function
        
        start: integer or float
            the point to start retrieving scans
            if integer, this will be a start scan number
            if float, this will be the start time
        end: (optional) integer or float
            the end point to stop retrieving scans
            same options as start
        fn: integer
            the function to pull scans from (default 1)
        mute: bool
            overrides the verbose setting of the mzml instance
        
        returns a list with each index corresponding to a scan, with two sublists for x and y data
        """
        # find spectrum indicies to extract between
        if fn not in self.functions:
            raise ValueError('The function "%d" is not in this mzml file.' %fn)
        start = self.scan_index(start,fn,bias='greater')
        end = self.scan_index(end,fn,bias='lesser')
        
        if self.ftt is False: # extract the timepoints and etc from the mzml
            self.function_timetic()
        out = []
        for spectrum in self.tree.getElementsByTagName('spectrum'): # go through each spectrum
            attr = self.attributes(spectrum)
            # func,proc,scan = self.fps(spectrum) # determine function, process, and scan numbers
            # p = self.cvparam(spectrum)
            if self.ks['verbose'] is True and mute is False:
                self.sys.stdout.write('\rExtracting scan data from spectrum #%d/%d  %.1f%%' %(attr['index']+1,self.nscans,float(attr['index']+1)/float(self.nscans)*100.))
            if attr['index'] >= start and attr['index'] <= end: # within the index bounds
                out.append(self.extractspectrum(spectrum))
        if self.ks['verbose'] is True and mute is False:
            self.sys.stdout.write(' DONE\n')
        return out
    
    def scan_index(self,scan=None,fn=1,bias='lesser'):
        """
        determines the index for a scan or timepoint in a given function
        scan: integer or float
            the scan number or time point to find
        fn: integer
            the function to look in
        bias: options dictated by locateinlist()
            bias of index finding
        """
        if fn not in self.functions:
            raise KeyError('The function %d is not in this mzML file.' %fn)
        if scan is None: # if no scan number is specified
            if bias == 'greater': # used for start point
                return self.functions[fn]['sr'][0]
            if bias == 'lesser': # used for end point
                return self.functions[fn]['sr'][1]
        if type(scan) is float: # timepoint
            if self.ftt is False:
                self.function_timetic()
            return self.locateinlist(self.functions[fn]['timepoints'],scan,bias=bias) + self.functions[fn]['sr'][0] # return located index plus start of the scan range
        elif type(scan) is int: # scan number
            if scan < 1:
                raise ValueError('The scan number must be greater or equal to 1 (specified: %d)' %scan)
            if scan > self.functions[fn]['nscans']:
                raise ValueError('The scan number %d exceeds the number of scans in function %d (%d)' %(scan,fn,self.functions[fn]['nscans']))
            return scan - 1 + self.functions[fn]['sr'][0] # return scan minus 1 (to shift into index domain) plus the start location index
        else:
            raise ValueError('An unexpected scan type was handed to the scan_index function ("%s", type: %s)' %(str(scan),type(scan)))
    
    def scanrange_from_timerange(self,start,end=None,fn=1):
        """
        OBSOLETE (use scan_index)
        determines the scan range implied by the bounds of a supplied time range
        start: start time
        end: end time (optional)
        fn: function
        
        returns [start scan, end scan] (if end time was not specified, end scan will be None)
        """
        self.sys.exit('The scanrange_from_timerange function is obsolete, use scan_index')
   
    def scanproperties(self,hand):
        """determines the scan properties of the provided spectrum"""
        mstypes = { # ms accession keys and their respective names (for spectrum identification)
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
        othertypes = { # other accession keys (non-MS)
        'MS:1000620': 'PDA spectrum',
        'MS:1000804': 'electromagnetic radiation spectrum',
        'MS:1000805': 'emission spectrum',
        'MS:1000806': 'absorption spectrum',
        }
        out = {}
        if isinstance(hand,self.cvparam): # handed a cvparam class object (expected)
            p = hand
        else: # handed a tree or branch (generate the cvparam class object)
            p = self.cvparam(hand)
        for acc in p:
            if acc in mstypes: # if scan is a type of mass spectrum
                out['acc'] = str(acc) # accession code
                out['name'] = mstypes[acc] # name of spectrum
                out['type'] = 'MS' # it is a mass spectrum
                out['level'] = p['MS:1000511'] # ms level
                out['window'] = [p['MS:1000501'],p['MS:1000500']] # scan window
                if p.has_key('MS:1000129'): # negative scan
                    out['mode'] = '-'
                elif p.has_key('MS:1000130'): # positive scan
                    out['mode'] = '+'
                if out['level'] == 2: # if msms
                    out['target'] = p['MS:1000827'] # isolation window target m/z
                elif out['level'] > 2: # if MSn > 2, not sure how to handle this (will have to be hard coded later as I have no examples)
                    raise ValueError('This script has not been coded to handle MSn > 2, please contact the author of the class')
                
            if acc in othertypes: # if the scan is something else
                out['acc'] = acc # accession code
                out['name'] = othertypes[acc] # name of spectrum
                if p.has_key('MS:1000804'): # if it is a UV-Vis
                    out['type'] = 'UV'
                else: # other other type (not handled by script)
                    raise KeyError('The script has not been coded to handle spectra types other than MS and UV-Vis. Please contact the authors to get this functionality included.')
        return out     
    
    def sum_scans(self,start=None,end=None,fn=1,dec=3,mute=False):
        """
        sums the specified scans together
        if the scan range moves into another function, an error is raised
        this function has a lower memory overhead than retrieve_scans()
        
        start: integer or float
            start point to begin summing
            integer is treated as scan number
            float is treated as a time point
        end: integer or float
            end point to finish summing
            as with start
        fn: function to look at
            default 1
        dec: int
            number of decimal places to track in the spectrum (lower values lower memory overhead)
            this is only relevant when summing spectra together
        mute: bool
            override for verbose toggle of mzml instance
        
        output: [xlist,ylist]
        """
        if self.functions[fn]['type'] != 'MS':
            raise ValueError('The sum_scans function does not have the functionality to sum non-mass spec scans.')
        start = self.scan_index(start,fn,'greater')
        end = self.scan_index(end,fn,'lesser')
        
        from _Spectrum import Spectrum
        spec = Spectrum(dec,self.functions[fn]['window'][0],self.functions[fn]['window'][1]) # create Spectrum object
        
        for spectrum in self.tree.getElementsByTagName('spectrum'): # go through each spectrum
            attr = self.attributes(spectrum) # get attributes
            if attr['index'] > end:
                break
            if self.ks['verbose'] is True and mute is False:
                self.sys.stdout.write('\rCombining spectrum #%d (scan range: %d-%d)  %.1f%%' %(attr['index']+1,start,end,(float(attr['index']-start))/(float(end-start))*100.))
            if attr['index'] >= start and attr['index'] <= end: # if within the specified bounds
                x,y = self.extractspectrum(spectrum) # pull spectrum
                spec.addspectrum(x,y) # add spectrum to Spectrum object
        out = spec.trim()
        if self.ks['verbose'] is True and mute is False:
            self.sys.stdout.write(' DONE\n')
        return out

    def trimspectrum(self,x,y,left,right):
        """trims a spectrum to the left and right bounds"""
        l,r = self.locateinlist(x,left,'greater'),self.locateinlist(x,right,'lesser') # find indicies
        return x[l:r],y[l:r] # trim spectrum
                    
if __name__ == '__main__':
    filename = 'BTM-42 UVVis'
    mzml = mzML(filename,verbose=True)
    #from _Spectrum import Spectrum
    #sp = {
    #'pos':{'bounds':[325,327],'affin':'+','spectrum':Spectrum(3),'raw':[]},
    #'neg':{'bounds':[348,350],'affin':'-','spectrum':Spectrum(3),'raw':[]},
    #'uv':{'bounds':[378,None],'affin':'UV','raw':[]}
    #}