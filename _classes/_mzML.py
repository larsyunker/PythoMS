"""
Class for interpreting and processing mzML files

CHANGELOG:
    renamed binarytolist to extractspectrum
    incorporated a kwarg into extractspectrum to also return the x and y unit on request
    fixed checkforfile (several errors and missed catches in the previous version)
    added __getitem__ to sum the specified scans (only sums MS1 spectra)
    class now accepts keyword arguments (only used for calling of pwconvert and verbose currently)
    reordered keywords to be assigned before referencing
    added obo class for loading and parsing an obo file (for accession codes that are not hard coded into the mzml file)
    cvparam now returns a dictionary of accession keys (not names) to make end user coding easier to specify exactly what they want to retrieve
    added accession key from name finder to obo (for ease of finding accession key if you know the name)
    removed loops from numberofthings (loops were unnecessary)
    updated str and repr outputs
    removed spectrumList and chromatogramList loops (these were redundant, as looping through the spectrum or chromatogram tags accomplishes the same thing)
    renamed scantype to spectrumtype (scan implies ms, where it is not necessarily)
    added _foreachscan and _foreachchrom decorator functions (which will apply the provided function to every scan or chromatogram in the mzml)
    moved cvparam method and units method into cvparam class
    ---2.3

to add:
    fix getitem to actually get the spectrum specified regardless of its ms level or whatever
    add plotspectrum call (which would use the tome_v02 plotter)
    identify spectrometer in softwareList
    obtain example files from different manufacturers and validate the script with those
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
        self.bisect = __import__('bisect')
        self.zlib = __import__('zlib')
        self.bl = self.bisect.bisect_left # for convenience of calls
        self.br = self.bisect.bisect_right
        self.sys.path.append(self.os.path.dirname(self.os.path.realpath(__file__))) # yes this is needed. don't be stupid.
        
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
        self.nscans,self.nchroms = self.numberofthings() # find number of scans and number of chromatograms
        if self.ks['verbose'] is True:
            self.sys.stdout.write(' DONE\n')
        
    def __str__(self):
        """The string that is returned when printed"""
        return 'loaded mzML file: "%s" (#spectra: %d, #chromatograms: %d)' %(self.filename, self.nscans, self.nchroms)
    
    def __repr__(self):
        """The representation that is returned"""
        return "%s('%s')" %(self.__class__.__name__,self.filename)
    
    def __getitem__(self,ind):
        """returns the summed scans with the supplied index/indicies"""
        if isinstance(ind,slice):
            if ind.start is None: # no start
                start = 0
            else:
                start = ind.start
            if ind.stop is None: # no stop
                stop = self.nscans
            else:
                stop = ind.stop
        else:
            start = ind
            stop = ind
        if start < 0:
            start = self.nscans+start+1
        if stop < 0:
            stop = self.nscans+stop+1
        return self.sumscans(sr=[start,stop])[0]
    
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
        
        def __getitem__(self,key):
            """
            returns the value of the accesssion key supplied
            (can also be handed an item name, provided that it is defined in this cvParameter set)
            """
            if isinstance(key,slice):
                raise ValueError('Slicing of a cvparam is not supported')
            if type(key) != str:
                print key, type(key)
                raise ValueError('Only accession keys or names may be retrieved from the cvparam object')
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
            out = {}
            for cvParam in branch.getElementsByTagName('cvParam'):
                acc = cvParam.getAttribute('accession') # accession key
                out[acc] = {}
                for attribute,value in cvParam.attributes.items(): # pull all the attributes
                    if attribute != 'accession':
                        out[acc][attribute] = value
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
                            
    def attributes(self,branch):
        """pulls all attributes of a supplied branch and creates a dictionary of them"""
        out = {}
        for pair in branch.attributes.items():
            out[pair[0]] = self.stringtodigit(pair[1])
        return out
    
    def autoresolution(self,n=10):
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
        from _Spectrum import Spectrum
        import scipy as sci
        ranges = [] # list of scan intervals
        for i in range(n): # generate 10 pseudo-random intervals to sample
            ran = random()
            if (ran*self.nscans)-10 >= 0 and (ran*self.nscans)+10 <= self.nscans:
                ranges.append([int((ran*self.nscans)-10),int((ran*self.nscans)+10)])
        summed = []
        for ind,rng in enumerate(ranges):
            if self.ks['verbose'] is True:
                self.sys.stdout.write('\rEstimating resolution of the instrument %.0f%%' %(float(ind+1)/float(n)*100.))
            spectra,sr,mzrange = self.pullspectra(rng,mute=True) # pull the spectra in the scan range
            spectrum = Spectrum(2,mzrange[0],mzrange[1]) # generate a Spectrum object
            for time in spectra: # add each spectrum to the Spectrum object
                spectrum.addspectrum(spectra[time]['x'],spectra[time]['y'])
            summed.append(spectrum.trim()) # append the summed spectra
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
        if units is not False: # appends the units onto out
            for unit in units:
                out.append(unit)
        return out
    
    def filepresent(self,fn,ty='file'):
        """
        checks for the presence of the specified file or directory in the current working directory
        ty specifies the type of thing to look for "file" for file or "dir" for directory
        """
        if ty == 'file':
            if self.os.path.isfile(fn) == False:
                return False
            else:
                return True
        if ty == 'dir':
            if self.os.path.isdir(fn) == False:
                return False
            else:
                return True
        
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
                if self.filepresent(fn+oopsg[key],'file') is True:
                    return fn+oopsg[key]
        for key in oopsx: # tries to complete mzml shortenings
            if fn.lower().endswith(key) is True:
                if self.filepresent(fn+oopsx[key],'file') is True:
                    return fn+oopsx[key]
        for key in oopsr: # tries to complete raw shortenings
            if fn.lower().endswith(key) is True:
                if self.filepresent(fn+oopsr[key],'dir') is True:
                    return fn+oopsr[key]
        if self.filepresent(fn+'.raw','dir') is True: # finally looks for raw file
            return fn+'.raw'
        raise IOError('The file "%s" could not be located in the current working directory'%(fn)) # if it can't be found, raise IOError
    
    def integ(self,name,start,end,x,y):
        """
        Function to integrate values between paired list indicies (e.g. a m/z list and an intensity list)
        
        name: name of the peak being integrated (only used for warning purposes)
        start: start x value
        end: end x value
        x: list of x values
        y: list of y values (paired with x)
        """
        if start > max(x) or start < min(x): # check that start is within the m/z bounds
            self.BE.warn(name,start,end,min(x),max(x))
        if end is None: # if only a start value is supplied, return closest to that value
            return y[self.takeclosest(x,start)]
        if end > max(x): # check that end is within the m/z bounds
            self.BE.warn(name,start,end,min(x),max(x))
        return sum(y[self.br(x,start):self.bl(x,end)]) # integrate using the nearest values inside the bounds        
    
    def numberofthings(self):
        """retrieves the number of scans and chromatograms in the file"""
        # I'm not sure if the following if statements will ever be true, but they'll break the class if they are
        if len(self.tree.getElementsByTagName('spectrumList')) > 1:
            raise ValueError("There are more than one set of scans, and this script can't handle it.\nGive this file to Lars so he can figure out how to fix it.")
        if len(self.tree.getElementsByTagName('chromatogramList')) > 1:
            raise ValueError("There are more than one set of chromatograms, and this script can't handle it.\nGive this file to Lars so he can figure out how to fix it.")
        nscans = int(self.tree.getElementsByTagName('spectrumList')[0].getAttribute('count')) # number of spectra
        nchroms = int(self.tree.getElementsByTagName('chromatogramList')[0].getAttribute('count')) # number of chromatograms
        return nscans,nchroms
    
    def pullchromdata(self):
        """Pulls mzML chromatograms"""
        chroms = {} #dictionary of chromatograms
        for chromatogram in self.tree.getElementsByTagName('chromatogram'):
            attr = self.attributes(chromatogram) # pull attributes
            #p = self.cvparam(chromatogram) # pull parameters
            if self.ks['verbose'] is True:
                self.sys.stdout.write('\rExtracting chromatogram #%s/%i  %.1f%%' %(attr['index']+1,self.nchroms,float(attr['index']+1)/float(self.nchroms)*100.))
                self.sys.stdout.flush()
            x,y,xunit,yunit = self.extractspectrum(chromatogram,units=True) # extract x list, y list, and units
            chroms[attr['id']] = {'x':x, 'y':y, 'xunit':xunit, 'yunit':yunit}
        
        if self.ks['verbose'] is True:
            self.sys.stdout.write(' DONE\n')
        return chroms

    def pullspeciesdata(self,sp,sumspec=False):
        """
        Extracts integrated data at every timepoint for all species specified in the sp dictionary
        
        sp: dictionary
            one key for every species to track
            format for each key: 'peak name':{'bounds':[peak start mz,peak end mz],'affin':['+' or '-' or 'UV'},'raw':[],'spectrum':[]
        sumspec: bool
            toggles summing of all spectra together (creates an additional output item)
        
        output:
            filled dictionary, dictionary of total ion currents, dictionary of retention times
            if sumspec is true, will also output an [x,y] list
        
        explicitly interprets full scan mass spectra and UV species
        """
        rtime = {} # generate empty lists required for data processing
        TIC = {}
        if sumspec is True:
            from _Spectrum import Spectrum
            spec = Spectrum(3)
        self.BE = self.BoundsError() # load warning instance for integration
        for spectrum in self.tree.getElementsByTagName('spectrum'):
            attr = self.attributes(spectrum) # get attributes
            p = self.cvparam(spectrum) # pull parameters of the scan
            mode,level = self.spectrumtype(p) # determine the scan type
            if self.ks['verbose'] is True:
                self.sys.stdout.write('\rExtracting species data from spectrum #%d/%d  %.1f%%' %(attr['index']+1,self.nscans,float(attr['index']+1)/float(self.nscans)*100.))
            if mode is not None and level < 2:
                modekey = 'raw'+mode # define dictionary key for current scan
                if rtime.has_key(modekey) is False: # create dictionary entry if not present
                    rtime[modekey] = []
                if TIC.has_key(modekey) is False:
                    TIC[modekey] = []
                TIC[modekey].append(int(p['MS:1000285'])) # append total ion current value
                rtime[modekey].append(float(p['MS:1000016'])) # append scan start time
                x,y = self.extractspectrum(spectrum) # generate spectrum
                if sumspec is True:
                    spec.addspectrum(x,y)
                for key in sp: # integrate each peak
                    if sp[key]['affin'] == mode: # if species has affinity to this spectrum type
                        if mode in ['+','-']: # if mass spectrum
                            sp[key]['raw'].append(self.integ(key,sp[key]['bounds'][0],sp[key]['bounds'][1],x,y)) # integrate
                            xt,yt = self.trimspectrum(x,y,sp[key]['bounds'][0],sp[key]['bounds'][1]) # trim spectrum
                            sp[key]['spectrum'].addspectrum(xt,yt) # add spectrum
                        if mode in ['UV']: # if UV spectrum
                            sp[key]['raw'].append(self.integ(key,sp[key]['bounds'][0],sp[key]['bounds'][1],x,y)/1000000.) # integrates and divides by 1 million bring it into au
        
        for key in sp: #remove mz/int values in spectrum that are None
            if sp[key]['affin'] in ['+','-']:
                sp[key]['spectrum'] = sp[key]['spectrum'].trim() # trim spectrum to remove Nonetypes        
        self.TIC = TIC
        self.rtime = rtime
        if self.ks['verbose'] is True:
            self.sys.stdout.write(' DONE\n')
        self.BE.printwarns() # print bounds warnings (if any)
        if sumspec is True:
            return sp,TIC,rtime,spec.trim()  
        else:
            return sp,TIC,rtime
    
    def pullspectra(self,sr=None,mzrange=None,mute=False):
        """
        iterates through a mzML tree and pulls full scan mass spectra

        sr: [start,end] or None
            scan range to sum (default 'all')
        mzrange: [start,end] or None
            limit the m/z range to extract (default None)
            None will keep whatever is present in the file
            this is a memory saving utility
        
        output:
            dictionary of timpoints with subkeys for m/z, intensity, and scan #
            also outputs the updated scan range and mz range
        
        exclusively pulls mass spectra with MS level 1 (for MSMS spectra use pullmsmsspectra)
        """
        if sr is None:
            sr = [1,self.nscans]
        if sr[1] < sr[0]:
            raise ValueError('The scan range is invalid: %d-%d' %(sr[0],sr[1]))
        speclist = {} # dictionary for spectra
        
        for spectrum in self.tree.getElementsByTagName('spectrum'): # go through each spectrum
            attr = self.attributes(spectrum) # get attributes
            p = self.cvparam(spectrum) # pull parameters of the scan
            if self.ks['verbose'] is True and mute is False:
                self.sys.stdout.write('\rExtracting mass spectrum #%s (scan range: %d-%d)  %.1f%%' %(attr['index']+1,sr[0],sr[1],(float(attr['index']-sr[0]+1))/(float(sr[1]-sr[0]))*100.))
            mode,level = self.spectrumtype(p) # determine the scan type
            if mode in ['+','-'] and level <2: # if type is full scan mass spectrum
                if attr['index']+1 >= sr[0] and attr['index']+1 <= sr[1]:
                    x,y = self.extractspectrum(spectrum)
                    if mzrange is not None: # if m/z range is specified, trim spectrum to specified range
                        x,y = self.trimspectrum(x,y,mzrange[0],mzrange[1])
                    speclist[float(p['MS:1000016']['value'])] = {'x':x,'y':y,'scan':attr['index']} # key by scan start time
        if mzrange is None: # determine m/z range if not specified
            minx = 1000000.
            maxx = 0.
            for scan in speclist:
                if min(speclist[scan]['x']) < minx:
                    minx = min(speclist[scan]['x'])
                if max(speclist[scan]['x']) > minx:
                    maxx = max(speclist[scan]['x'])
            mzrange = [minx,maxx]
        if self.ks['verbose'] is True and mute is False:
            self.sys.stdout.write(' DONE\n')
        return speclist,sr,mzrange
    
    def pullmsmsspectra(self):
        """
        Extracts MSMS spectra from the mzML file
        
        Exclusively finds mass spectra with level greater or equal to 2
        groups spectra by key defined by thge isolation window target
        
        each dictionary has subkeys for each time point which each has subkeys containing the information of that scan
        
        also returns a dictionary of m/z limits using the same keys as msms
        """
        msms = {}
        limits = {}
        for spectrum in self.tree.getElementsByTagName('spectrum'): # go through each spectrum
            attr = self.attributes(spectrum)
            p = self.cvparam(spectrum)
            if self.ks['verbose'] is True:
                self.sys.stdout.write('\rExtracting species data from spectrum #%d/%d  %.1f%%' %(attr['index']+1,self.nscans,float(attr['index']+1)/float(self.nscans)*100.))
            mode,level = self.spectrumtype(p)
            if level >= 2: # if it is a msms spectrum
                tic = float(p['MS:1000285']) # total ion current
                t = float(p['MS:1000016']) # scan start time
                ce = float(p['MS:1000045']) # collision energy
                target = float(p['MS:1000827']) # isolation window target m/z
                lowmz = float(p['MS:1000501']) # scan window lower limit
                highmz = float(p['MS:1000500']) # scan window upper limit
                
                if limits.has_key(target) is False:
                    limits[target] = [lowmz,highmz]
                x,y = self.extractspectrum(spectrum)
                if msms.has_key(target) is False:
                    msms[target] = {}
                msms[target][t] = {'CE':ce,'TIC':tic,'x':list(x),'y':list(y)}
        if self.ks['verbose'] is True:
            self.sys.stdout.write(' DONE\n')
        return msms,limits
    
    def pulluvspectra(self):
        """
        extracts UV spectra from the mzML file
        
        returns a list of timepoints, a list of wavelengths, and a list of lists of intensity corresponding to timepoints
        """
        rtime = [] # generate empty lists required for data processing
        uvlambda = None
        uvint = []
        for spectrum in self.tree.getElementsByTagName('spectrum'): # go through each spectrum
            attr = self.attributes(spectrum)
            p = self.cvparam(spectrum)
            if self.ks['verbose'] is True:
                self.sys.stdout.write('\rExctracting UV spectrum #%i/%i  %.1f%%' %(attr['index']+1,self.nscans,float(attr['index']+1)/float(self.nscans)*100.))
            mode,level = self.spectrumtype(p)
            if mode == 'UV': # if type is UV-Vis
                rtime.append(p['MS:1000016']) # scan start time
                x,y = self.extractspectrum(spectrum)
                if uvlambda is None: # if wavelength region has not yet been defined
                    uvlambda = list(x)
                for ind,val in enumerate(y): # normalize y value by 1 million to bring value into a.u.
                    y[ind] = val/1000000.
                uvint.append(y) # append intensity values
        if self.ks['verbose'] is True:
            self.sys.stdout.write(' DONE\n')
        return rtime,uvlambda,uvint    

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
        
        If you use this python script, you should cite the paper of the folks who wrote the program
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
            raise ValueError('ProteoWizard conversion was called with an invalid bit precision "%s".' %str(bit))
        
        if compression is True: # call for compression
            callstring += ' --zlib'
        
        if gzip is True: # call to compress entire mzml
            callstring += ' --gzip'
            outname += '.gz'
        
        import subprocess
        if self.ks['verbose'] is True:
            callstring += ' --verbose'
            self.sys.stdout.write('Generating *.mzML file from *.raw')
            self.sys.stdout.flush()
            subprocess.call(callstring)
            self.sys.stdout.write(' DONE\n')
            self.sys.stdout.flush()  
        else:
            subprocess.call(callstring)
        return outname
    
    def spectrumtype(self,hand):
        """determines the scan type of the provided spectrum"""
        if isinstance(hand,self.cvparam): # handed a cvparam class object (expected)
            p = hand
        else: # handed a tree or branch (generate the cvparam class object)
            p = self.cvparam(hand)
        if p.has_key('MS:1000579') or p.has_key('MS:1000580'): # MS1 spectrum (full scan MS spectrum) or MSn spectrum (MS/MS)
            level = int(p['MS:1000511']) # ms level
            if p.has_key('MS:1000129'): # negative scan
                mode = '-'
            elif p.has_key('MS:1000130'): # positive scan
                mode = '+'
            return mode,level
        elif p.has_key('MS:1000804'): # otherwise is it a UV spectrum
            return 'UV',0
        else: # type that is not currently handled by the script
            return None,None       
    
    def stringtodigit(self,string):
        """attempts to convert a unicode string to float or integer"""
        try:
            value = int(string)
        except ValueError:
            try:
                value = float(string)
            except ValueError:
                value = string
        return value    
    
    def sumscans(self,sr=None,mzrange=[50.,2000.],dec=3,mute=False):
        """
        sums all ms1 full spectrum scans together
        this function has a lower memory overhead than pullscans()
        
        sr: [start,end] or None
            scan range to sum (default None)
            None will sum all scans
        mzrange: [start,end]
            mz range to sum between
        dec: int
            number of decimal places to track in the spectrum (lower values lower memory overhead)
        """
        from _Spectrum import Spectrum
        spec = Spectrum(dec, startmz=mzrange[0], endmz=mzrange[1])
        if sr == None:
            sr = [1,self.nscans]
        if sr[1] < sr[0]:
            raise ValueError('The supplied scan range is invalid: %d-%d' %(sr[0],sr[1]))
        for spectrum in self.tree.getElementsByTagName('spectrum'): # go through each spectrum
            attr = self.attributes(spectrum) # get attributes
            if attr['index']+1 > sr[1]:
                break
            p = self.cvparam(spectrum) # pull parameters of the scan
            mode,level = self.spectrumtype(p)
            if self.ks['verbose'] is True and mute is False:
                if sr[1]-sr[0] != 0:
                    self.sys.stdout.write('\rCombining mass spectrum #%d (scan range: %d-%d)  %.1f%%' %(attr['index']+1,sr[0],sr[1],(float(attr['index']-sr[0]+1))/(float(sr[1]-sr[0]))*100.))
                else:
                    self.sys.stdout.write('\rExtracting scan #%d' %sr[0])
            if mode in ['+','-'] and level <2: # if type is full scan mass spectrum
                if attr['index']+1 >= sr[0] and attr['index']+1 <= sr[1]:
                    x,y = self.extractspectrum(spectrum)
                    if mzrange is not None: # if m/z range is specified, trim spectrum to specified range
                        x,y = self.trimspectrum(x,y,mzrange[0],mzrange[1])
                    spec.addspectrum(x,y)
        out = spec.trim()
        if self.ks['verbose'] is True and mute is False:
            self.sys.stdout.write(' DONE\n')
        return out,sr

    def takeclosest(self,lst, value):
        """
        Finds index in a sorted list of the value closest to a given value
        
        If two numbers are equally close, return the smallest number.
        based on http://stackoverflow.com/questions/12141150/from-list-of-integers-get-number-closest-to-a-given-value
        
        input:
            lst list of values
            value number to find
        """
        pos = self.bl(lst, value)
        
        if pos == 0: # if at start of list
            return pos
        elif pos == len(lst):
            return pos # if at end of list
        else:
            prev = abs(lst[pos-1] - value) # difference between value and previous index
            cur = abs(lst[pos] - value) # difference between value and current index
            if prev < cur: # if it is closer to previous value, return that index
                return pos-1
            else:
                return pos
    
    def trimspectrum(self,x,y,left,right):
        """trims a spectrum to the left and right bounds"""
        l,r = self.bl(x,left),self.br(x,right) # find indicies
        return x[l:r],y[l:r] # trim spectrum
    
    #def units(self,p):
    #    """
    #    takes a cvparam dictionary and returns the appropriate unit
    #    (expects to be handed the cvparam dictionary of a binaryDataArray)
    #    all units as of 2016-07-12 are defined here, but all possible units can be found in
    #    https://raw.githubusercontent.com/HUPO-PSI/psi-ms-CV/master/psi-ms.obo
    #    """
    #    unitkeys = {
    #    # x units
    #    'MS:1000040':'m/z',
    #    'UO:0000010':'second',
    #    'UO:0000017':'micrometer',
    #    'UO:0000018':'nanometer',
    #    'UO:0000028':'millisecond',
    #    'UO:0000031':'minute',
    #    'UO:0000221':'dalton',
    #    'UO:0000222':'kilodalton',
    #    # y units
    #    'MS:1000131':'number of detector counts',
    #    'MS:1000132':'percent of base peak',
    #    'MS:1000814':'counts per second',
    #    'MS:1000905':'percent of base peak times 100',
    #    'UO:0000187':'percent',
    #    'UO:0000269':'absorbance unit',
    #    # other
    #    'MS:1000807':'Th/s',
    #    'UO:0000002':'mass unit',
    #    'UO:0000008':'meter',
    #    'UO:0000012':'kelvin',
    #    'UO:0000021':'gram',
    #    'UO:0000027':'degree Celsius',
    #    'UO:0000098':'milliliter',
    #    'UO:0000106':'hertz',
    #    'UO:0000110':'pascal',
    #    'UO:0000112':'joule',
    #    'UO:0000166':'parts per notation unit',
    #    'UO:0000169':'parts per million',
    #    'UO:0000175':'gram per liter',
    #    'UO:0000185':'degree',
    #    'UO:0000218':'volt',
    #    'UO:0000228':'tesla',
    #    'UO:0000266':'electronvolt',
    #    'UO:0000268':'volt per meter',
    #    }
    #    for key in p:
    #        if p[key].has_key('unitAccession'): # find unitAccession (the accession code defines the unit)
    #            if p[key]['unitAccession'] in unitkeys: # if it is defined, return unit code
    #                return unitkeys[p[key]['unitAccession']]
    #            else:
    #                if self.__dict__.has_key('loadedobo') is False:
    #                    self.loadedobo = self.obo(self.ks['obo']) # load obo file
    #                return self.loadedobo[p[key]['unitAccession']] # return accession key name as defined in obo file
                    
if __name__ == '__main__':
    filename = 'HZ-140516_HOTKEYMSMS 1376 II'
    mzml = mzML(filename,verbose=True)
    #from _Spectrum import Spectrum
    #sp = {
    #'pos':{'bounds':[325,327],'affin':'+','spectrum':Spectrum(3),'raw':[]},
    #'neg':{'bounds':[348,350],'affin':'-','spectrum':Spectrum(3),'raw':[]},
    #'uv':{'bounds':[378,None],'affin':'UV','raw':[]}
    #}