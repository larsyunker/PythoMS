"""
Class for interpreting and processing mzML files

CHANGELOG:

new:
    incorporated conversion to mzml into the class
    moved base64 and struct imports into init
    added verbose call to class (and moved sys import to init)
    moved find_all into pwconvert
    ---1.1---
    moved bisect import into init (avoided use of arrays because searchsorted is an order of magnitude slower than bisect when floats are involved)
    removed .mzXML functionality (very outdated)
    commented, cleaned up, and lower-cased takeClosest
    added catches for incomplete file extensions
    changed sys.exit to exceptions in pwconvert
    created pullmsmsspectra which searches for and pulls spectra with MS level >=2
    ---1.2---
    cleaned up and fixed checkforfile (I think most exceptions and eventualities have been accounted for now)
    added trimspectrum() method
    added sumscans() method (sums all scans together)
    ---1.3---
    changed integrate to warn instead of exception raise in the event that the bounds exceed the m/z range
    created a BoundsError subclass to handling bounds warnings
    added number of scans retrieval (used for auto resolution)
    added automatic resolution calculator
    ---1.4---
    created cvparams and attributes methods to pull all parameters and attributes easily
    consolidated gettext and decode into a single function that pulls binary strings and converts them
    rewrote pull functions to use the consolidated methods
    elaborated output of auto resolution to give progress
    ---2.0---

to add:
    update pullspeciesdata docstring to better represent the newest version of the supplied dictionary
"""

class mzML(object):
    def __init__(self,filename,verbose=True):
        """
        Class for interpreting an mzML (mass spectrum) file
        """
        self.v = verbose
        self.sys = __import__('sys')
        import os
        self.sys.path.append(os.path.dirname(os.path.realpath(__file__)))
        self.filename = self.checkforfile(filename)
        self.b64 = __import__('base64')
        self.st = __import__('struct')
        self.bisect = __import__('bisect')
        self.bl = self.bisect.bisect_left # for convenience of calls
        self.br = self.bisect.bisect_right
        if self.v is True:
            self.sys.stdout.write('Loading %s into memory' %self.filename)
            self.sys.stdout.flush()
        import xml.dom.minidom
        try:
            self.tree = xml.dom.minidom.parse(self.filename) # full mzML file
        except:
            raise IOError('The mzML file "%s" could not be loaded. The file is either corrupt or incomplete.' %self.filename)
        self.BE = self.BoundsError() # load warning instance for integration
        self.nscans,self.nchroms = self.numberofthings() # find number of scans and number of chromatograms
        if self.v is True:
            self.sys.stdout.write(' DONE\n')
        
    def __str__(self):
        """The string that is returned when printed"""
        return "The loaded mzML file is '{}'".format(self.filename)
    
    def __repr__(self):
        """The representation that is returned"""
        return "{}('{}')".format(self.__class__.__name__,self.filename)
    
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
            if self.v is True:
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
        if self.v is True:
            self.sys.stdout.write(' DONE\n')
        res = [y for y in res if y is not None] # removes None values (below S/N)
        return sum(res)/len(res) # return average
        
    def binarytolist(self,spectrum):
        """pulls and converts binary data to list"""
        def decode(line,extension='.mzML'):
            """
            decodes base32 strings contained in mzML files
            Currently only functions with 32-bit precision
            """
            decoded = self.b64.decodestring(line) #decodes base64 string 
            unpack_format = "<%dL" % speclen # little-endian, number of values, unsigned long
            values = [] 
            for tmp in self.st.unpack(unpack_format,decoded):
                tmp_i = self.st.pack("I",tmp) #pack as unsigned integer I
                tmp_f = self.st.unpack("f",tmp_i)[0] #unpack as float f
                values.append(float(tmp_f))
            return values
        
        def gettext(nodelist):
            """gets text from a simple XML object"""
            rc = []
            for node in nodelist:
                if node.nodeType == node.TEXT_NODE:
                    rc.append(node.data)
            return ''.join(rc)
        
        speclen = int(spectrum.getAttribute('defaultArrayLength')) # spectrum length (defined in the spectrum attricubes)
        B32 = []
        for binary in spectrum.getElementsByTagName('binaryDataArray'): # pull both binary strings
            B32.append(gettext(binary.getElementsByTagName('binary')[0].childNodes))
        return [decode(B32[0]),decode(B32[1])] # return decoded binary data
    
    def checkforfile(self,fn):
        """checks for file and converts if necessary"""
        if fn.endswith('.raw') is True:
            if self.filepresent(fn[:-4]+'.mzML') is False:
                if self.filepresent(fn,'dir') is True:
                    self.pwconvert(fn)
                    return fn[:-4]+'.mzML'
                else:
                    raise IOError('The raw file "%s" is not in the current working directory. Please check your input' %fn)
            else:
                return fn[:-4]+'.mzML'
        elif fn.lower().endswith('.mzml') is True:
            if self.filepresent(fn) is True:
                return fn
            else:
                raise IOError('The mzML file "%s" is not in the current working directory. Please check your input' %fn)
        else:
            fn = self.fixextension(fn)
            if fn.endswith('.raw') is True:
                self.pwconvert(fn)
                return fn[:-4]+'.mzML'
            return fn
    
    def cvparam(self,branch):
        """
        retrieves the values of each cvParam in the branch
        
        # there is currently no need to retrieve all attributes of the cvParams,
        # but they can be extracted in dictionary format
        for cvParam in branch.getElementsByTagName('cvParam'):
            out[cvParam.getAttribute('name')] = {}
            for attribute,value in cvParam.attributes.items():
                out[cvParam.getAttribute('name')][attribute] = value
        """
        out = {}
        for cvParam in branch.getElementsByTagName('cvParam'):
            out[cvParam.getAttribute('name')] = self.stringtodigit(cvParam.getAttribute('value'))
        return out
    
    def filepresent(self,fn,ty='file'):
        """
        checks for the presence of the specified file or directory in the current working directory
        ty specifies the type of thing to look for "file" for file or "dir" for directory
        """
        import os
        if ty == 'file':
            if os.path.isfile(fn) == False:
                return False
            else:
                return True
        if ty == 'dir':
            if os.path.isdir(fn) == False:
                return False
            else:
                return True
        
    def fixextension(self,fn):
        """tries to fix invalid file extensions"""
        oopsx = {'.mzm':'l','.mz':'ml','.m':'zml','.':'mzml'} # incomplete mzml extensions
        oopsr = {'.ra':'w','.r':'aw','.':'raw'} # incomplete raw extionsions
        for key in oopsx: # tries to complete mzml shortenings
            if fn.lower().endswith(key) is True:
                if self.filepresent(fn+oopsx[key],'file') is True:
                    return fn+oopsx[key]
        for key in oopsr: # tries to complete raw shortenings
            if fn.lower().endswith(key) is True:
                if self.filepresent(fn+oopsr[key],'dir') is True:
                    return fn+oopsr[key]
        if self.filepresent(fn+'.mzml') is True: # tries to add the extension
            return fn+'.mzml'
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
        nscans = 0
        nchroms = 0
        if len(self.tree.getElementsByTagName('spectrumList')) > 1:
            raise ValueError("There are more than one set of scans, and this script can't handle it.\nGive this file to Lars so he can figure out how to fix it.")
        for spectrumList in self.tree.getElementsByTagName('spectrumList'): #for each spectrum list
            if nscans < int(spectrumList.getAttribute('count')):
                nscans = int(spectrumList.getAttribute('count'))
        for chromatogramList in self.tree.getElementsByTagName('chromatogramList'):
            if nchroms < int(chromatogramList.getAttribute('count')):
                nchroms = int(chromatogramList.getAttribute('count'))
        return nscans,nchroms
    
    def pullchromdata(self):
        """Pulls mzML chromatograms"""
        chroms = {} #dictionary of chromatograms
        for chromatogramList in self.tree.getElementsByTagName('chromatogramList'):
            for chromatogram in chromatogramList.getElementsByTagName('chromatogram'):
                attr = self.attributes(chromatogram)
                if self.v is True:
                    self.sys.stdout.write('\rExtracting chromatogram #%s/%i  %.1f%%' %(attr['index']+1,self.nchroms,float(attr['index']+1)/float(self.nchroms)*100.))
                    self.sys.stdout.flush()
                x,y = self.binarytolist(chromatogram)
                """
                currently the x and y units are hard-coded (there is no easy way to tell which cvParam corresponds to which spectrum)
                since no chromatograms seem to have units other than these, it should provide no problems, but it could be changed in the following line
                """
                chroms[attr['id']] = {'x':x, 'y':y, 'xunit':'minute', 'yunit':'number of counts'}
        
        if self.v is True:
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
        
        for spectrumList in self.tree.getElementsByTagName('spectrumList'):
            for spectrum in spectrumList.getElementsByTagName('spectrum'):
                attr = self.attributes(spectrum) # get attributes
                p = self.cvparam(spectrum) # pull parameters of the scan
                mode,level = self.scantype(p) # determine the scan type
                if self.v is True:
                    self.sys.stdout.write('\rExtracting species data from spectrum #%d/%d  %.1f%%' %(attr['index']+1,self.nscans,float(attr['index']+1)/float(self.nscans)*100.))
                if mode is not None and level < 2:
                    modekey = 'raw'+mode # define dictionary key for current scan
                    if rtime.has_key(modekey) is False: # create dictionary entry if not present
                        rtime[modekey] = []
                    if TIC.has_key(modekey) is False:
                        TIC[modekey] = []
                    TIC[modekey].append(int(p['total ion current'])) # append TIC
                    rtime[modekey].append(float(p['scan start time'])) # append scan time
                    x,y = self.binarytolist(spectrum) # generate spectrum
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
        if self.v is True:
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
        
        for spectrumList in self.tree.getElementsByTagName('spectrumList'): #for each spectrum list
            for spectrum in spectrumList.getElementsByTagName('spectrum'): # go through each spectrum
                attr = self.attributes(spectrum) # get attributes
                p = self.cvparam(spectrum) # pull parameters of the scan
                if self.v is True and mute is False:
                    self.sys.stdout.write('\rExtracting mass spectrum #%s (scan range: %d-%d)  %.1f%%' %(attr['index']+1,sr[0],sr[1],(float(attr['index']-sr[0]+1))/(float(sr[1]-sr[0]))*100.))
                mode,level = self.scantype(p) # determine the scan type
                if mode in ['+','-'] and level <2: # if type is full scan mass spectrum
                    if attr['index']+1 >= sr[0] and attr['index']+1 <= sr[1]:
                        x,y = self.binarytolist(spectrum)
                        if mzrange is not None: # if m/z range is specified, trim spectrum to specified range
                            x,y = self.trimspectrum(x,y,mzrange[0],mzrange[1])
                        speclist[float(p['scan start time'])] = {'x':x,'y':y,'scan':attr['index']}
        if mzrange is None: # determine m/z range if not specified
            minx = 1000000.
            maxx = 0.
            for scan in speclist:
                if min(speclist[scan]['x']) < minx:
                    minx = min(speclist[scan]['x'])
                if max(speclist[scan]['x']) > minx:
                    maxx = max(speclist[scan]['x'])
            mzrange = [minx,maxx]
        if self.v is True and mute is False:
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
        for spectrumList in self.tree.getElementsByTagName('spectrumList'): #for each spectrum list
            for spectrum in spectrumList.getElementsByTagName('spectrum'): # go through each spectrum
                attr = self.attributes(spectrum)
                p = self.cvparam(spectrum)
                if self.v is True:
                    self.sys.stdout.write('\rExtracting species data from spectrum #%d/%d  %.1f%%' %(attr['index']+1,self.nscans,float(attr['index']+1)/float(self.nscans)*100.))
                mode,level = self.scantype(p)
                if level >= 2: # if it is a msms spectrum
                    tic = p['total ion current']
                    t = p['scan start time']
                    ce = p['collision energy']
                    target = p['isolation window target m/z']
                    lowmz = p['scan window lower limit']
                    highmz = p['scan window upper limit']
                    
                    if limits.has_key(target) is False:
                        limits[target] = [lowmz,highmz]
                    x,y = self.binarytolist(spectrum)
                    if msms.has_key(target) is False:
                        msms[target] = {}
                    msms[target][t] = {'CE':ce,'TIC':tic,'x':list(x),'y':list(y)}
        if self.v is True:
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
        for spectrumList in self.tree.getElementsByTagName('spectrumList'): #for each spectrum list
            for spectrum in spectrumList.getElementsByTagName('spectrum'): # go through each spectrum
                attr = self.attributes(spectrum)
                p = self.cvparam(spectrum)
                if self.v is True:
                    self.sys.stdout.write('\rExctracting UV spectrum #%i/%i  %.1f%%' %(attr['index']+1,self.nscans,float(attr['index']+1)/float(self.nscans)*100.))
                mode,level = self.scantype(p)
                if mode == 'UV': # if type is UV-Vis
                    rtime.append(p['scan start time'])
                    x,y = self.binarytolist(spectrum)
                    if uvlambda is None: # if wavelength region has not yet been defined
                        uvlambda = list(x)
                    for ind,val in enumerate(y): # normalize y value by 1 million to bring value into a.u.
                        y[ind] = val/1000000.
                    uvint.append(y) # append intensity values
        if self.v is True:
            self.sys.stdout.write(' DONE\n')
        return rtime,uvlambda,uvint    

    def pwconvert(self,filename):
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
            import os
            locations = []
            for root,dirs,files in os.walk(path):
                if fname in files:
                    locations.append(os.path.join(root,fname))                   
            return locations
        
        import subprocess
        if self.sys.platform != 'win32':
            raise OSError('The conversion function of the mzML class is limited to Windows operating systems.\nYou can attempt to manually convert to *.mzML using the proteowizard standalone package (32-bit binary encoding precision)')
        locs = []
        for val in ['c:\\program files\\proteowizard','c:\\program files (x86)\\proteowizard']: #searches for msconvert.exe in expected folders
            locs.extend(find_all('msconvert.exe',val))
                    
        if len(locs)==0: #exits if script cannot find msconvert.exe
            raise IOError('The python script could not find msconvert.exe\nPlease ensure that ProteoWizard is installed in either:\nc:\\program files\\proteowizard\nor\nc:\\program files (x86)\\proteowizard')
        
        if self.v is True:
            self.sys.stdout.write('Generating *.mzML file from *.raw...')
            self.sys.stdout.flush()
            subprocess.call(locs[-1]+' "'+filename+'" --mzML --32 -v')
            self.sys.stdout.write(' DONE\n')
            self.sys.stdout.flush()  
        else:
            subprocess.call(locs[-1]+' "'+filename+'" --mzML --32')
    
    def scantype(self,hand):
        """determines the scan type of the provided spectrum"""
        if type(hand) == dict: # handed a parameters dictionary
            p = hand
        else: # handed a tree or branch
            p = self.cvparam(hand)
        
        if p.has_key('MS1 spectrum'): # normal MS spectrum
            MS = True
        if p.has_key('MSn spectrum'): # MSMS spectrum
            MS = True
        if MS is True: # if MS, determine level and mode
            level = int(p['ms level'])
            if p.has_key('negative scan'):
                mode = '-'
            if p.has_key('positive scan'):
                mode = '+'
            return mode,level
        elif p.has_key('electromagnetic radiation spectrum'): # otherwise is it a UV spectrum
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
    
    def sumscans(self,sr=None,mzrange=[50.,2000.],dec=3):
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
        for spectrumList in self.tree.getElementsByTagName('spectrumList'): #for each spectrum list
            for spectrum in spectrumList.getElementsByTagName('spectrum'): # go through each spectrum
                attr = self.attributes(spectrum) # get attributes
                p = self.cvparam(spectrum) # pull parameters of the scan
                mode,level = self.scantype(p)
                if self.v is True:
                    self.sys.stdout.write('\rCombining mass spectrum #%i/%i  %.1f%%' %(attr['index']+1,self.nscans,float(attr['index']+1)/float(self.nscans)*100.))
                if mode in ['+','-'] and level <2: # if type is full scan mass spectrum
                    if attr['index']+1 >= sr[0] and attr['index']+1 <= sr[1]:
                        x,y = self.binarytolist(spectrum)
                        if mzrange is not None: # if m/z range is specified, trim spectrum to specified range
                            x,y = self.trimspectrum(x,y,mzrange[0],mzrange[1])
                        spec.addspectrum(x,y)
        out = spec.trim()
        if self.v is True:
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

if __name__ == '__main__':
    filename = 'LY-2015-09-15 06'
    mzml = mzML(filename,verbose=True)
    #from _Spectrum import Spectrum
    #sp = {
    #'pos':{'bounds':[325,327],'affin':'+','spectrum':Spectrum(3),'raw':[]},
    #'neg':{'bounds':[348,350],'affin':'-','spectrum':Spectrum(3),'raw':[]},
    #'uv':{'bounds':[378,None],'affin':'UV','raw':[]}
    #}