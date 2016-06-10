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

to add:
    remove pwconvert and filepresent from tome (or generalize filepresent)
    update pullspeciesdata docstring to better represent the newest version of the supplied dictionary
"""


class mzML(object):
    def __init__(self,filename,verbose=True):
        """
        Class for interpreting an mzML (mass spectrum) file
        """
        self.v = verbose
        self.sys = __import__('sys')
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
        self.tree = xml.dom.minidom.parse(self.filename) # full mzML file
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
    
    def checkforfile(self,fn):
        """checks for file and converts if necessary"""
        if fn.lower().endswith('.mzml') is False:
            if fn.endswith('.raw') is True:
                self.filepresent(fn,ty='dir')
                if self.filepresent(fn[:-4]+'.mzML') is False:
                    self.pwconvert(fn) # generate mzml if false
                return fn[:-4]+'.mzML'
            else:
                fn = self.fixextension(fn) # try to fix extension
            
        else:
            self.filepresent(fn)
            return fn
   
    def decode(self,line,extension='.mzML'):
        """
        Function to decode the base64 strings contained in XML files
        Currently only functions with 32-bit precision
        
        Module requirements: base64, struct
        """
        decoded = self.b64.decodestring(line) #decodes base64 string 
        tmp_size = len(decoded)/4 #number of values in string
        if extension == ".mzML": # mzML files have separated string pairs
            unpack_format = "<%dL" % tmp_size # little-endian, number of values, unsigned long
            values = [] 
            for tmp in self.st.unpack(unpack_format,decoded):
                tmp_i = self.st.pack("I",tmp) #pack as unsigned integer I
                tmp_f = self.st.unpack("f",tmp_i)[0] #unpack as float f
                values.append(float(tmp_f))
            return values
    
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
                if self.filepresent(fn+oopsx[key],ty='file') is True:
                    return fn+oopsx[key]
        for key in oopsr: # tries to complete raw shortenings
            if fn.lower().endswith(key) is True:
                if self.filepresent(fn+oopsr[key],ty='dir') is True:
                    return fn+oopsr[key]
        if self.filepresent(fn+'.mzml') is True: # tries to add the extension
            return fn+'.mzml'
        raise IOError('The file "%s" could not be located in the current working directory'%(fn)) # if it can't be found, raise IOError
    
    def foreachscan(self):
        """
        test function, no particular use at this time
        """
        n = 0
        for spectrumList in self.tree.getElementsByTagName('spectrumList'):
            for spectrum in spectrumList.getElementsByTagName('spectrum'):
                n+=1
        print n
    
    def getText(self,nodelist):
        """
        function to get text from a simple XML object
        snippet written by Dennis Hore
        """
        rc = []
        for node in nodelist:
            if node.nodeType == node.TEXT_NODE:
                rc.append(node.data)
        return ''.join(rc)

    def integ(self,start,end,x,y):
        """
        Function to integrate values between paired list indicies (e.g. a m/z list and an intensity list)
        
        input (start,end,x,y)
        start:
            start x value
        end:
            end x value
        x:
            list of x values
        y:
            list of y values (paired with x)
        """
        if start < min(x):
            raise ValueError('\nThe provided start integration value {} is less than the minimum value of the spectrum being summed {}.\nCheck your input.'.format(start,min(x)))
        if end is None: # if there a range is not supplied, return a single value
            return y[self.takeclosest(x,start)]
        if end is not None:
            if start > max(x):
                raise ValueError('\nThe provided start integration value {} is greater than the maximum value of the spectrum being summed {}.\nCheck your input.'.format(start,max(x)))
            elif end > max(x):
                raise ValueError('\nThe provided end integration value {} is greater than the maximum value of the spectrum being summed {}.\nCheck your input.'.format(end,max(x)))
            return sum(y[self.br(x,start):self.bl(x,end)]) #take nearest value inside the bounds    
    
    def pullchromdata(self):
        """
        Pulls mzML chromatograms
        """
        chroms = {} #dictionary of chromatograms
        for chromatogramList in self.tree.getElementsByTagName('chromatogramList'):
            nchroms = int(chromatogramList.getAttribute('count')) # pull total number of scans
            for chromatogram in chromatogramList.getElementsByTagName('chromatogram'):
                curchrom = int(chromatogram.getAttribute('index'))+1 #get current chromatogram number
                if self.v is True:
                    self.sys.stdout.write('\rExtracting chromatogram #%i/%i  %.1f%%' %(curchrom,nchroms,float(curchrom)/float(nchroms)*100.))
                    self.sys.stdout.flush()
                units = []
                B64 = []
                for binaryList in chromatogram.getElementsByTagName('binaryDataArrayList'): #pull array data for mz and intensity
                    for cvParam in binaryList.getElementsByTagName('cvParam'): #find units
                        if len(cvParam.getAttribute('unitName')) >0:
                            units.append(str(cvParam.getAttribute('unitName'))) 
                    for binary in binaryList.getElementsByTagName('binaryDataArray'): # pull both binary strings in each chromatogram 
                        B64.append(self.getText(binary.getElementsByTagName('binary')[0].childNodes))
                    x,y = self.decode(B64[0]),self.decode(B64[1]) #decode binary data to a list of values
                    chroms[str(chromatogram.getAttribute('id'))] = {'x':x,'y':y,'xunit':units[0],'yunit':units[1]} #add dictionary entry for each chromatogram
        if self.v is True:
            self.sys.stdout.write(' DONE\n')
        return chroms

    def pullspeciesdata(self,sp):
        """
        Iterates through the loaded mzML file and extracts data for each of the provided species in the dictionary
        
        input:
            dictionary of selected peaks
            format {'peak name':{'bounds':[peak start mz,peak end mz],'affin':['+' or '-' or 'UV'},'raw':[],'spectrum':[]}
        output:
            filled dictionary, dictionary of TICs, dictionary of retention times
        
        explicitly interprets full scan mass spectra and UV species
        """
        rtime = {} # generate empty lists required for data processing
        TIC = {}
        for spectrumList in self.tree.getElementsByTagName('spectrumList'): #for each spectrum list
            nscans = int(spectrumList.getAttribute('count')) # pull total number of scans
            mode = None
            for spectrum in spectrumList.getElementsByTagName('spectrum'): # go through each spectrum
                curspec = int(spectrum.getAttribute('index'))+1 #get current spectrum number
                if self.v is True:
                    self.sys.stdout.write('\rExtracting species data from spectrum #%i/%i  %.1f%%' %(curspec,nscans,float(curspec)/float(nscans)*100.))
                mode,level = self.scantype(spectrum)
                if mode is not None and level < 2:
                    modekey = 'raw'+mode # define dictionary key for current scan
                    if rtime.has_key(modekey) is False: # create dictionary entry if not present
                        rtime[modekey] = []
                    if TIC.has_key(modekey) is False:
                        TIC[modekey] = []
                
                    for cvParam in spectrum.getElementsByTagName('cvParam'): #find TIC and time
                        if cvParam.getAttribute('name') == 'total ion current':
                            try:
                                TIC[modekey].append(int(cvParam.getAttribute('value')))
                            except ValueError:
                                TIC[modekey].append(float(cvParam.getAttribute('value')))
                        if cvParam.getAttribute('name') == 'scan start time':
                            rtime[modekey].append(float(cvParam.getAttribute('value')))
                
                    B64 = []
                    for binary in spectrum.getElementsByTagName('binaryDataArray'): #pull array data for mz and intensity
                        B64.append(self.getText(binary.getElementsByTagName('binary')[0].childNodes))
                    x,y = self.decode(B64[0]),self.decode(B64[1]) #decode binary data to a list of values
                    
                    for key in sp: # integrate each peak
                        if sp[key]['affin'] == mode: # if species has affinity to this spectrum type
                            if mode in ['+','-']: # if mass spectrum
                                sp[key]['raw'].append(self.integ(sp[key]['bounds'][0],sp[key]['bounds'][1],x,y))
                                # change this to bisect and .addspectrum()
                                for ind,mz in enumerate(x): # for each mz
                                    if mz >= sp[key]['bounds'][0] and mz <= sp[key]['bounds'][1]: # skips mz not in desired range
                                        sp[key]['spectrum'].addvalue(mz,y[ind])
                            if mode in ['UV']: # if UV spectrum
                                sp[key]['raw'].append(self.integ(sp[key]['bounds'][0],sp[key]['bounds'][1],x,y)/1000000.) # integrates and divides by 1 million bring it into au
                

        for key in sp: #remove mz/int values in spectrum that are None
            if sp[key]['affin'] in ['+','-']:
                sp[key]['spectrum'] = sp[key]['spectrum'].trim() # trim spectrum to remove Nonetypes        
        self.TIC = TIC
        self.rtime = rtime
        if self.v is True:
            self.sys.stdout.write(' DONE\n')
        return sp,TIC,rtime    
    
    def pullspectra(self,sr='all',mzrange=None):
        """
        iterates through a mzML tree and pulls full scan mass spectra

        sr:
            scan range to sum (default 'all')
        mzrange:
            m/z range to keep (None keeps the entire range)
            this is a memory saving utility
        
        output:
            dictionary of timpoints with subkeys for m/z, intensity, and scan #
            also outputs the updated scan range and mz range
        
        exclusively pulls mass spectra with MS level 1 (for MSMS spectra use pullmsmsspectra)
        """
        speclist = {} # list of dictionaries
        for spectrumList in self.tree.getElementsByTagName('spectrumList'): #for each spectrum list
            nscans = int(spectrumList.getAttribute('count')) # pull total number of scans
            if sr == 'all':
                sr = [1,nscans]
            for spectrum in spectrumList.getElementsByTagName('spectrum'): # go through each spectrum
                mode,level = self.scantype(spectrum)
                curspec = int(spectrum.getAttribute('index'))+1 #get current spectrum number
                if self.v is True:
                    self.sys.stdout.write('\rExtracting mass spectrum #%i/%i  %.1f%%' %(curspec,nscans,float(curspec)/float(nscans)*100.))
                if mode in ['+','-'] and level <2: # if type is full scan mass spectrum
                    if curspec >= sr[0] and curspec <= sr[1]:
                        for cvParam in spectrum.getElementsByTagName('cvParam'):
                            if cvParam.getAttribute('name') == 'scan start time': # scan time
                                t = float(cvParam.getAttribute('value'))
                        B64 = []
                        for binary in spectrum.getElementsByTagName('binaryDataArray'): #pull array data for mz and intensity
                            B64.append(self.getText(binary.getElementsByTagName('binary')[0].childNodes))
                        x,y = self.decode(B64[0]),self.decode(B64[1]) #decode binary data to a list of values
                        if mzrange is not None: # if m/z range is specified, trim spectrum to specified range
                            l,r = self.bl(x,mzrange[0]),self.br(x,mzrange[1])
                            x = x[l:r]
                            y = y[l:r]
                        speclist[t] = {'x':x,'y':y,'scan':curspec}
        if mzrange is None: # determine m/z range if not specified
            minx = 1000000.
            maxx = 0.
            for scan in speclist:
                if min(speclist[scan]['x']) < minx:
                    minx = min(speclist[scan]['x'])
                if max(speclist[scan]['x']) > minx:
                    maxx = max(speclist[scan]['x'])
            mzrange = [minx,maxx]
        if self.v is True:
            self.sys.stdout.write(' DONE\n')
        return speclist,sr,mzrange
    
    def pullmsmsspectra(self):
        """
        Iterates through a mzML tree and extracts any MSMS spectra
        
        Exclusively finds mass spectra with level greater or equal to 2
        groups spectra by key defined by thge isolation window target
        
        returns a dictionary containing all selected ions
        each dictionary has subkeys for each time point which each has subkeys containing the information of that scan
        
        also returns a dictionary of m/z limits using the same keys as msms
        """
        msms = {}
        limits = {}
        for spectrumList in self.tree.getElementsByTagName('spectrumList'): #for each spectrum list
            nscans = int(spectrumList.getAttribute('count')) # pull total number of scans
            mode = None
            for spectrum in spectrumList.getElementsByTagName('spectrum'): # go through each spectrum
                curspec = int(spectrum.getAttribute('index'))+1 #get current spectrum number
                if self.v is True:
                    self.sys.stdout.write('\rExtracting species data from spectrum #%i/%i  %.1f%%' %(curspec,nscans,float(curspec)/float(nscans)*100.))
                mode,level = self.scantype(spectrum)
                if level >= 2: # if it is a msms spectrum
                    for cvParam in spectrum.getElementsByTagName('cvParam'): #find TIC and time
                        if cvParam.getAttribute('name') == 'total ion current':
                            try:
                                tic = int(cvParam.getAttribute('value'))
                            except ValueError:
                                tic = float(cvParam.getAttribute('value'))
                        if cvParam.getAttribute('name') == 'scan start time':
                            t = float(cvParam.getAttribute('value'))
                        if cvParam.getAttribute('name') == 'collision energy': # collision energy applied
                            ce = float(cvParam.getAttribute('value'))
                        if cvParam.getAttribute('name') == 'isolation window target m/z': # target m/z value
                            target = cvParam.getAttribute('value')
                        if cvParam.getAttribute('name') == 'scan window lower limit': # target m/z value
                            lowmz = float(cvParam.getAttribute('value'))
                        if cvParam.getAttribute('name') == 'scan window upper limit': # target m/z value
                            highmz = float(cvParam.getAttribute('value'))
                    if limits.has_key(target) is False:
                        limits[target] = [lowmz,highmz]
                    B64 = []
                    for binary in spectrum.getElementsByTagName('binaryDataArray'): #pull array data for mz and intensity
                        B64.append(self.getText(binary.getElementsByTagName('binary')[0].childNodes))
                    x,y = self.decode(B64[0]),self.decode(B64[1]) #decode binary data to a list of values
                    if msms.has_key(target) is False:
                        msms[target] = {}
                    msms[target][t] = {'CE':ce,'TIC':tic,'x':list(x),'y':list(y)}
        return msms,limits
    
    def pullUVspectra(self):
        """
        Iterates through a mzML tree and extracts UV-Vis spectra
        
        returns a list of timepoints, a list of wavelengths, and a list of lists of intensity corresponding to timepoints
        """
        rtime = [] # generate empty lists required for data processing
        uvlambda = None
        uvint = []
        for spectrumList in self.tree.getElementsByTagName('spectrumList'): #for each spectrum list
            nscans = int(spectrumList.getAttribute('count')) # pull total number of scans
            for spectrum in spectrumList.getElementsByTagName('spectrum'): # go through each spectrum
                curspec = int(spectrum.getAttribute('index'))+1 #get current spectrum number
                if self.v is True:
                    self.sys.stdout.write('\rExctracting UV spectrum #%i/%i  %.1f%%' %(curspec,nscans,float(curspec)/float(nscans)*100.))
                mode,level = self.scantype(spectrum)
                if mode == 'UV': # if type is UV-Vis
                    for cvParam in spectrum.getElementsByTagName('cvParam'): #find time
                        if cvParam.getAttribute('name') == 'scan start time':
                            rtime.append(float(cvParam.getAttribute('value')))
                    B64 = []
                    for binary in spectrum.getElementsByTagName('binaryDataArray'): #pull array data for wavelength and intensity
                        B64.append(self.getText(binary.getElementsByTagName('binary')[0].childNodes))
                    x,y = self.decode(B64[0]),self.decode(B64[1]) #decode binary data to a list of values
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
    
    
    def scantype(self,spectrum):
        """
        determines the scan type of the provided spectrum
        """
        MS = False
        level = 0
        UV = False
        mode = None
        for cvParam in spectrum.getElementsByTagName('cvParam'): # determine if spec or uv vis and set affinity
            if cvParam.getAttribute('name') == 'MS1 spectrum':
                MS = True
            if cvParam.getAttribute('name') == 'MSn spectrum':
                MS = True
            if cvParam.getAttribute('name') == 'ms level':
                level = int(cvParam.getAttribute('value'))
            if cvParam.getAttribute('name') == 'electromagnetic radiation spectrum':
                UV = True
                break
            if cvParam.getAttribute('name') == 'negative scan':
                mode = '-'
            if cvParam.getAttribute('name') == 'positive scan':
                mode = '+'
        if UV is True: # electromagnetic radiation spectrum
            return 'UV',level
        elif MS is True: # MS scan
            return mode,level
        #elif MS is True and level is 1: # full scan
        #    return mode,level
        #elif MS is True and level > 1: # MSMS of some type
        #    return mode,level
        else: # type that is not interpretable
            return None,None        

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

if __name__ == '__main__':
    filename = 'EJ-UVTQ-056-MRM-NL.m'
    mzml = mzML(filename,verbose=True)
    from _Spectrum import Spectrum
    sp = {
    'pos':{'bounds':[325,327],'affin':'+','spectrum':Spectrum(3),'raw':[]},
    'neg':{'bounds':[348,350],'affin':'-','spectrum':Spectrum(3),'raw':[]},
    'uv':{'bounds':[378,None],'affin':'UV','raw':[]}
    }