"""
Spectrum class for combining spectra with different dimensions but are in the same domain

CHANGELOG
new:
    ---2.3---
    fixed addition (arrays weren't working with addspectrum)
    trim method can now be restricted to output between specific x bounds (allows for easy narrowing of range for extracting segments from the parent spectrum)
    ---2.4

to add/fix:
    update subtract method (reference add for what to change)
"""

#from _ScriptTime import ScriptTime
#st = ScriptTime(profile=True)
    

class Spectrum(object):
    def __init__(self,decpl,startmz=50.,endmz=2000.,specin=None):
        """
        class for generating a full spectrum with specified number of decimal places and None values for intensity
        used for summing mass spectra where dimensions are unequal
        decpl: (int)
            number of decimal places to track
        startmz: (float)
            start m/z of list
            default 50
        endmz: (float)
            end m/z of list
            default 2000
        specin: [[x values],[y values]]
            can supply an original spectrum on initialization
            list of lists of index-matched x and y values
        """
        self.decpl = decpl
        self.startmz = round(startmz,self.decpl)
        self.endmz = round(endmz,self.decpl)
        self.sp = __import__('scipy')
        self.x,self.yimm = self.fullspeclist(self.startmz,self.endmz) # m/z and intensity lists (these are intended to remain immutable)
        self.y = list(self.yimm) # create the list that will be actively modified
        if specin is not None:
            self.addspectrum(specin[0],specin[1])
    
    def __str__(self):
        return 'Full spectrum list from {} to {} keeping {} decimal places'.format(self.startmz,self.endmz,self.decpl)
    
    def __repr__(self):
        return '{}(decpl={},startmz={},endmz={})'.format(self.__class__.__name__,self.decpl,self.startmz,self.endmz)
    
    def __len__(self):
        return len(self.x)
    
    def __getitem__(self,ind):
        return [self.x[ind],self.y[ind]]
    
    def __add__(self,x):
        """
        Since addition to this class requires generating a complete copy of the class then addition,
        using the built in addition methods is recommended
        e.g. to add a single value use .addvalue()
        to add a spectrum use .addspectrum()
        """
        if isinstance(x,self.__class__) is True: # if it is another Spectrum instance
            if x.decpl != self.decpl:
                raise ValueError('The decimal places of the two spectra to be added are not equal. Addition is not supported')
            newstart = min(min(self.x),min(x.x)) # find new start m/z
            newend = max(max(self.x),max(x.x)) # find new end m/z
            tempnsp = Spectrum(self.decpl,newstart,newend) # temporary instance
            specin = self.trim()
            tempnsp.addspectrum(specin[0],specin[1]) # add self spectrum
            specin = x.trim()
            tempnsp.addspectrum(specin[0],specin[1]) # add input spectrum
            return tempnsp
        elif type(x) is int: # add this integer to every m/z
            tempnsp = Spectrum(self.decpl,self.startmz,self.endmz)
            tempnsp.addspectrum(self.x,self.y)
            for mz in self.x:
                tempnsp.addvalue(mz,x)
            return tempnsp
        elif len(x) == 2 and len(x[0]) == len(x[1]): # if it is a list of paired lists (another spectrum)
            tempnsp = Spectrum(self.decpl,self.startmz,self.endmz)
            specin = self.trim()
            tempnsp.addspectrum(specin[0],specin[1])
            tempnsp.addspectrum(x[0],x[1])
            return tempnsp
        else:
            return 'Addition of %s to the Spectrum class is unsupported' %`x`
    
    def __sub__(self,x):
        if isinstance(x,self.__class__) is True: # if it is another Spectrum instance
            if x.decpl != self.decpl:
                raise ValueError('The decimal places of the two spectra to be added are not equal. Subtraction is not supported')
            newstart = min(min(self.x),min(x.x)) # find new start m/z
            newend = max(max(self.x),max(x.x)) # find new end m/z
            tempnsp = Spectrum(self.decpl,newstart,newend) # temporary instance
            for ind,mz in enumerate(self.x): # add original list
                tempnsp.addvalue(mz,self.y[ind])
            for ind,mz in enumerate(x.x): # subtract new values
                tempnsp.addvalue(mz,-x.y[ind])
            return tempnsp
        elif type(x) is int: # add this integer to every m/z
            tempnsp = Spectrum(self.decpl,self.startmz,self.endmz)
            for mz in self.x:
                tempnsp.addvalue(mz,-x)
            return tempnsp
        else:
            return 'Subtraction of %s from the Spectrum class is unsupported' %`x`
    
    def __mul__(self,x):
        return 'Multiplication of the Spectrum class is unsupported'
    def __div__(self,x):
        return 'Division of the Spectrum class is unsupported'    
    
    def addvalue(self,xval,yval):
        """
        adds an intensity value to the mutable y list
        """
        try:
            index = self.index(xval)
            try:
                self.y[index] += yval # try to add value
            except TypeError:
                self.y[index] = yval # if None, then set to value
        except ValueError: # if index is not in spectrum
            pass # do nothing (the value will not be added to the spectrum)
        #if self.y[index] < 0: # catch for negative intensities while subtracting
        #    return ValueError('The intensity value for m/z %f is now negative (%f)'%(xval,self.y[index]))
    
    def addspectrum(self,x,y,subtract=False):
        """
        adds an entire spectrum to the Spectrum object
        assumes that x and y are paired values
        (the x values do not need to be sorted)
        
        subtract tells the method to subtract values instead of adding them (useful for comparing spectra)
        """
        if len(x) != len(y):
            raise ValueError('The addspectrum() method only supports two lists of the same dimension')
        if subtract is True: # if subtraction is called for
            for ind,mz in enumerate(x):
                self.addvalue(mz,-y[ind])
        else:
            for ind,mz in enumerate(x):
                self.addvalue(mz,y[ind])
    
    def checknone(self):
        """counts the number of not-None values in the current y list"""
        count = 0
        for val in self.y:
            if val is not None:
                count += 1
        return count
    
    def cp(self):
        """
        returns a list (clone) of the empty spectrum
        """
        return [list(self.x),list(self.yimm)]
    
    def cpfilled(self):
        """
        returns a list (clone) of the filled spectrum
        """
        return [list(self.x),list(self.y)]
    
    def fullspeclist(self,startmz,endmz):
        """
        Generates two paired lists (one m/z, one None) from start to end with a specified number of decimal places
        """
        x = self.sp.arange(startmz,endmz+10**-self.decpl,10**-self.decpl) # generate x values using arange and convert to list
        # endmz + increment ensures that the end value is present in the list
        y = [None]*len(x) # generate y list of equal length
        return x,y
    
    #@st.profilefn
    def fillzeros(self):
        """takes the current spectrum and replaces None with zeros"""
        for ind,inten in enumerate(self.y):
            if inten is None:
                self.y[ind] = 0
        return self.y
    
    #@st.profilefn    
    def index(self,mzval):
        """
        Calculates index of a given mz value in the object's list
        mzval: m/z value to find index of
            float
        """
        if mzval > self.endmz or mzval < self.startmz:
            raise ValueError('the m/z value ({}) is outside of the m/z range of this spectrum ({}-{})'.format(mzval,self.startmz,self.endmz))
        else:
            # uses array search; where() does not account for differences arising from floating point rounding
            # searchsorted is slightly less efficient than calculating the index (~100 ns slower)
            return self.sp.searchsorted(self.x,round(mzval,self.decpl)-10**-self.decpl)
            
            ## index can be instead calculated for small performance upgrade (~ 1 order of magnitude faster on average)
            ## calculates mismatches ~4.5% of the time, with maximum differences placing the value in the wrong index
            #return int(round((mzval-self.startmz)*(10**self.decpl))) # rounds after multiplication
    
    def normalize(self,top=100.):
        """normalizes the spectrum to the specified value"""
        m = max(self.y)
        for ind,inten in enumerate(self.y):
            if inten is not None:
                self.y[ind] = inten/m*top
    
    #@st.profilefn        
    def resety(self):
        """
        resets the y list to None
        """
        self.y = list(self.yimm)
        return 'intensity list was reset'
    
    def sum(self):
        """sums the y values"""
        out = 0
        for val in self.y:
            if val is not None:
                out += val
        return out
    
    def threshold(self,thresh):
        """
        trims the spectrum to a particular threshold value
        all intensity values below this threshold will be dropped
        """
        for ind,inten in enumerate(self.y):
            if inten < thresh:
                self.y[ind] = None
    
    #@st.profilefn
    def trim(self,zeros=False,xbounds=None):
        """
        trims pairs that have None intensity
        zeros specifies whether there should be zeros at the startmz and endmz (for generating continuous spectra across the range)
        the zeros will not overwrite existing intensity at that m/z
        """
        if xbounds is None:
            xbounds = [self.startmz,self.endmz]
        xout = []
        yout = []
        for ind,inten in enumerate(self.y):
            if self.x[ind] >= xbounds[0] and self.x[ind] <= xbounds[1]: # if within the x bounds
                if inten is not None:
                    xout.append(round(self.x[ind],self.decpl)) # rounded to avoid array floating point weirdness
                    yout.append(inten)
                elif zeros is True: # if zeros at the edges of spectrum are desired
                    if self.x[ind] == xbounds[0] or self.x[ind] == xbounds[1]: # at the edges of the output spectrum
                        xout.append(self.x[ind])
                        yout.append(0)
        return [xout,yout]
    
       
if __name__ == '__main__':
    nsp = Spectrum(3)
    print nsp
    
    #for i in range(10):
    #    nsp.fillzeros()
    #    nsp.trim(zeros=True)
    #    nsp.resety()
    #
    #st.printprofiles()
    