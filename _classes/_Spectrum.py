"""
This class is designed to efficiently combine, add to, or otherwise manipulate
spectra together whose dimensions are not equal.
For example, combining mass spectra together where resolution to the 3rd 
decimal (not to the 10th decimal) is desired.
Upon initialization, specify the number of decimal places desired. 
Start and end values for the x bounds may also be specified, and
an input spectrum can be provided (this spectrum will be added to 
the object on initialization).
When adding a value to the Spectrum object, it will find the closest x value
with the decimal place specified and add the y value to that x in the object.
e.g. if the decimal place is 3, adding x=545.34898627,y=10 will add 10 to x=545.349

Once the desired spectrum has been constructed, calling Spectrum.trim() will return
an [[x values],[y values]] list with only the x values that have intensities. Other
manipulations are available, see below for details.

CHANGELOG
new:
    ---2.3---
    fixed addition (arrays weren't working with addspectrum)
    trim method can now be restricted to output between specific x bounds (allows for easy narrowing of range for extracting segments from the parent spectrum)
    addvalue now ignores nonetypes instead of overwriting existing values
    add and subtract methods are updated to be more efficient (no longer need to trim)
    elaborated get to accept float values
    mul, div, and pow attempts now raise an AttributeError
    replaced startmz and endmz with start and end
    ---2.4
"""

class Spectrum(object):
    def __init__(self,decpl,start=50.,end=2000.,specin=None):
        """
        A class for manipulating a spectrum with the specified number of decimal places
        
        decpl: (int) REQUIRED
            number of decimal places to track
        start: (float)
            start m/z of list
            default 50
        end: (float)
            end m/z of list
            default 2000
        specin: [[x values],[y values]]
            can supply an original spectrum on initialization
            list of lists of index-matched x and y values
        """
        self.decpl = decpl
        self.start = round(start,self.decpl)
        self.end = round(end,self.decpl)
        self.sp = __import__('scipy')
        self.x,self.yimm = self.fullspeclist(self.start,self.end) # m/z and intensity lists (these are intended to remain immutable)
        self.y = list(self.yimm) # create the list that will be actively modified
        if specin is not None:
            self.addspectrum(specin[0],specin[1])
    
    def __str__(self):
        return 'Full spectrum list from {} to {} keeping {} decimal places'.format(self.start,self.end,self.decpl)
    
    def __repr__(self):
        return '{}(decpl={},start={},end={})'.format(self.__class__.__name__,self.decpl,self.start,self.end)
    
    def __len__(self):
        return len(self.x)
    
    def __getitem__(self,ind):
        """
        if supplied index is an integer, return the x and y value of that index in the list
        if a float, return the intensity of that m/z
        """
        if type(ind) is int:
            return [self.x[ind],self.y[ind]]
        elif type(ind) is float: # returns the intensity value of the specified m/z
            if ind < self.start or ind > self.end:
                raise IndexError('The supplied float %f is outside of the m/z range of this Spectrum instance (%.3f -%.3f)' %(ind,self.start,self.end))
            return self.y[self.index(ind)]
    
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
            newstart = min(self.start, x.start) # find new start m/z
            newend = max(self.end, x.end) # find new end m/z
            tempnsp = Spectrum(self.decpl,newstart,newend,[self.x,self.y]) # temporary instance with self as specin
            tempnsp.addspectrum(x.x,x.y) # add incoming spectrum
            return tempnsp
        elif type(x) is int: # add this integer to every m/z
            tempnsp = Spectrum(self.decpl,self.start,self.end,[self.x,self.y])
            for mz in tempnsp.x:
                tempnsp.addvalue(mz,x)
            return tempnsp
        elif len(x) == 2 and len(x[0]) == len(x[1]): # if it is a list of paired lists (another spectrum)
            tempnsp = Spectrum(self.decpl,self.start,self.end,[self.x,self.y])
            tempnsp.addspectrum(x[0],x[1])
            return tempnsp
        else:
            return 'Addition of %s to the Spectrum class is unsupported' %`x`
    
    def __sub__(self,x):
        if isinstance(x,self.__class__) is True: # if it is another Spectrum instance
            if x.decpl != self.decpl:
                raise ValueError('The decimal places of the two spectra to be added are not equal. Subtraction is not supported')
            newstart = min(self.start, x.start) # find new start m/z
            newend = max(self.end, x.end) # find new end m/z
            tempnsp = Spectrum(self.decpl,newstart,newend,[self.x,self.y]) # temporary instance
            tempnsp.addspectrum(x.x,x.y,True) # subtract incoming spectrum
            return tempnsp
        elif type(x) is int: # add this integer to every m/z
            tempnsp = Spectrum(self.decpl,self.start,self.end,[self.x,self.y])
            for mz in self.x:
                tempnsp.addvalue(mz,-x)
            return tempnsp
        elif len(x) == 2 and len(x[0]) == len(x[1]): # if it is a list of paired lists (another spectrum)
            tempnsp = Spectrum(self.decpl,self.start,self.end,[self.x,self.y])
            tempnsp.addspectrum(x[0],x[1],True) # subtract the incoming spectrum
            return tempnsp
        else:
            return 'Subtraction of %s from the Spectrum class is unsupported' %`x`
    
    def __mul__(self,x):
        raise AttributeError('Multiplication of the Spectrum class is unsupported')
    def __div__(self,x):
        raise AttributeError('Division of the Spectrum class is unsupported') 
    def __pow__(self,x):
        raise AttributeError('Raising a Spectrum instance to a power is unsupported.\nAlso... really?!')
    
    def addvalue(self,xval,yval,subtract=False):
        """adds an intensity value to the mutable y list"""
        if yval is not None: # if handed an actual value
            try: # try indexing
                index = self.index(xval)
                if subtract is True: # set sign based on input
                    sign = -1
                else:
                    sign = 1
                try:
                    self.y[index] += yval*sign # try to add value
                except TypeError:
                    self.y[index] = yval*sign # if None, then set to value
            except ValueError: # if index is not in spectrum
                pass # do nothing (the value will not be added to the spectrum)
    
    def addspectrum(self,x,y,subtract=False):
        """
        adds an entire spectrum to the Spectrum object
        assumes that x and y are paired values
        (the x values do not need to be sorted)
        
        subtract tells the method to subtract values instead of adding them (useful for comparing spectra)
        """
        if len(x) != len(y):
            raise ValueError('The addspectrum() method only supports two lists of the same dimension')
        for ind,mz in enumerate(x):
            self.addvalue(mz,y[ind],subtract)
    
    def checknone(self):
        """counts the number of not-None values in the current y list"""
        count = 0
        for val in self.y:
            if val is not None:
                count += 1
        return count
    
    def cp(self):
        """returns a list (clone) of the empty spectrum"""
        return [list(self.x),list(self.yimm)]
    
    def cpfilled(self):
        """returns a list (clone) of the filled spectrum"""
        return [list(self.x),list(self.y)]
    
    def fullspeclist(self,start,end):
        """
        Generates two paired lists (one m/z, one None) from start to end with a specified number of decimal places
        """
        x = self.sp.arange(start,end+10**-self.decpl,10**-self.decpl) # generate x values using arange and convert to list
        # end + increment ensures that the end value is present in the list
        y = [None]*len(x) # generate y list of equal length
        return x,y
    
    def fillzeros(self):
        """takes the current spectrum and replaces None with zeros"""
        for ind,inten in enumerate(self.y):
            if inten is None:
                self.y[ind] = 0
        return self.y
    
    def index(self,mzval,method='search'):
        """
        Calculates index of a given mz value in the object's list
        mzval: float
            m/z value to find index of
        method: 'search' or 'calculate'
            search uses an array search
            calculate attempts to calculate the index based on the start and end values
        """
        if mzval > self.end or mzval < self.start:
            raise ValueError('the m/z value ({}) is outside of the m/z range of this spectrum ({}-{})'.format(mzval,self.start,self.end))
        else:
            """
            details regarding the method:
            where() does not account for differences arising from floating point rounding
            searchsorted is slightly less efficient than calculating the index (~100 ns slower)
            
            calculation can be ~ 1 order of magnitude faster on average than searching, but
            calculates mismatches ~4.5% of the time, with maximum differences placing the value in the wrong index
            """
            if method == 'search': # uses array search; 
                return self.sp.searchsorted(self.x,round(mzval,self.decpl)-10**-self.decpl)
            elif method == 'calculate': # calculate the location
                return int(round((mzval-self.start)*(10**self.decpl))) # rounds after multiplication
    
    def normalize(self,top=100.):
        """normalizes the spectrum to the specified value"""
        m = max(self.y)
        for ind,inten in enumerate(self.y):
            if inten is not None:
                self.y[ind] = inten/m*top
    
    def resety(self):
        """
        resets the y list to None
        this is substantially faster than creating a new spectrum instance and is
        recommended if the same spectrum object is used repeatedly
        """
        self.y = list(self.yimm)
        return 'intensity list was reset'
    
    def sum(self):
        """returns the sum of all y values"""
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
    
    def trim(self,zeros=False,xbounds=None):
        """
        trims pairs that have None intensity
        zeros specifies whether there should be zeros at the start and end (for generating continuous spectra across the range)
        the zeros will not overwrite existing intensity at that m/z
        """
        if xbounds is None:
            xbounds = [self.start,self.end]
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
    spec = Spectrum(3)
    print spec