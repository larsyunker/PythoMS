"""
Molecule class (previously "isotope pattern generator" and "MolecularFormula")

The output of this builder has been validated against values calculated by ChemCalc (www.chemcalc.org)
Negligable differences are attributed to different low value discarding techniques
(ChemCalc keeps the top 5000 peaks, this script drops values less than a threshold 5 orders of magnitude below the maximum value)

CHANGELOG:
    ---0.1---BETA
    rebuilt rawisotopepattern to be more efficient and retain all decimals
    fixed error in calcindex (within tome_v02; it was returning an index of one less than intended)
    ---0.2---BETA
    rewrote barisotopepattern to group m/z values and determine mass by weighted average
    added support for addition and subtraction from the MF class
    added returns (of unsupport) for multiplication and division
    ---0.3---
    enhanced formula interpreter (now handles single brackets and common abbreviations)
    added support for common functional groups in _formabbrvs
    ---1.0---
    added support for specification of isotopes
    added support for multiplication and division by integers
    supports nested brackets (even of the same type)
    ---2.0---
    fixed gaussian addition by using NoneSpectrum instances
    removed generation of the gaussian isotope pattern from the automatic execution
    calling of plotgaus() and gaussianisotopepattern() generates the gaussian pattern
    fixed add,sub,mul,div to raise errors and made the errors more explicit
    added reset functionality to reset the instance to its original values from when the instance was called
    ---2.1---
    renamed to Molecule
    added molecular weight calculator
    added percent composition calculator
    fixed molecular formula generator to exclude implicit 1's
    tweaked instance checks in __add__ and __sub__ to look for the class name
    added a percent composition printer (formatted output)
    added amino acids to the abbreviations dictionary (can now interpret polypeptides)
    removed all sys.exit calls and changed them to error raising
    incorporated detail printing into the class as a callable function
    removed resolution from the class call and incorporated it specifically into the functions which require it
    updated dostrings
    ---2.3---
    reinstated resolution in the class call for use with PyRSIM
    creating bounds function for integration with PyRSIM
    shuffled fwhm and sigma calculations into a separate function and call them in __init__
    fixed gaussian isotope pattern generation to work with NoneSpectrum addition method (now works correctly)
    added compare function to compare an experimental pattern to the predicted one
    bounds no longer includes peaks below a certain threshold
    created dictionary from the crc handbook masses (nearly identical masses)
    ---2.4---
    significantly improved the generation of the raw spectrum (now uses the built in methods of the Spectrum class -- formerly the NoneSpectrum class)
    added verbose call to raw spectrum
    significantly improved the efficiency of gaussianisotopepattern by calling a single Spectrum object and adding the subsequent spectra to that
    centered bar isotope pattern
    added plotraw() in case a raw check is desired
    fixed isotope handling in molecularweight()
    ---2.5---
    added charge interpreter into both the charge input and the string input (the string input will override the charge input)
    added a catch for not closing a bracket
    exactmass() can account for the mass of an electron, but this is commented out because of the extremely minute difference this makes
    removed __pow__ method because it's just so farfetched that someone might use it
    changed print calls to sys.stdout.write calls
    corrected generation of the formula if the class is added/subtracted/multiplied/divided
    ---2.6 building

to add:
    
"""

#from _ScriptTime import ScriptTime
#st = ScriptTime(profile=True)

class Molecule(object):
    def __init__(self,string,charge=1,res=5000):
        """
        Determines many properties of a given molecule
        
        string: (string) the molecule to interpret
        charge: (int or string) the charge of the molecule (for mass spectrometric applications)
        res: (int) the resolution of the mass spectrometer
        
        supported input:
        common abbreviations can be predefined in _formabbrvs.py
        brackets to signify multiples of a given component (nested brackets are supported)
        specification of isotopes (e.g. carbon 13 could be specified as "(13C)" )
        charge can be specified either in the string by enclosing it in brackets (e.g. "(2+)" or in the charge kwarg)
        """
        from _nist_mass import nist_mass as mass # masses from the NIST database
        #from _crc_mass import crc_mass as mass # masses from the CRC Handbook of Chemistry and Physics
        self.md = mass # mass dictionary that the script will use
        self.formula = string # input formula
        self.charge,self.sign = self.interpretcharge(charge) # charge
        self.res = res # the resolution of the mass spectrometer
        self.comp = self.composition(self.formula) # determine composition from formula
        self.checkinnist(self.comp) # checks that all the composition keys are valid
        self.calculate()
        self.default()
    
    def __str__(self):
        return "Molecule {}".format(self.formula)
    
    def __repr__(self):
        return "{}('{}')".format(self.__class__.__name__,self.formula)
    
    def __add__(self,x):
        if type(x) is str:
            addition = self.composition(x)
        else:
            if isinstance(x,self.__class__) is True:
                addition = x.comp
            else:
                raise ValueError('Addition of {} to Molecule object {} is invalid'.format(x,self.formula))
        for key in addition:
            try:
                self.comp[key] += addition[key]
            except KeyError:
                self.comp[key] = addition[key]
        self.calculate() # recalculates mass and isotope patterns
        self.formula = self.sf
        return self.sf
    
    def __sub__(self,x):
        if type(x) is str:
            addition = self.composition(x)
        else:
            if isinstance(x,self.__class__) is True:
                addition = x.comp
            else:
                raise ValueError('Subtraction of {} from Molecule object {} is invalid'.format(x,self.formula))
        for key in addition:
            try:
                if self.comp[key] - addition[key] >0:
                    self.comp[key] -= addition[key]
                elif self.comp[key] - addition[key] is 0:
                    del self.comp[key]
                else:
                    raise ValueError('Subtracting %d number of element %s from %s would yield a negative amount.' %(addition[key],key,self.sf))
            except KeyError:
                return 'The molecular formula {} does not have an element in {}'.format(self.formula,x)
        self.calculate() # recalculates mass and isotope patterns
        self.formula = self.sf
        return self.sf
    
    def __mul__(self,x):
        if type(x) != int:
            raise ValueError('Non-integer multiplication of a Molecule object is unsupported')
        for key in self.comp:
            self.comp[key] = self.comp[key]*x
        self.calculate() # recalculates mass and isotope patterns
        self.formula = self.sf
        return self.sf
    
    def __div__(self,x):
        if type(x) != int:
            raise ValueError('Non-integer division of a Molecule object is unsupported')
        for key in self.comp:
            tempval  = float(self.comp[key])/float(x)
            if tempval.is_integer() is False:
                raise ValueError('Division of %s (%s) by %d yielded a non-integer number of element %s' %(self.formula,self.sf,x,key))
            else:
                self.comp[key] = int(tempval)
        self.calculate() # recalculates mass and isotope patterns
        self.formula = self.sf
        return self.sf
    
    def barisotopepattern(self,rawip,charge,dec=3):
        """
        generates an isotope pattern for use in bar plots
        effectively this consolidates all mass defects into a single peak determined from the exact mass
        """
        def groupmasses(ip,delta=0.5):
            """
            groups masses in an isotope pattern
            looks for differences in m/z greater than the specified delta
            expects a paired list of [[mz values],[intensity values]]
            """
            num = 0
            out = [[[],[]]]
            for ind,val in enumerate(ip[0]):
                out[num][0].append(ip[0][ind])
                out[num][1].append(ip[1][ind])
                try:
                    if ip[0][ind+1]-ip[0][ind] > delta:
                        num+=1
                        out.append([[],[]])
                except IndexError:
                    continue
            return out
        
        def ipgrouptopoint(ipgroup):
            """
            Takes a group of mz values and consolidates them into a single value
            determines point by weighted average of the m/z values
            """
            s = 0
            for ind,val in enumerate(ipgroup[0]): # sum mz*int pairs
                s+= val*ipgroup[1][ind]
            return s/sum(ipgroup[1]),sum(ipgroup[1]) # return weighted m/z, summed intensity
        
        groupedip = groupmasses(rawip)
        out = [[],[]]
        for group in groupedip:
            x,y = ipgrouptopoint(group) # determine weighted mass and summed intensity
            out[0].append(x)
            out[1].append(y)
        maxint = max(out[1])
        for ind,val in enumerate(out[1]): 
            out[0][ind] = out[0][ind]/abs(charge)
            out[1][ind] = val/maxint*100. # normalize to 100
        return out
    
    def bounds(self,conf=0.95,perpeak=False,threshold=1):
        """
        calculates bounds based on a set confidence interval and the bar isotope pattern
        for use with RSIM calculations
        conf: (float) the confidence interval to use
        perpeak: (bool) toggle for whether the function should return a dictionary of 
        boundaries for each peak, or a single pair of bounds that covers the entire isotope pattern
        threshold: (int) minimum threshold for peaks to be included in bounds
        """
        from scipy import stats
        tempip = [[],[]]
        for ind,inten in enumerate(self.barip[1]): # checks for intensities above threshold
            if inten >= threshold:
                tempip[0].append(self.barip[0][ind])
                tempip[1].append(self.barip[1][ind])
        if perpeak is True: # if per-peak bounds are called for
            out = {}
            for mz in tempip[0]:
                out[str(mz)] = {}
                out[str(mz)]['bounds'] = stats.norm.interval(conf,mz,scale=self.sigma)
        else: # a general range that covers the entire isotope pattern
            out = [None,None]
            out[0] = stats.norm.interval(conf,tempip[0][0],scale=self.sigma)[0]
            out[1] = stats.norm.interval(conf,tempip[0][-1],scale=self.sigma)[1]
        return out
        
    
    def calculate(self):
        """calls the calculation functions"""
        self.sf = self.molecularformula() # generates a string version of the molecular formula
        self.em = self.exactmass(self.comp,charge=self.charge) # monoisotopic mass (will not work for large number of carbons)
        self.fwhm,self.sigma = self.sigmafwhm(self.res,self.em)
        self.mw,self.pcomp = self.molecularweight() # molecular weight and elemental percent composition
        self.rawip = self.rawisotopepattern(self.comp,dec=2,verbose=False) # generates a raw isotope pattern (charge of 1)
        self.barip = self.barisotopepattern(self.rawip,self.charge) # bar isotope pattern based on the generated raw pattern
        #self.gausip = self.gaussianisotopepattern(self.barip,self.em,res=self.res) # simulated normal distribution of the bar isotope pattern
        ## gausip is rarely used and is a bit more time consuming, so it now requires a specific call. If call is desired on generation, uncomment the above line
    
    def checkinnist(self,comp):
        """checks for each 'element' in the nist dictionary and returns an error if not found"""
        for key in comp:
            if self.md.has_key(key) is False:
                ele,iso = self.isotope(key)
                if self.md[ele].has_key(iso) is False:
                    raise ValueError('The element "%s" does not have a defined isotope "%d" in the NIST element database, please check your input' %(ele,iso))
                    
    def compare(self,exp):
        """
        compares a provided real spectrum to the simulated gaussian isotope pattern
        exp: paired list of lists of mz and intensity values
        
        returns the standard error of the regression (lower is better)
        (this a measure of the average distance between the experimental and predicted lines)
        """
        def sumsquare(lst):
            """calculates the sum of squares"""
            ss = 0
            for val in lst:
                ss += val**2
            return ss
        
        if self.__dict__.has_key('gausip') is not True: # generate gaussian isotope pattern if not already generated
            self.gaussianisotopepattern()
        yvals = []
        res = []
        #tot = []
        maxy = float(max(exp[1]))
        for ind,val in enumerate(exp[1]): # normalize y values
            yvals.append(float(val)/maxy*100.)
        #avgy = sum(exp[1])/len(exp[1])
        for ind,mz in enumerate(exp[0]):
            if mz > min(self.gausip[0]) and mz < max(self.gausip[0]): # if within isotope pattern
                nspind = self.nsp.index(mz) # calculate index
                if self.nsp.y[nspind] is not None: # if the predicted intensity is not None
                    res.append(yvals[ind]-self.nsp.y[nspind]) # difference between observed and predited (residuals)
                    #tot.append(self.nsp.y[nspind]-avgy) # difference between predicted and mean
        #rsqrd = 1-(sumsquare(res)/sumsquare(tot)) # r-squared value (apparently not applicable to non-linear fits)
        from math import sqrt
        return sqrt(sumsquare(res)/len(res))
    
    def composition(self,formula):
        """
        works through a formula string to determine the elemental composition
        """
        sbrack = ['(','{','['] # start brackets
        ebrack = [')','}',']'] # closing brackets
        def interpret(block):
            """
            interprets an element block, breaking it into element and number of that element
            """
            
            if block[0].isdigit() is True: # if isotope number is encountered
                return {block:1}
            else:
                ele = block[0]
                i = 0
                num = ''
                while i < len(block)-1:
                    i+=1
                    if block[i].isdigit() is True: # add digits
                        num += block[i]
                    else:
                        ele += block[i]
                if num == '':
                    num = 1
                else:
                    num = int(num)
                return {ele:num}
        
        def bracket(form):
            """
            finds the string block contained within a bracket and determines the formula within that bracket
            """
            bracktype = sbrack.index(form[0]) # sets bracket type (so close bracket can be identified)
            bnum = '' # number of things indicated in the bracket
            block = '' # element block
            nest = 1 # counter for nesting brackets
            for loc in range(len(form)): # look for close bracket
                if loc == 0:
                    continue
                elif form[loc] == sbrack[bracktype]: # if a nested bracket is encountered
                    nest += 1
                    block += form[loc]
                elif form[loc] == ebrack[bracktype]: # if close bracket is encountered
                    nest -= 1
                    if nest == 0:
                        i = loc+1 # index of close bracket
                        break
                    else:
                        block += form[loc]
                else:
                    block += form[loc]
            
            try: # look for digits outside of the bracket
                while form[i].isdigit() is True: 
                    bnum += form[i]
                    i+=1
            except IndexError: # if i extends past the length of the formula
                pass
            except UnboundLocalError: # if a close bracket was not found, i will not be defined
                raise ValueError('A close bracket was not encountered for the "%s" bracket in the formula segment "%s". Please check your input molecular formula.' %(form[0],form))
            
            lblock = len(block)+len(bnum)+2 # length of the internal block + the length of the number + 2 for the brackets
            if bnum == '': # if no number is specified
                bnum = 1
            else:
                bnum = int(bnum)
            outdict = {}
            while len(block) > 0: # chew through bracket
                ftemp,tempdict = chewformula(block)
                for key in tempdict:
                    try:
                        outdict[key] += tempdict[key]*bnum
                    except KeyError:
                        outdict[key] = tempdict[key]*bnum
                block = ftemp
            return form[lblock:],outdict # returns remaining formula and composition of the block
        
        def chewformula(formula):
            """
            Iterates through provided formula, extracting blocks, interpreting the blocks, and returning the formula minus the blocks
            """
            if formula[0].isupper() is True: # element is recognized by an uppercase letter
                block = formula[0] # element block
                for loc in range(len(formula)):
                    if loc == 0:
                        continue
                    if formula[loc].isupper() is True: # if an uppercase character is encountered
                        break
                    elif formula[loc] in sbrack: # if a bracket is encountered
                        break
                    else:
                        block += formula[loc]
                return formula[len(block):],interpret(block) # return remaining formula and the interpreted block
            elif formula[0] in sbrack:
                return bracket(formula)
            elif formula[0].isdigit() is True: # either isotope or charge
                for ind,val in enumerate(formula):
                    if formula[ind].isalpha() is True: # if isotope encountered, return that isotope with n=1
                        return '',{formula:1}
                self.charge,self.sign = self.interpretcharge(formula) # otherwise, interpret as charge and return empty dict
                return '',{}
                
        
        def abbreviations(dic):
            """looks for predefined common abbreviations"""
            from _formabbrvs import abbrvs # import dictionary of common abbreviations
            comptemp = {}
            for key in dic:
                if key in abbrvs: # if a common abbreviation is found in formula
                    for subkey in abbrvs[key]:
                        try:
                            comptemp[subkey] += abbrvs[key][subkey]*dic[key]
                        except KeyError:
                            comptemp[subkey] = abbrvs[key][subkey]*dic[key]
                else:
                    try:
                        comptemp[key] += dic[key]
                    except KeyError:
                        comptemp[key] = dic[key]
            return comptemp
        
        def aas(dic):
            """
            looks for one-letter amino acid keys in the formula and converts them to three letter keys
            this function is currently unused and broken (and also likely not useful)
            there are some one-letter keys which are the same as element symbols, which will lead to unexpected behaviour
            The author recommend generating a separate class to handle amino acids (or use pyteomics)
            """
            from _formabbrvs import aminoacids # import dictionary of one-letter amino acid abbreviations
            comptemp = {}
            for key in dic:
                if key in aminoacids:
                    comptemp[aminoacids[key]] = dic[key]
                else:
                    comptemp[key] = dic[key]
            return comptemp
        
        comp = {}
        while len(formula) > 0: # chew through formula
            ftemp,nomdict = chewformula(formula) # find the next block   
            for ele in nomdict:
                try:
                    comp[ele] += nomdict[ele]
                except KeyError:
                    comp[ele] = nomdict[ele]
            formula = ftemp
        comp = abbreviations(comp) # look for common abbreviations    
        return comp
    
    def default(self):
        """saves the original values when the class was called"""
        self.original = dict(self.__dict__)
    
    def exactmass(self,comp,charge=1):
        """
        a quick estimation of the exact mass given a molecular formula
        This may not be exact for high mass species
        """
        em = 0.
        for key in comp:
            try:
                em += self.md[key][0][0]*comp[key]
            except KeyError:
                ele,iso = self.isotope(key)
                em += self.md[ele][iso][0]*comp[key]
        ## accounts for the mass of an electron (barely affects the mass)
        #if self.sign == '+': 
        #    em -= (9.10938356*10**-28)*charge
        #if self.sign == '_':
        #    em += (9.10938356*10**-28)*charge
        return em/charge
    
    def gaussianisotopepattern(self,verbose=False):
        """
        simulates the isotope pattern obtained in a mass spectrometer by applying a gaussian distribution to a bar isotope pattern with a given resolution
        """
        import numpy as np
        import matplotlib.mlab as mlab
        from _Spectrum import Spectrum
        
        def normaldist(center,fwhm,height,step=0.001):
            """
            generates a normal distribution about the center with the full width at half max specified
            y values will be normalized to the height specified
            
            requires:
            import numpy as np
            import matplotlib.mlab as mlab
            """
            x = np.arange(center-fwhm*2,center+fwhm*2,step)
            y = mlab.normpdf(x,center,self.sigma) # generate normal distribution
            y /= max(y) #normalize
            y *= height #scale to height
            return [x.tolist(),y.tolist()]
        
        if verbose is True:
            import sys
        self.nsp = Spectrum(3,startmz=min(self.barip[0])-self.fwhm*2,endmz=max(self.barip[0])+self.fwhm*2) # generate Spectrum object to encompass the entire region
        for ind,val in enumerate(self.barip[0]): # generate normal distributions for each peak
            if verbose is True:
                sys.stdout.write('\rSumming m/z %.3f %d/%d' %(val,ind+1,len(self.barip[0])))
            nd = normaldist(val,self.fwhm,self.barip[1][ind]) # generate normal distribution for that peak
            self.nsp.addspectrum(nd[0],nd[1]) # add the generated spectrum to the total spectrum
        self.nsp.normalize() # normalize
        self.gausip = self.nsp.trim() # trim None values and output
        return self.gausip 
    
    def interpretcharge(self,string):
        """interprets a charge string and sets values"""
        modes = ['+','-']
        value = ''
        sign = '+' # default value for sign
        if type(string) is int:
            return string,sign
        for ind,val in enumerate(string):
            if val in modes: # if val sets mode
                sign = val
            else: # number
                value += val
        return int(value),sign
    
    def isotope(self,string):
        """tries to interpret an undefined key as an isotope and raises an error if it fails"""
        iso = string[0]
        ele = ''
        i = 1
        try:
            while i < len(string):
                if string[i].isdigit() is True:
                    iso += string[i]
                    i += 1
                if string[i].isalpha() is True:
                    ele += string[i]
                    i += 1
            return ele,int(iso)
        except ValueError:
            raise ValueError('The element "%s" could not be found in either the predefined common abbreviations, in the NIST element database, nor interpreted as an isotope, please check your input.' %(string))

    def molecularformula(self):
        """generates the molecular formula as the string"""
        out = ''
        for key,val in sorted(self.comp.items()):
            out += key
            if self.comp[key] > 1:
                out += str(self.comp[key])
        return out
    
    def molecularweight(self):
        """determines the molecular weight from natural abundances"""
        mwout = 0
        pcompout = {}
        for element in self.comp:
            try:
                for isotope in self.md[element]:
                    if isotope == 0:
                        continue
                    mwout += self.md[element][isotope][0]*self.md[element][isotope][1]*self.comp[element] # add every isotope times its natural abundance times the number of that element
                    try:
                        pcompout[element] += self.md[element][isotope][0]*self.md[element][isotope][1]*self.comp[element]
                    except KeyError:
                        pcompout[element] = self.md[element][isotope][0]*self.md[element][isotope][1]*self.comp[element]
            except KeyError: # if isotope
                ele,iso = self.isotope(element)
                mwout += self.md[ele][iso][0]*self.comp[element] # assumes 100% abundance if specified
                pcompout[ele] = self.md[ele][iso][0]*self.comp[element]
        for element in pcompout: # determines the percent composition of each element
            pcompout[element] = pcompout[element]/mwout
        return mwout,pcompout
    
    def printdetails(self):
        """prints the details of the generated molecule"""
        if self.__dict__.has_key('sys') is False:
            self.sys = __import__('sys')
        self.sys.stdout.write('%s\n' %self)
        self.sys.stdout.write('exact mass: %.5f\n' %round(self.em,5))
        self.sys.stdout.write('molecular weight: %.6f\n' %round(self.mw,6))
        self.sys.stdout.write('formula: %s\n' %self.sf)
        self.printpercentcomposition()
    
    def printpercentcomposition(self):
        """prints the percent composition in a readable format"""
        if self.__dict__.has_key('sys') is False:
            self.sys = __import__('sys')
        self.sys.stdout.write('\nelemental percent composition:\n')
        for key,val in sorted(self.pcomp.items()):
            self.sys.stdout.write('%3s: %7.3f %%\n' %(key,self.pcomp[key]*100))
    
    def plotbar(self):
        """quickly plots a bar plot of the isotope bar pattern"""
        import pylab as pl
        fwhm = self.em/self.res
        pl.bar(self.barip[0], self.barip[1], width=fwhm, align='center')
        pl.xlabel('m/z', style='italic')
        pl.ylabel('normalized intensity')
        pl.ticklabel_format(useOffset=False)
        pl.show()
    
    def plotgaus(self,exp=None):
        """quickly plots the simulated gaussian isotope pattern"""
        import pylab as pl
        try:
            pl.plot(self.gausip[0],self.gausip[1],linewidth=1)
        except AttributeError:
            self.gausip = self.gaussianisotopepattern(verbose=False)
            pl.plot(self.gausip[0],self.gausip[1],linewidth=1)
        if exp is not None: # plots experimental if supplied
            y = []
            maxy = max(exp[1])
            for val in exp[1]: # normalize
                y.append(val/maxy*100)
            comp = self.compare(exp)
            pl.plot(exp[0],exp[1])
            pl.text(max(exp[0]),95,'SER: '+`comp`)
            #pl.fill_between(x,self.gausip[1],exp[1],where= exp[1] =< self.gausip[1],interpolate=True, facecolor='red')
        pl.fill(self.gausip[0],self.gausip[1],facecolor='blue',alpha=0.25)
        pl.xlabel('m/z', style='italic')
        pl.ylabel('normalized intensity')
        pl.ticklabel_format(useOffset=False)
        pl.show()
    
    def plotraw(self):
        """quickly plots the raw isotope pattern (with mass defects preserved)"""
        import pylab as pl
        pl.bar(self.rawip[0],self.rawip[1],width=0.0001)
        pl.xlabel('m/z', style='italic')
        pl.ylabel('normalized intensity')
        pl.ticklabel_format(useOffset=False)
        pl.show()
    
    def rawisotopepattern(self,comp,thresh=0.01,dec=5,verbose=False):
        """
        generates an isotope pattern given a molecular formula
        operates using values obtained from the provided mass dictionary in __init__ (default NIST database)
        thresh defines the intensity threshold above which peaks will be tracked (where the max peak height is 100)
        
        returns an uncharged isotope pattern (z will be 1) with all mass defects preserved (infinite resolution)

        supported mass dictionary format is:
        dict = {'element':{0:(monoisotopic mass,1.0),
        isotope#:(exact mass,natural abundance (normalized to 1)),
        ...other isotopes...},
        'next element':...
        ...}
        """
        from _Spectrum import Spectrum
        if verbose is True:
            import sys
            sys.stdout.write('Generating raw isotope pattern.\n')
        out = [[0.],[100.]]
        for key in comp: # for each element
            if self.md.has_key(key) is True: # if natural abundance
                bnds = []
                for mass in self.md[key]:
                    bnds.append(self.md[key][mass][0]) # pull all the masses (used for Spectrum object generation)
                for n in range(comp[key]): # for n number of atoms of each element
                    if verbose is True:
                        sys.stdout.write('\rProcessing element %s %d/%d' %(key,n+1,comp[key]))
                    spec = Spectrum(dec,startmz=min(out[0])+min(bnds),endmz=max(out[0])+max(bnds)) # generate spectrum object
                    for mass in self.md[key]:# for each mass of that element 
                        if mass != 0:
                            if self.md[key][mass][1] != 0: # if intensity is nonzero
                                for ind,val in enumerate(out[0]): # for every current mass in the building isotope pattern
                                    spec.addvalue(out[0][ind]+self.md[key][mass][0],out[1][ind]*self.md[key][mass][1]) # add intensity at the appropriate mass
                    spec.normalize(top=100.) # normalize spectrum
                    spec.threshold(thresh) # drop values below threshold
                    out = spec.trim()
            else: # if specific isotope
                temp = []
                ele,iso = self.isotope(key)
                for ind,val in enumerate(out[0]):
                    temp.append(out[0][ind]+self.md[ele][iso][0])
                out = [list(temp),out[1]]
            if verbose is True:
                sys.stdout.write('\n')
        return out
    
    def reset(self):
        """resets values to when the instance was created"""
        self.__dict__ = self.original
    
    def sigmafwhm(self,res,em):
        """determines the full width at half max and sigma for a normal distribution"""
        import math
        fwhm = em/res
        sigma = fwhm/(2*math.sqrt(2*math.log(2))) # based on the equation FWHM = 2*sqrt(2ln2)*sigma
        return fwhm,sigma
    
    
if __name__ == '__main__': # for testing and troubleshooting
    string = 'L2PdAr+I'
    charge = 1
    res = 5000
    mol = Molecule(string,charge=charge,res=res)
    mol.printdetails()