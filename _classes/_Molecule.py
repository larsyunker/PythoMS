"""
IGNORE:
Molecule class (previously "isotope pattern generator" and "MolecularFormula")

The output of this builder has been validatemwd against values calculated by ChemCalc (www.chemcalc.org)
Negligable differences are attributed to different low value discarding techniques
(ChemCalc keeps the top 5000 peaks, this script drops values less than a threshold 5 orders of magnitude below the maximum value)

CHANGELOG:
- added exact mass comparison
- separated fwhm calculation from sigma
- fwhm calculation now uses monoisotopic mass
- barisotope pattern now groups using the full width at half max
- gaussian isotope pattern generation now works off of rawip by default
- updated to use Progress class
- updated gaussian isotope pattern generator to automatically determine the appropriate decimal places
---2.9 INCOMPATIBLE WITH SPECTRUM v2.4 or older
- moved charge application to raw isotope pattern function
- fixed bug in validation function for charged molecules
- added support for and enabled auto-saving of molecule instances (loading and saving to .mol files)
IGNORE
"""

from PythoMS._classes._ScriptTime import ScriptTime

st = ScriptTime(profile=True)


class Molecule(object):
    def __init__(self, string, **kwargs):
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
        """
        Determines several mass-related properties of a molecular formula
        
        **Parameters**
        
        string: *string*
            The molecular formula to interpret
        
        
        **Examples**
        
        A molecule object can be defined by calling the class with a provided molecular formula.
        
        ::
        
            >>> mol = Molecule('C61H51IP3Pd')
            >>> mol.mimass # monoisotopic mass
            1109.12832052557
            >>> mol.em # estimated exact mass
            1109.1303776683571
            >>> mol.mw # molecular weight
            1110.300954822616
            >>> mol.print_percent_composition() # prints the percent composition

            elemental percent composition:
            C:  65.987 %
            H:   4.630 %
            I:  11.430 %
            P:   8.369 %
            Pd:   9.584 %
        
        
        Common abbreviations can be specified in _formabbrvs.py located in the same
        directory as this class. The abbreviations should be in dictionary format, 
        with the molecular composition of that component. The abbreviations must begin
        with a capital letter, and must not be an elemental formula. e.g.
        
        ::
        
            abbrvs = {
            'Ar+':{'C':25,'H':21,'P':1}, # phosphonium tagged aryl group
            'L':{'C':18,'H':15,'P':1}, # triphenylphosphine
            }
        
        The molecule class could then be called with those abbreviations. 
        
        ::
        
            >>> mol = Molecule('L2PdAr+I')
            >>> mol.em
            1109.1303776683571
        
        Both brackets and bracket nesting are supported. Supported bracket types are
        (,),[,],{,}. e.g.
        
        ::
        
            >>> mol = Molecule('Pd{Ar+}I(P(C6H5)3)')
            >>> mol.em
            1109.1303776683571
        
        Isotopes can be specified with the isotope preceeding the element symbol with 
        both enclosed in a bracket. The mass of that isotope must be defined in the 
        mass dictionary being used by the script (default NIST mass, _nist_mass.py). e.g.
        
        ::
        
            >>> mol1 = Molecule('N(CH2CH3)4')
            >>> mol1.em
            130.159574
            >>> mol2 = Molecule('N([13C]H2CH3)4')
            >>> mol2.em
            134.1729934
        
        The charge of the molecule can be specified either in brackets in the molecular 
        formula, or as a keyword argument (see below). The charge in the molecular formula 
        will override the keyword argument charge. e.g.
        
        ::
            >>> mol = Molecule ('L2PdAr+(2+)')
            >>> mol.em
            491.11295233417843
        
        
        **\*\*kwargs**
    
        charge: 1
            The charge of the molecule. This will affect any properties related to the mass to charge 
            ratio. If the charge is specified in the input molecular formula, this will be 
            overridden. Options: integer. 
        
        consolidate: 3
            When using the consolidate drop method, consolidate peaks within 10^-*consolidate* 
            of each other. See *dropmethod* for more details. Options: integer. 
        
        criticalerror: 3*10**-6 (3 parts per million)
            The critical error value used for warning the user of a potential calculation error. 
            This only affects the ``print_details()`` function output. 
        
        decpl: 7
            The number of decimal places to track while calculating the isotope pattern. 
            Decreasing this will improve efficiency but decrease accuracy. Options: integer. 
        
        dropmethod: None
            The peak drop method to use if desired. Using a peak dropping method will improve 
            calculation times, but decrease the accuracy of the calculated isotope pattern. 
            Options: None, 'threshold', 'npeaks', 'consolidate'. 
            Threshold drops all peaks below a specified threshold value (specified using the 
            *threshold* keyword argument). 
            npeaks keeps the top *n* peaks, specified by the *npeaks* keyword argument. 
            Consolidate combines the intensity of peaks below the threshold value into the 
            nearest peak (within the delta specified by the *consolidate* keyword argument). 
            The new peak *m/z* value is determined by the weighted average of the combined 
            peaks. 
            This will be repeated until the peak is above the threshold or there are no 
            near peaks. 
            The consolidate method is the most accurate. 
        
        emptyspec: True
            Whether to use an empty spectrum obect. Disable this for very large molecules to 
            improve calculation time. Options: bool. 
            
        groupmethod: 'weighted'
            The grouping method to use when calculating the bar isotope pattern from the raw 
            isotope pattern. 
            Options: 'weighted' or 'centroid'. 
            Weighted calculates the peak locations using the weighted average of the *m/z* 
            and intensity values. 
            Centroid finds the center *m/z* value of a group of peaks. 
        
        ipmethod: 'multiplicative'
            The method to use for determining the isotope pattern. 
            'multiplicative' multiplies the existing list of intensities by each element
            'combinatorial' uses combinatorics and iterators to calculate each possible combination
            'hybrid' uses combinatorics to calcuate the pattern from each element, then multiplies those together
        
        keepall: False
            Whether to keep all peaks calculated in the isotope pattern. 
            When false, this will drop all intensities below 0.0001 after calculating the isotope pattern. 
        
        npeaks: 5000
            The number of peaks to keep if *dropmethod* is 'npeaks'. See *dropmethod* for more 
            details. Options: integer. 
        
        res: 5000
            The resolution of the instrument to simulate when generating the gaussian isotope pattern. 
            This also affects the bounds() method. Options: integer or float. 
        
        threshold: 0.01
            The threshold value determining whether or not to drop a peak. Only has an effect 
            if *dropmethod* is not ``None``. See *dropmethod* for more details. Options: float. 
        
        verbose: False
            Verbose output. Mostly useful when calculating for large molecules or while debugging. 
            Options: bool. 
        
        """
        self.kw = {  # default keyword arguments
            'verbose': False,  # toggle verbose
            'decpl': 7,  # number of decimal places to track while generating the raw isotope pattern
            'res': 5000,  # resolution of the instrument being matched
            'charge': 1,  # charge of the molecule (this can also be specified in the formula)
            'emptyspec': True,  # use an empty spectrum object (disable this for massive molecules)
            'ipmethod': 'multiplicative',  # method of calculating the raw isotope pattern
            'dropmethod': None,  # method to drop negligible peaks (either 'threshold' or 'npeaks')
            'threshold': 0.01,
            # threshold to use when dropping negligible peaks (only applies if dropmethod is threshold)
            'npeaks': 5000,
            # number of peaks to keep when dropping negligible peaks (only applies if dropmethod is npeaks)
            'consolidate': 3,
            # when using the consolidate drop method, consolidate values within 10**-value of each other
            'groupmethod': 'weighted',
            # the peak grouping method applied when generating the bar isotope pattern ('weighted' or 'centroid')
            'criticalerror': 3 * 10 ** -6,
            # the critical error at which the script warns regarding a potential mismatch (default 3 ppm)
            'keepall': False,
            # whether to keep all peaks calculated in the isotope pattern (otherwise drops below 0.0001
            'precalculated': None,  # path to precalculated molecule pickle file (optional, this can be disabled by setting this to False)
            'save': True,  # whether to save values
        }
        if set(kwargs.keys()) - set(self.kw.keys()):  # check for invalid keyword arguments
            string = ''
            for i in set(kwargs.keys()) - set(self.kw.keys()):
                string += ' %s' % i
            raise KeyError('Unsupported keyword argument(s): %s' % string)
        self.kw.update(kwargs)  # update defaults with provided keyword arguments

        if self.kw['verbose'] is True:
            self.sys = __import__('sys')
            self.sys.stdout.write('Generating molecule object from input "%s"\n' % string)

        self.np = __import__('numpy')
        from PythoMS._classes._nist_mass import nist_mass as mass  # masses from the NIST database
        # from PythoMS._classes._crc_mass import crc_mass as mass # masses from the CRC Handbook of Chemistry and Physics
        self.md = mass  # mass dictionary that the script will use
        self.formula = string  # input formula
        self.kw['charge'], self.kw['sign'] = self.interpret_charge(self.kw['charge'])  # charge
        self.comp = self.composition(self.formula)  # determine composition from formula
        self.check_in_mass_dict(self.comp)  # checks that all the composition keys are valid
        self.calculate()
        self.default()
        if self.kw['verbose'] is True:
            self.print_details()

    def __str__(self):
        return "Molecule {}".format(self.formula)

    def __repr__(self):
        return "{}('{}')".format(self.__class__.__name__, self.formula)

    def __add__(self, x):
        """
        Several supported addition methods:
        If a valid molecular formula string is provided, that string will be added.
        If another Molecule class instance is provided, the provided instance will be
        added to the current instance.
        """
        if type(x) is str:
            addition = self.composition(x)
        else:
            if isinstance(x, self.__class__) is True:
                addition = x.comp
            else:
                raise ValueError('Addition of {} to Molecule object {} is invalid'.format(x, self.formula))
        for key in addition:
            try:
                self.comp[key] += addition[key]
            except KeyError:
                self.comp[key] = addition[key]
        self.calculate()  # recalculates mass and isotope patterns
        self.formula = self.sf
        return self.sf

    def __sub__(self, x):
        """
        See __add__ for details.
        Subtract has a catch for a negative number of a given element
        (the minimum that can be reached is zero).
        """
        if type(x) is str:
            addition = self.composition(x)
        else:
            if isinstance(x, self.__class__) is True:
                addition = x.comp
            else:
                raise ValueError('Subtraction of {} from Molecule object {} is invalid'.format(x, self.formula))
        for key in addition:
            try:
                if self.comp[key] - addition[key] > 0:
                    self.comp[key] -= addition[key]
                elif self.comp[key] - addition[key] is 0:
                    del self.comp[key]
                else:
                    raise ValueError('Subtracting %d number of element %s from %s would yield a negative amount.' % (
                        addition[key], key, self.sf))
            except KeyError:
                return 'The molecular formula {} does not have an element in {}'.format(self.formula, x)
        self.calculate()  # recalculates mass and isotope patterns
        self.formula = self.sf
        return self.sf

    def __mul__(self, x):
        """allows integer multiplication of the molecular formula"""
        if type(x) != int:
            raise ValueError('Non-integer multiplication of a Molecule object is unsupported')
        for key in self.comp:
            self.comp[key] = self.comp[key] * x
        self.calculate()  # recalculates mass and isotope patterns
        self.formula = self.sf
        return self.sf

    def __div__(self, x):
        """allows integer division of the molecular formula"""
        if type(x) != int:
            raise ValueError('Non-integer division of a Molecule object is unsupported')
        for key in self.comp:
            tempval = float(self.comp[key]) / float(x)
            if tempval.is_integer() is False:
                raise ValueError('Division of %s (%s) by %d yielded a non-integer number of element %s' % (
                    self.formula, self.sf, x, key))
            else:
                self.comp[key] = int(tempval)
        self.calculate()  # recalculates mass and isotope patterns
        self.formula = self.sf
        return self.sf

    def bar_isotope_pattern(self, rawip, delta=0.5):
        """
        Converts a raw isotope pattern into a bar isotope pattern. This groups mass defects
        that are within a given difference from each other into a single *m/z* value and
        intensity.

        **Parameters**

        rawip: *list*
            The raw isotope pattern of the class object.

        delta: *float*
            The *m/z* difference to check around a peak when grouping it into a single *m/z* value.
            The script will look delta/2 from the peak being checked


        **Returns**

        return item: *list*
            Paired lists in ``[[m/z values],[intensity values]]`` format.

        """

        def groupmasses(ip, dm):
            """
            groups masses in an isotope pattern
            looks for differences in m/z greater than the specified delta
            expects a paired list of [[mz values],[intensity values]]
            """
            num = 0
            out = [[[], []]]
            for ind, val in enumerate(ip[0]):
                out[num][0].append(ip[0][ind])
                out[num][1].append(ip[1][ind])
                try:
                    if ip[0][ind + 1] - ip[0][ind] > dm:
                        num += 1
                        out.append([[], []])
                except IndexError:
                    continue
            return out

        def centroid(ipgroup):
            """
            takes a group of mz and intensity values and finds the centroid
            this method results in substantially more error compared to the weighted_average method (~9 orders of magnitude for C61H51IP3Pd)
            """
            return sum(ipgroup[0]) / len(ipgroup[0]), sum(ipgroup[1]) / len(ipgroup[1])

        if self.kw['verbose'] is True:
            self.sys.stdout.write('Generating bar isotope pattern')
        groupedip = groupmasses(rawip, delta / 2)
        out = [[], []]
        for group in groupedip:
            if self.kw['groupmethod'] == 'weighted':
                x, y = self.weighted_average(group)  # determine weighted mass and summed intensity
            elif self.kw['groupmethod'] == 'centroid':
                x, y = centroid(group)
            out[0].append(x)
            out[1].append(y)
        maxint = max(out[1])
        for ind, val in enumerate(out[1]):
            out[0][ind] = out[0][ind]  #/ abs(charge)
            out[1][ind] = val / maxint * 100.  # normalize to 100
        if self.kw['verbose'] is True:
            self.sys.stdout.write(' DONE\n')
        return out

    def bounds(self, conf=0.95, perpeak=False, threshold=0.01):
        """
        calculates bounds of the isotope pattern based on a confidence interval and the bar isotope pattern

        conf: (float) the confidence interval to use
        perpeak: (bool) toggle for whether the function should return a dictionary of
        boundaries for each peak, or a single pair of bounds that covers the entire isotope pattern
        threshold: (int/float) minimum threshold as a percentage of the maximmum for peaks to be included in bounds
        """
        """
        Calculates the *m/z* bounds of the isotope pattern of the molecule object based 
        on a confidence interval and the *m/z* values of the bar isotope pattern. 
        This can be used to automatically determine the integration bounds required to 
        contain XX% of the counts associated with that molecule in a mass spectrum. 
        
        
        **Parameters**
        
        conf: *float*
            The confidence interval to use for calculating the bounds. 
            e.g. *0.95* corresponds to a 95% confidence interval. 
        
        perpeak: *bool*
            Whether or not to return the bounds required to integrate each 
            peak of the isotope pattern individually. 
            This can be useful in a very noisy mass spectrum to avoid 
            baseline noise within the integration interval. 
        
        threshold: *float*
            The threshold used to determine whether a peak should be 
            included in the bounds. 
        
        **Returns**
        
        bounds: *list* or *dictionary*
            If *perpeak* is False, this will return a two item list 
            corresponding to the start and end *m/z* bounds. 
            If *perpeak* is True, returns a dictionary of bounds with 
            the key format of 
            ``dict[parent m/z value]['bounds'] = [start m/z, end m/z]``
        
        
        **Examples**
        
        To determine the integration bounds of C61H51IP3Pd: 
            
        ::
        
            >>> mol = Molecule('C61H51IP3Pd')
            >>> mol.bounds(0.95)
            [1104.9458115053008, 1116.3249999321531]
            
            >>> mol.bounds(0.99)
            [1104.8877964620444, 1116.3830149754094]
            
            >>> mol.bounds(0.95,True)
            {'1105.1304418': {'bounds': (1104.9458115053008, 1105.3150720946992)},
            '1106.13382235': {'bounds': (1105.9491920547823, 1106.3184526441808)},
            '1107.12903188': {'bounds': (1106.9444015896975, 1107.3136621790959)},
            '1108.13051519': {'bounds': (1107.9458848935217, 1108.3151454829201)},
            '1109.13037767': {'bounds': (1108.9457473736579, 1109.3150079630564)},
            '1110.13288962': {'bounds': (1109.9482593265234, 1110.3175199159218)},
            '1111.13024042': {'bounds': (1110.9456101206658, 1111.3148707100643)},
            '1112.13263766': {'bounds': (1111.9480073654438, 1112.3172679548422)},
            '1113.13193341': {'bounds': (1112.9473031156144, 1113.3165637050129)},
            '1114.13415503': {'bounds': (1113.9495247326277, 1114.3187853220261)},
            '1115.13715205': {'bounds': (1114.9525217596001, 1115.3217823489986)},
            '1116.14036964': {'bounds': (1115.9557393427547, 1116.3249999321531)}}
        
        """
        if self.kw['verbose'] is True:
            self.sys.stdout.write('Calculating bounds from simulated gaussian isotope pattern')
        threshold = threshold * max(self.barip[1])
        from scipy import stats
        tempip = [[], []]
        for ind, inten in enumerate(self.barip[1]):  # checks for intensities above threshold
            if inten >= threshold:
                tempip[0].append(self.barip[0][ind])
                tempip[1].append(self.barip[1][ind])
        if perpeak is True:  # if per-peak bounds are called for
            out = {}
            for mz in tempip[0]:
                out[str(mz)] = {}
                out[str(mz)]['bounds'] = stats.norm.interval(conf, mz, scale=self.sigma)
        else:  # a general range that covers the entire isotope pattern
            out = [stats.norm.interval(conf, tempip[0][0], scale=self.sigma)[0],
                   stats.norm.interval(conf, tempip[0][-1], scale=self.sigma)[1]]
        if self.kw['verbose'] is True:
            if perpeak is False:
                self.sys.stdout.write(': %.3f-%.3f' % (out[0], out[1]))
            self.sys.stdout.write(' DONE\n')
        return out

    def calculate(self):
        """This is a function which calls all the calculation functions in the class"""
        if self.kw['precalculated'] is not False and self.load_from_pickle(self.kw['precalculated']) is True:
            return None  # try to load from pickle
        self.sf = self.molecular_formula()  # generates a string version of the molecular formula
        self.mimass = self.monoisotopic_mass()  # monoisotopic mass
        self.fwhm = self.mimass / self.kw['res']  # calculate the full width at half max
        self.mw, self.pcomp = self.molecular_weight()  # molecular weight and elemental percent composition
        if self.kw['ipmethod'] == 'combinatorics':
            self.rawip = self.isotope_pattern_combinatorics(self.comp)
        elif self.kw['ipmethod'] == 'multiplicative':
            self.rawip = self.isotope_pattern_multiplicative(self.comp)  # generates a raw isotope pattern (charge of 1)
        elif self.kw['ipmethod'] == 'hybrid':
            self.rawip = self.isotope_pattern(self.comp)
        self.barip = self.bar_isotope_pattern(self.rawip, self.fwhm)  # bar isotope pattern based on the generated raw pattern
        self.em = self.estimated_exact_mass()
        self.sigma = self.standard_deviation(self.fwhm)
        self.error = self.validate_pattern(self.mw, self.rawip)
        if self.kw['save'] is True:  # save if specified
            self.save_to_pickle(self.kw['precalculated'])
        # self.gausip = self.gaussian_isotope_pattern(self.barip,self.em,res=self.kw['res']) # simulated normal distribution of the bar isotope pattern
        ## gausip is rarely used and is a bit more time consuming, so it now requires a specific call. If call is desired on generation, uncomment the above line

    def check_in_mass_dict(self, comp):
        """checks for each 'element' in the nist dictionary and returns an error if not found"""
        for key in comp:
            if key not in self.md:
                ele, iso = self.isotope(key)
                if iso not in self.md[ele]:
                    raise ValueError(
                        'The element "%s" does not have a defined isotope "%d" in the NIST element database, please check your input' % (
                            ele, iso))

    def compare(self, exp):
        """
        Compares a provided mass spectrum (experimental) to the simulated gaussian
        isotope pattern. Returns a standard error of the regression as an assessment
        of the goodness of fit.

        **Parameters**

        exp: *list*
            The experimentally acquired mass spectra provided as a paired list of lists
            ``[[m/z values],[intensity values]]``


        **Returns**

        Standard error of the regression: *float*
            A measure of the average distance between the experimental and predicted
            values. Lower is better, although this is a qualitative assessment.

        """

        def sumsquare(lst):
            """calculates the sum of squares"""
            ss = 0
            for val in lst:
                ss += val ** 2
            return ss

        if 'gausip' not in self.__dict__:  # generate gaussian isotope pattern if not already generated
            self.gaussian_isotope_pattern(self.rawip)
        yvals = []
        res = []
        maxy = float(max(exp[1]))
        if maxy == 0.:
            return 'could not calculate'
        for ind, val in enumerate(exp[1]):  # normalize y values
            yvals.append(float(val) / maxy * 100.)
        # avgy = sum(exp[1])/len(exp[1])
        for ind, mz in enumerate(exp[0]):
            if mz > min(self.gausip[0]) and mz < max(self.gausip[0]):  # if within isotope pattern
                nspind = self.spec.index(mz)  # calculate index
                if self.spec.y[nspind] is not None:  # if the predicted intensity is not None
                    res.append(yvals[ind] - self.spec.y[nspind])  # difference between observed and predited (residuals)
                    # tot.append(self.spec.y[nspind]-avgy) # difference between predicted and mean
        # rsqrd = 1-(sumsquare(res)/sumsquare(tot)) # r-squared value (apparently not applicable to non-linear fits)
        return self.np.sqrt(sumsquare(res) / len(res))

    def compare_exact_mass(self, mass, use='est'):
        """
        Compares the provided mass to the exact mass of the calculated molecule.

        **Parameters**

        mass: *float*
            experimental mass to compare

        use: est or mi (optional)
            Whether to compare the estimated exact mass or the monoisotopic
            mass to the provided value. Default: est

        **Returns**

        relative error: *float*
            The relative error of the provided mass to the exact mass
        """
        if use == 'est':
            delta = mass - self.em
            return delta / self.em * 10 ** 6
        elif use == 'mi':
            delta = mass - self.mimass
            return delta / self.mimass * 10 ** 6

    def composition(self, formula):
        """
        works through a formula string to determine the elemental composition
        """
        """
        Interprets a provided string as a molecular formula. 
        Supports nested brackets, charges, and isotopes (see __init__ docstring 
        for more details). 
        
        **Parameters**
        
        provided formula: *string*
            The provided molecular formula. 
        
        
        **Returns**
        
        composition dictionary: *dictionary*
            A dictionary where each key is an element or isotope with its value 
            being the number of each of the elements or isotopes. e.g. the 
            molecule CH4 would have the composition ``comp = {'C':1, 'H':4}
        
        **See Also**
        
        See the __init__ docstring for more details. 
        
        """
        sbrack = ['(', '{', '[']  # start brackets
        ebrack = [')', '}', ']']  # closing brackets

        def interpret(block):
            """interprets an element block, breaking it into element and number of that element"""
            if block[0].isdigit() is True:  # if isotope number is encountered
                return {block: 1}
            else:
                ele = block[0]
                i = 0
                num = ''
                while i < len(block) - 1:
                    i += 1
                    if block[i].isdigit() is True:  # add digits
                        num += block[i]
                    else:
                        ele += block[i]
                if num == '':
                    num = 1
                else:
                    num = int(num)
                return {ele: num}

        def bracket(form):
            """finds the string block contained within a bracket and determines the formula within that bracket"""
            bracktype = sbrack.index(form[0])  # sets bracket type (so close bracket can be identified)
            bnum = ''  # number of things indicated in the bracket
            block = ''  # element block
            nest = 1  # counter for nesting brackets
            for loc in range(len(form)):  # look for close bracket
                if loc == 0:
                    continue
                elif form[loc] == sbrack[bracktype]:  # if a nested bracket is encountered
                    nest += 1
                    block += form[loc]
                elif form[loc] == ebrack[bracktype]:  # if close bracket is encountered
                    nest -= 1
                    if nest == 0:
                        i = loc + 1  # index of close bracket
                        break
                    else:
                        block += form[loc]
                else:
                    block += form[loc]

            try:  # look for digits outside of the bracket
                while form[i].isdigit() is True:
                    bnum += form[i]
                    i += 1
            except IndexError:  # if i extends past the length of the formula
                pass
            except UnboundLocalError:  # if a close bracket was not found, i will not be defined
                raise ValueError(
                    'A close bracket was not encountered for the "%s" bracket in the formula segment "%s". Please check your input molecular formula.' % (
                        form[0], form))

            lblock = len(block) + len(
                bnum) + 2  # length of the internal block + the length of the number + 2 for the brackets
            if bnum == '':  # if no number is specified
                bnum = 1
            else:
                bnum = int(bnum)
            outdict = {}
            while len(block) > 0:  # chew through bracket
                ftemp, tempdict = chewformula(block)
                for key in tempdict:
                    try:
                        outdict[key] += tempdict[key] * bnum
                    except KeyError:
                        outdict[key] = tempdict[key] * bnum
                block = ftemp
            return form[lblock:], outdict  # returns remaining formula and composition of the block

        def chewformula(formula):
            """
            Iterates through provided formula, extracting blocks, interpreting the blocks,
            and returning the formula minus the blocks
            """
            if formula[0].isupper() is True:  # element is recognized by an uppercase letter
                block = formula[0]  # element block
                for loc in range(len(formula)):
                    if loc == 0:
                        continue
                    if formula[loc].isupper() is True:  # if an uppercase character is encountered
                        break
                    elif formula[loc] in sbrack:  # if a bracket is encountered
                        break
                    else:
                        block += formula[loc]
                return formula[len(block):], interpret(block)  # return remaining formula and the interpreted block
            elif formula[0] in sbrack:
                return bracket(formula)
            elif formula[0].isdigit() is True:  # either isotope or charge
                for ind, val in enumerate(formula):
                    if formula[ind].isalpha() is True:  # if isotope encountered, return that isotope with n=1
                        return '', {formula: 1}
                self.kw['charge'], self.kw['sign'] = self.interpret_charge(
                    formula)  # otherwise, interpret as charge and return empty dict
                return '', {}

        def abbreviations(dic):
            """
            looks for predefined common abbreviations
            These can be defined in _formabbrvs.py
            """
            from PythoMS._classes._formabbrvs import abbrvs  # import dictionary of common abbreviations
            comptemp = {}
            for key in dic:
                if key in abbrvs:  # if a common abbreviation is found in formula
                    for subkey in abbrvs[key]:
                        try:
                            comptemp[subkey] += abbrvs[key][subkey] * dic[key]
                        except KeyError:
                            comptemp[subkey] = abbrvs[key][subkey] * dic[key]
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
            from PythoMS._classes._formabbrvs import \
                aminoacids  # import dictionary of one-letter amino acid abbreviations
            comptemp = {}
            for key in dic:
                if key in aminoacids:
                    comptemp[aminoacids[key]] = dic[key]
                else:
                    comptemp[key] = dic[key]
            return comptemp

        if self.kw['verbose'] is True:
            self.sys.stdout.write('Determining composition from supplied molecular formula')
        comp = {}
        while len(formula) > 0:  # chew through formula
            ftemp, nomdict = chewformula(formula)  # find the next block
            for ele in nomdict:
                try:
                    comp[ele] += nomdict[ele]
                except KeyError:
                    comp[ele] = nomdict[ele]
            formula = ftemp
        comp = abbreviations(comp)  # look for common abbreviations
        if self.kw['verbose'] is True:
            self.sys.stdout.write(' DONE\n')
        return comp

    def default(self):
        """saves the original values when the class was called"""
        self.original = dict(self.__dict__)

    def estimated_exact_mass(self):
        """determines the precise exact mass from the bar isotope pattern"""
        ind = self.barip[1].index(100.)
        return self.barip[0][ind]

    def gaussian_isotope_pattern(self, barip=None):
        """
        Simulates the isotope pattern that would be observed in a mass
        spectrometer with the resolution specified in the keyword arguments
        of the class.

        **Parameters**

        barip: *list*
            The isotope pattern to be simulated. This can be either the bar isotope
            pattern or the raw isotope pattern (although this will be substantially
            slower for large molecules). If unspecified, uses self.rawip


        **Returns**

        gaussian isotope pattern: *list*
            The predicted gaussian isotope pattern in the form of a paired list
            ``[[m/z values],[intensity values]]``
        """
        import numpy as np
        import matplotlib.mlab as mlab
        from PythoMS._classes._Spectrum import Spectrum

        def autodec(fwhm):
            "automatically calculates the appropriate decimal place to track"""
            shift = fwhm
            n = 0
            while shift < 1.:
                n += 1
                shift = fwhm * 10 ** n
            return n + 1  # track 1 higher

        dec = autodec(self.fwhm)

        def normaldist(center, fwhm, height, step=10 ** -dec):
            """
            generates a normal distribution about the center with the full width at half max specified
            y values will be normalized to the height specified

            requires:
            import numpy as np
            import matplotlib.mlab as mlab
            """
            x = np.arange(center - fwhm * 2, center + fwhm * 2, step)
            y = mlab.normpdf(x, center, self.sigma)  # generate normal distribution
            y /= max(y)  # normalize
            y = y * height
            return [x.tolist(), y.tolist()]

        if self.kw['verbose'] is True:
            self.sys.stdout.write('Generating simulated isotope pattern using a resolution of %.0f' % self.kw['res'])

        if barip is None:
            barip = self.rawip

        self.spec = Spectrum(
            dec,
            start=min(barip[0]) - self.fwhm * 2,
            end=max(barip[0]) + self.fwhm * 2,
            empty=False,  # self.kw['emptyspec'], # whether or not to use emptyspec
            filler=0.,  # fill with zeros, not None
        )  # generate Spectrum object to encompass the entire region
        for ind, val in enumerate(barip[0]):  # generate normal distributions for each peak
            # if verbose is True:
            #    sys.stdout.write('\rSumming m/z %.3f %d/%d' %(val,ind+1,len(self.barip[0])))
            nd = normaldist(val, self.fwhm, barip[1][ind])  # generate normal distribution for that peak
            self.spec.addspectrum(nd[0], nd[1])  # add the generated spectrum to the total spectrum
        self.spec.normalize()  # normalize
        self.gausip = self.spec.trim()  # trim None values and output
        if self.kw['verbose'] is True:
            self.sys.stdout.write(' DONE\n')
        return self.gausip

    def interpret_charge(self, string):
        """interprets a charge string and sets values"""
        modes = ['+', '-']
        value = ''
        sign = '+'  # default value for sign
        if type(string) is int:
            return string, sign
        for ind, val in enumerate(string):
            if val in modes:  # if val sets mode
                sign = val
            else:  # number
                value += val
        return int(value), sign

    def isotope(self, string):
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
            return ele, int(iso)
        except ValueError:
            raise ValueError(
                'The element "%s" could not be found in either the predefined common abbreviations, '
                'in the NIST element database, nor interpreted as an isotope, please check your input.' % (string))

    def isotope_pattern(self, comp):
        """
        A hybrid isotope pattern calculator which calculates the isotope pattern from
        each element, then multiplies the list
        """
        from PythoMS._classes._Spectrum import Spectrum
        eleips = {}  # dictionary for storing the isotope patterns of each element
        for element in comp:
            eleips[element] = self.isotope_pattern_combinatorics(
                {element: comp[element]})  # calculate the isotope pattern for each element

        sortlist = []
        for element in eleips:
            sortlist.append((len(eleips[element][0]), element))
        sortlist = sorted(sortlist)  # sorted list of elements based on the length of their isotope patterns
        sortlist.reverse()

        if self.kw['verbose'] is True:
            from PythoMS._classes._Progress import Progress
            prog = Progress(percent=False, fraction=False)

        spec = None
        for lenlist, element in sortlist:
            if self.kw['verbose'] is True:
                prog.updatestring('Adding element %s to isotope pattern' % element)
                prog.write()
            if spec is None:
                spec = Spectrum(
                    self.kw['decpl'],  # decimal places
                    start=None,  # minimum mass
                    end=None,  # maximum mass
                    specin=eleips[element],  # supply masses and abundances as initialization spectrum
                    empty=self.kw['emptyspec'],  # whether or not to use emptyspec
                    filler=0.,  # fill with zeros, not None
                )
                if self.kw['verbose'] is True:
                    prog.fin()

                continue
            spec.addelement(eleips[element][0], eleips[element][1])
            spec.normalize(100.)  # normalize spectrum object
            if self.kw['dropmethod'] == 'threshold':  # drop values below threshold
                spec.threshold(self.kw['threshold'])
            elif self.kw['dropmethod'] == 'npeaks':  # keep top n number of peaks
                spec.keeptopn(self.kw['npeaks'])
            elif self.kw['dropmethod'] == 'consolidate':  # consolidate values being dropped
                spec.consolidate(self.kw['threshold'], 3 * 10 ** -self.kw['consolidate'])
            if self.kw['verbose'] is True:
                self.sys.stdout.write(' DONE\n')
        if self.kw['keepall'] is False:
            spec.threshold(0.0001)  # drop very low intensity
        out = spec.trim()
        return out

    @st.profilefn
    def isotope_pattern_combinatorics(self, comp):
        """
        Calculates the raw isotope pattern of a given molecular formula with mass defects preserved.
        Uses a combinatorial method to generate isotope formulae

        """
        from itertools import combinations_with_replacement as cwr
        import sympy as sym
        from PythoMS._classes._Spectrum import Spectrum

        class Reiterable(object):
            def __init__(self, isos, number):
                """a reiterable version of combinations with replacments iterator"""
                self.isos = isos  # isotopoes group
                self.number = number  # number of atoms of the element

            def __iter__(self):
                return cwr(self.isos, self.number)

        @st.profilefn
        def num_permu(lst, isos):
            """
            calculates the number of unique permutations of the given set of isotopes
            the calculation is generated as a sympy function before evaluation
            numpy factorial is limited in the size of factorials that are calculable
            """
            counts = [lst.count(x) for x in isos]  # counts the number of each isotope in the set
            num = sym.factorial(len(lst))  # numerator is the factorial of the length of the list
            denom = 1  # denominator is the product of the factorials of the counts of each isotope in the list
            for count in counts:
                denom *= sym.factorial(count)
            return (num / denom).evalf()  # divide, evaluate, and return

        @st.profilefn
        def product(*iterables):
            """
            cartesian product of iterables
            from http://stackoverflow.com/questions/12093364/cartesian-product-of-large-iterators-itertools
            """
            if len(iterables) == 0:
                yield ()
            else:
                it = iterables[0]
                for item in it() if callable(it) else iter(it):
                    for items in product(*iterables[1:]):
                        yield (item,) + items

        @st.profilefn
        def numberofcwr(n, k):
            """
            calculates the number of combinations with repitition
            n: number of things to choose from
            k: choose k of them
            """
            fn = sym.factorial(n + k - 1)
            fn /= sym.factorial(k)
            fn /= sym.factorial(n - 1)
            return fn.evalf()

        def listproduct(iterable):
            """returns the product of a list"""
            prod = 1
            for n in iterable:
                prod *= n
            return prod

        speciso = False  # set state for specific isotope
        isos = {}  # isotopes dictionary
        isosets = {}  # set of isotopes for each element
        iterators = []  # list of iterators
        nk = []
        for element in comp:  # for each element
            if element in self.md:
                isosets[element] = []  # set of isotopes
                for isotope in self.md[element]:  # for each isotope of that element in the mass dictionary
                    if isotope != 0 and self.md[element][isotope][1] != 0:  # of the intensity is nonzero
                        isosets[element].append(isotope)  # track set of isotopes
                        isos[isotope] = element  # create isotope,element association for reference
                iterators.append(
                    Reiterable(  # create iterator instance
                        isosets[element],
                        comp[element]
                    )
                )
                if self.kw['verbose'] is True:
                    nk.append([  # track n and k for list length output
                        len(isosets[element]),
                        comp[element]
                    ])
            else:  # if it's an isotope
                speciso = True

        spec = Spectrum(  # initiate spectrum object
            self.kw['decpl'],  # decimal places
            start=None,  # no minimum mass
            end=None,  # no maximum mass
            empty=True,  # whether or not to use emptyspec
            filler=0.,  # fill with zeros, not None
        )
        if self.kw['verbose'] is True:
            counter = 0  # create a counter
            iterations = listproduct([numberofcwr(n, k) for n, k in nk])  # number of iterations
            from PythoMS._classes._Progress import Progress
            prog = Progress(string='Adding isotope combination to queue', last=iterations)  # create a progress instance

        for comb in product(*iterators):
            if self.kw['verbose'] is True:
                counter += 1
                # remaining = st.progress(counter,iterations,'combinations')
                # self.sys.stdout.write('\rProcessing isotope combination #%d/%d %.1f%% (remaining: %s)' %(counter,iterations,float(counter)/float(iterations)*100.,remaining))
                if len(comp) == 1:
                    string = 'for %s%d ' % (comp.keys()[0], comp[comp.keys()[0]])
                else:
                    string = ''
                try:
                    self.sys.stdout.write('\rProcessing isotope combination #%d/%d %.1f%% %s' % (
                        counter, iterations, float(counter) / float(iterations) * 100., string))
                except ValueError:
                    pass
            num = 1  # number of combinations counter
            x = 0.  # mass value
            y = 1.  # intensity value
            for tup in comb:  # for each element combination
                element = isos[tup[0]]  # associate isotope to element
                # counts = [tup.count(x) for x in isosets[element]] # count the number of occurances of each isotope
                # num *= num_permu(tup,counts) # determine the number of permutations of the set
                # for ind,isotope in enumerate(isosets[element]):
                #    x += self.md[element][isotope][0] * counts[ind]
                #    y *= self.md[element][isotope][1] ** counts[ind]
                num *= num_permu(tup, isosets[element])  # multiply the number by the possible permutations
                for isotope in tup:  # for each isotope
                    x += self.md[element][isotope][0]  # shift x
                    y *= self.md[element][isotope][1]  # multiply intensity
            spec.addvalue(x,
                          y * num)  # add the x and y combination factored by the number of times that combination will occur

        if speciso is True:  # if an isotope was specified
            for element in comp:
                if element not in self.md:  # if an isotope
                    ele, iso = self.isotope(element)  # determine element and isotope
                    spec.shiftx(self.md[ele][iso][0] * comp[element])  # shift the x values by the isotopic mass
        if self.kw['charge'] > 1:  # if the magnitude of the charge is greater than 1 divide m/z values
            spec.applycharge(self.kw['charge'])
        spec.normalize()  # normalize the spectrum object
        if self.kw['keepall'] is not False:
            spec.threshold(0.0001)  # drop very low intensity
        out = spec.trim()  # trim to x,y lists
        if self.kw['verbose'] is True:
            self.sys.stdout.write('DONE\n')
        return out

    def isotope_pattern_cuda(self, comp):
        """
        Uses CUDA acceleration to calculate the isotope pattern
        :param comp:
        :return:

        notes to self:
        base this on the combinatorial approach
        1) store mass dictionary in CUDA memory as lookup?
        2) create generators as in the combinatorial approach
        3) create list chunks from the generators and send them to PyCUDA, which will process each tuple
            - determine element
            - calculate the number of permutations (need a good way to do this)
            - add element and value
        3b) permutations could be calcuated by extending the factorial, and evaluating at the end of the tuple set
        4) store of isotopes in CUDA memory
        5) create sorter in CUDA?

        """

    @st.profilefn
    def isotope_pattern_multiplicative(self, comp):
        """
        Calculates the raw isotope pattern of a given molecular formula with mass defects preserved.


        **Parameters**

        comp: *dictionary*
            The molecular composition dictionary. See ``Molecule.composition()`` for more details.

        dec: *integer*
            The number of decimal places to track. This is normally controlled by the keyword
            arguments of the class, but can be specified if called separately.


        **Returns**

        raw isotope pattern: *list*
            Returns the isotope pattern with mass defects preserved (referred to as the 'raw'
            isotope pattern in this script).


        **Additional details**

        The pattern is calculated using masses and abundances retrieved from a mass dictionary
        (e.g. _nist_mass.py or _crc_mass.py)
        The dictionary to be used can be changed by modifying the import line in __init__, and
        the supported mass dictionary format is:

        ::

            dict = {'element':{0:(monoisotopic mass,1.0),
            isotope#:(exact mass,natural abundance (normalized to 1)),
            ...other isotopes...},
            'next element':...
            ...}

        """
        from PythoMS._classes._Spectrum import Spectrum
        spec = None  # initial state of spec
        if self.kw['verbose'] is True:
            from PythoMS._classes._Progress import Progress
            self.sys.stdout.write('Generating raw isotope pattern.\n')
        for key in comp:  # for each element
            if key in self.md:  # if not a single isotope
                if self.kw['verbose'] is True:
                    prog = Progress(string='Processing element %s' % (key), first=0, last=comp[key])
                masses = []  # list for masses of each isotope
                abunds = []  # list for abundances
                for mass in self.md[key]:
                    if mass != 0:
                        if self.md[key][mass][1] > 0:  # if abundance is nonzero
                            masses.append(self.md[key][mass][0])
                            abunds.append(self.md[key][mass][1])
                for n in range(comp[key]):  # for n number of each element
                    if spec is None:  # if spectrum object has not been defined
                        spec = Spectrum(
                            self.kw['decpl'],  # decimal places
                            start=min(masses) - 10 ** -self.kw['decpl'],  # minimum mass
                            end=max(masses) + 10 ** -self.kw['decpl'],  # maximum mass
                            specin=[masses, abunds],  # supply masses and abundances as initialization spectrum
                            empty=self.kw['emptyspec'],  # whether or not to use emptyspec
                            filler=0.,  # fill with zeros, not None
                        )
                        continue
                    if self.kw['verbose'] is True:
                        prog.write(n + 1)
                    spec.addelement(masses, abunds)  # add the element to the spectrum object
                    spec.normalize(100.)  # normalize spectrum
                    if self.kw['dropmethod'] == 'threshold':  # drop values below threshold
                        spec.threshold(self.kw['threshold'])
                    elif self.kw['dropmethod'] == 'npeaks':  # keep top n number of peaks
                        spec.keeptopn(self.kw['npeaks'])
                    elif self.kw['dropmethod'] == 'consolidate':  # consolidate values being dropped
                        spec.consolidate(self.kw['threshold'], 3 * 10 ** -self.kw['consolidate'])
            else:  # if specific isotope

                ele, iso = self.isotope(key)  # find element and isotope
                if self.kw['verbose'] is True:
                    prog = Progress(string='Processing isotope %s' % (key), fraction=False, percent=False)
                if spec is None:  # if spectrum object has not been defined
                    spec = Spectrum(
                        self.kw['decpl'],  # decimal places
                        start=(self.md[ele][iso][0] * float(comp[key])) - 10 ** -self.kw['decpl'],  # minimum mass
                        end=(self.md[ele][iso][0] * float(comp[key])) + 10 ** -self.kw['decpl'],  # maximum mass
                        specin=[[self.md[ele][iso][0] * float(comp[key])], [1.]],
                        # supply masses and abundances as initialization spectrum
                        empty=self.kw['emptyspec'],  # whether or not to use emptyspec
                        filler=0.  # fill with zeros, not None
                    )
                    continue
                spec.shiftx(self.md[ele][iso][0])  # offset spectrum object by the mass of that
            if self.kw['verbose'] is True:
                prog.fin(' ')
        if self.kw['charge'] > 1:  # if the magnitude of the charge is greater than 1 divide m/z values
            spec.applycharge(self.kw['charge'])
        spec.normalize()
        if self.kw['keepall'] is False:
            spec.threshold(0.0001)  # drop very low intensity
        out = spec.trim()  # trim to x,y lists
        if self.kw['verbose'] is True:
            self.sys.stdout.write('DONE\n')
        return out

    def load_from_pickle(self, customfile=None):
        """loads data from pickle"""
        # TODO specify hierachy and pull if better method than specified
        import os
        import pickle
        if customfile is None:  # if no directory was specified, use current working directory
            customfile = os.path.join(
                os.getcwd(),
                'molecules',
                self.molecular_formula() + '.mol',
            )
        if os.path.isfile(customfile) is True:
            if self.kw['ipmethod'] in ['multiplicative', ]:
                key = 'multiplicative'
            elif self.kw['ipmethod'] in ['combinatorics', ]:
                key = 'combinatorics'
            if self.kw['dropmethod'] is not None:
                key += ' %s' % self.kw['dropmethod']
            subkey = self.kw['decpl']  # decimal places
            with open(customfile, 'rb') as targetfile:
                incoming = pickle.load(targetfile)
                if key in incoming and subkey in incoming[key]:
                    items = incoming[key][subkey]
                    strcharge = '%s%d' % (self.kw['sign'], self.kw['charge'])
                    if items['charge'] == strcharge:  # if the charge combination matches
                        print('Loading data from saved file %s' % customfile)
                        self.barip = items['bar isotope pattern']
                        self.rawip = items['raw isotope pattern']
                        self.gausip = items['gaussian isotope pattern']
                        self.mw = items['mw']
                        self.mimass = items['monoisotopic mass']
                        self.em = items['estimated exact mass']
                        self.pcomp = items['percent composition']
                        self.error = items['error']
                        self.fwhm = items['full width at half max']
                        self.sigma = items['standard deviation']
                        self.sf = self.molecular_formula()
                        return True
        return False  # if the exact match was not found, False

    def molecular_formula(self):
        """
        generates the molecular formula as the string
        The output will be ordered with carbon, then hydrogen, then
        the remaining elements alphabetically.
        """
        if self.kw['verbose'] is True:
            self.sys.stdout.write('Molecular formula: ')
        out = ''
        if 'C' in self.comp:  # carbon and hydrogen first according to hill formula
            out += 'C'
            if self.comp['C'] > 1:
                out += str(self.comp['C'])
        if 'H' in self.comp:
            out += 'H'
            if self.comp['H'] > 1:
                out += str(self.comp['H'])
        for key, val in sorted(self.comp.items()):  # alphabetically otherwise
            if key != 'C' and key != 'H':
                if key in self.md:
                    out += key
                    if self.comp[key] > 1:
                        out += str(self.comp[key])
                else:  # if an isotope
                    ele, iso = self.isotope(key)
                    out += '(%d%s)' % (iso, ele)
                    if self.comp[key] > 1:
                        out += str(self.comp[key])
        if self.kw['verbose'] is True:
            self.sys.stdout.write('%s DONE\n' % out)
        return out

    def molecular_weight(self):
        """determines the molecular weight from natural abundances"""
        if self.kw['verbose'] is True:
            self.sys.stdout.write('Calculating molecular weight')
        mwout = 0
        pcompout = {}  # percent composition dictionary
        for element in self.comp:
            try:
                for isotope in self.md[element]:
                    if isotope == 0:
                        continue
                    mwout += self.md[element][isotope][0] * self.md[element][isotope][1] * self.comp[
                        element]  # add every isotope times its natural abundance times the number of that element
                    if element not in pcompout:
                        pcompout[element] = 0
                    pcompout[element] += self.md[element][isotope][0] * self.md[element][isotope][1] * self.comp[
                        element]  # add mass contributed by that element
            except KeyError:  # if isotope
                ele, iso = self.isotope(element)
                mwout += self.md[ele][iso][0] * self.comp[element]  # assumes 100% abundance if specified
                pcompout[str(iso) + ele] = self.md[ele][iso][0] * self.comp[element]
        for element in pcompout:  # determines the percent composition of each element
            try:
                pcompout[element] = pcompout[element] / mwout
            except ZeroDivisionError:
                pcompout[element] = 0.
        if self.kw['verbose'] is True:
            self.sys.stdout.write(': %.6f DONE\n' % mwout)
        return mwout, pcompout

    def monoisotopic_mass(self):
        """
        a quick estimation of the exact mass given a molecular formula
        This may not be accurate for high mass species
        """
        if self.kw['verbose'] is True:
            self.sys.stdout.write('Estimating exact mass: ')
        em = 0.
        for key in self.comp:
            try:
                em += self.md[key][0][0] * self.comp[key]
            except KeyError:
                ele, iso = self.isotope(key)
                em += self.md[ele][iso][0] * self.comp[key]
        # # accounts for the mass of an electron (uncomment if this affects your data)
        # if self.kw['sign'] == '+':
        #    em -= (9.10938356*10**-28)*charge
        # if self.kw['sign'] == '-':
        #    em += (9.10938356*10**-28)*charge
        if self.kw['verbose'] is True:
            self.sys.stdout.write('%.5f DONE\n' % (em / self.kw['charge']))
        return em / self.kw['charge']

    def nominal_mass(self):
        """calculates the nominal mass of the molecule"""
        return int(round(self.em))

    def print_details(self):
        """prints the details of the generated molecule"""
        if 'sys' not in self.__dict__:
            self.sys = __import__('sys')
        self.sys.stdout.write('%s\n' % self)
        self.sys.stdout.write('formula: %s\n' % self.sf)
        self.sys.stdout.write('molecular weight: %.6f\n' % round(self.mw, 6))
        self.sys.stdout.write('monoisotopic mass: %.5f\n' % round(self.mimass, 5))
        self.sys.stdout.write('estimated exact mass: %.5f\n' % round(self.em, 5))
        self.sys.stdout.write('error: %.2g\n' % self.error)
        if abs(self.error) > self.kw['criticalerror']:
            self.sys.stdout.write('WARNING: Error is greater than %.2g!\n' % self.kw['criticalerror'])
        self.print_percent_composition()

    def print_percent_composition(self):
        """prints the percent composition in a reader-friendly format"""
        if 'sys' not in self.__dict__:
            self.sys = __import__('sys')
        self.sys.stdout.write('\nelemental percent composition:\n')
        for key, val in sorted(self.pcomp.items()):
            self.sys.stdout.write('%3s: %7.3f %%\n' % (key, self.pcomp[key] * 100))

    def plot_bar_pattern(self):
        """plots and shows the isotope bar pattern"""
        import pylab as pl
        fwhm = self.em / self.kw['res']
        pl.bar(self.barip[0], self.barip[1], width=fwhm, align='center')
        pl.xlabel('m/z', style='italic')
        pl.ylabel('normalized intensity')
        pl.ticklabel_format(useOffset=False)
        pl.show()

    def plot_gaussian_pattern(self, exp=None):
        """plots and shows the simulated gaussian isotope pattern"""
        import pylab as pl
        try:
            pl.plot(self.gausip[0], self.gausip[1], linewidth=1)
        except (AttributeError, TypeError):
            self.gausip = self.gaussian_isotope_pattern(self.rawip)
            pl.plot(self.gausip[0], self.gausip[1], linewidth=1)
        if exp is not None:  # plots experimental if supplied
            y = []
            maxy = max(exp[1])
            for val in exp[1]:  # normalize
                y.append(val / maxy * 100)
            comp = self.compare(exp)
            pl.plot(exp[0], exp[1])
            pl.text(max(exp[0]), 95, 'SER: ' + str(comp))
            # pl.fill_between(x,self.gausip[1],exp[1],where= exp[1] =< self.gausip[1],interpolate=True, facecolor='red')
        pl.fill(self.gausip[0], self.gausip[1], alpha=0.25)  # ,facecolor='blue')
        pl.xlabel('m/z', style='italic')
        pl.ylabel('normalized intensity')
        pl.ticklabel_format(useOffset=False)
        pl.show()

    def plot_raw_pattern(self):
        """plots and shows the raw isotope pattern (with mass defects preserved)"""
        import pylab as pl
        pl.bar(self.rawip[0], self.rawip[1], width=self.fwhm)
        pl.xlabel('m/z', style='italic')
        pl.ylabel('normalized intensity')
        pl.ticklabel_format(useOffset=False)
        pl.show()

    def reset(self):
        """resets values to when the instance was created"""
        self.__dict__ = self.original

    def save_to_jcamp(self, name=None):
        """
        Saves the bar isotope pattern to JCAMP-DX file format
        Output type roughly based on the output from ChemCalc.org
        see http://www.jcamp-dx.org/protocols.html for details on the JCAMP-DX specifications.

        :param name: optional name for the output file (default is {supplied formula} ({molecular formula}).jdx)
        """
        mf = self.molecular_formula()
        if name is None:  # if no name supplied, auto generate
            name = self.formula
            if self.formula != mf:
                name += ' (%s)' % mf
            name += '.jdx'
        elif name.lower().endswith('.jdx') is False:
            name += '.jdx'
        from datetime import datetime

        header = [  # comment lines to put before data
            # header items
            'TITLE= %s (%s)' % (self.formula, mf) if self.formula != mf else 'TITLE= %s' % self.formula,
            'JCAMP-DX= 5.01',
            'DATA TYPE= MASS SPECTRUM',
            'DATA CLASS= PEAK TABLE',
            'ORIGIN= Calculated spectrum from PythoMS Molecule class https://github.com/larsyunker/PythoMS',
            'OWNER= User of PythoMS',
            'SPECTROMETER/DATA SYSTEM= Molecule class %s method' % self.kw['ipmethod'],
            'CREATION DATE= %s' % datetime.now().astimezone(),
            'XUNITS= M/Z',
            'YUNITS= RELATIVE ABUNDANCE',
            'NPOINTS= %d' % len(self.barip[0]),
            'FIRSTX= %f' % self.barip[0][0],
            'LASTX= %f' % self.barip[0][-1],

            # user defined labels
            '$Molecular weight= %f' % self.molecular_weight()[0],
            '$Resolution= %f' % self.kw['res'],
            '$Threshold= %f' % self.kw['threshold'] if self.kw['threshold'] is not None else '',
            '$Error= %.2g' % self.error,
            '$Nominal mass = %d' % self.nominal_mass(),
            '$Monoisotopic mass= %f' % self.mimass,
            '$Estimated exact mass= %f' % self.em,
        ]
        with open(name, 'wt') as outfile:
            for line in header:  # write header lines
                if len(line) != 0:
                    outfile.write('##%s\n' % line)
            outfile.write('##PEAK TABLE= (XY..XY)\n')
            for ind, val in enumerate(self.barip[0]):  # write data lines
                outfile.write('%f, %f\n' % (val, self.barip[1][ind]))
            outfile.write('##END=\n')

    def save_to_pickle(self, customdir=None):
        """
        Saves the molecule's properties to pickle
        """
        def buildout():
            """retrieves items from self"""
            nonlocal self
            return {
                'bar isotope pattern': self.barip,
                'raw isotope pattern': self.rawip,
                'gaussian isotope pattern': self.gausip if 'gausip' in self.__dict__ else None,
                'mw': self.mw,
                'monoisotopic mass': self.mimass,
                'estimated exact mass': self.em,
                # 'nominal mass': self.nominal_mass(),
                'percent composition': self.pcomp,
                'error': self.error,
                'charge': '%s%d' % (self.kw['sign'], self.kw['charge']),
                'full width at half max': self.fwhm,
                'standard deviation': self.sigma,
            }

        import pickle
        import os
        if customdir is None:  # if no directory was specified, use current working directory
            customdir = os.path.join(
                os.getcwd(),
                'molecules'
            )
        if os.path.isdir(customdir) is False:  # if directory does not exist
            if input('The directory "%s" does not exist, create it? (Y/N) ' % customdir).lower() in ['y','yes']:
                os.makedirs(customdir)
            else:
                print('Directory not created, data was not saved.')
                return None

        fp = os.path.join(
            customdir,
            self.molecular_formula() + '.mol',
        )
        if self.kw['ipmethod'] in ['multiplicative', ]:
            key = 'multiplicative'
        elif self.kw['ipmethod'] in ['combinatorics', ]:
            key = 'combinatorics'
        if self.kw['dropmethod'] is not None:
            key += ' %s' % self.kw['dropmethod']
        subkey = self.kw['decpl']  # decimal places
        if os.path.isfile(fp):
            with open(fp, 'r+b') as targetfile:
                out = pickle.load(targetfile)  # load existing values
                thisinstance = buildout()
                if key in out and subkey in out[key]:
                    if out[key][subkey] == thisinstance:  # if the output is equal to
                        return None
                    out[key][subkey] = buildout()
                pickle.dump(out, targetfile)
        with open(fp, 'wb') as targetfile:
            pickle.dump(
                {key: {subkey: buildout()}},
                targetfile
            )

    def standard_deviation(self, fwhm):
        """determines the standard deviation for a normal distribution with the full width at half max specified"""
        sigma = fwhm / (2 * self.np.sqrt(2 * self.np.log(2)))  # based on the equation FWHM = 2*sqrt(2ln2)*sigma
        return sigma

    def standard_deviation_comp(self, comp):
        """
        cacluates the standard deviation of the isotope pattern of the supplied composition
        this calculation is based on Rockwood and Van Orden 1996 doi: 10.1021/ac951158i
        """
        stdev = 0
        for key in comp:
            meanmass = 0
            for mass in self.md[key]:
                if mass != 0:
                    meanmass += self.md[key][mass][1] * self.md[key][mass][0]  # weighted average mass
            eledev = 0  # elemental deviation
            for mass in self.md[key]:
                if mass != 0:
                    eledev += self.md[key][mass][1] * (self.md[key][mass][0] - meanmass) ** 2
            stdev += eledev * comp[key]
        return self.np.sqrt(stdev)

    def standard_deviation_ip(self, ip):
        """
        cacluates the standard deviation of the isotope pattern of the supplied composition
        this calculation is based on Rockwood and Van Orden 1996 doi: 10.1021/ac951158i
        """
        import scipy as sp
        stdev = 0
        for ind, val in enumerate(ip[0]):
            stdev += ip[1][ind] * (val - self.pmw) ** 2  # weighted distance from the estimated molecular weight
        return sp.sqrt(stdev)

    def validate_pattern(self, mw, pattern):
        """
        Validates a computed isotope pattern by comparing the molecular weight from 
        the computation to the molecular weight given by the formula and average masses. 
        
        **Parameters**
        
        mw: *float*
            The molecular weight given by the molecular formula and average masses. 
        
        pattern: *list*
            The computed isotope pattern. 
        
        
        **Returns**
        
        relative difference: *float*
            The relative difference of the calculated molecular formula to the 
            actual molecular formula. 
            Typically a difference of 3 parts per million (3*10^-6) is deemed acceptable 
            error. 
        
        """
        if self.kw['verbose'] is True:
            self.sys.stdout.write('Estimated molecular weight from isotope pattern: ')
        self.pmw = 0  # pattern molecular weight
        for ind, mz in enumerate(pattern[0]):
            self.pmw += mz * pattern[1][ind] * self.kw['charge']
        self.pmw /= sum(pattern[1])
        if self.kw['verbose'] is True:
            self.sys.stdout.write('%f DONE\n' % (self.pmw))
        return (self.pmw - mw) / mw

    def weighted_average(self, ipgroup):
        """
        Determines the weighted average of a group of x and y values.

        **Parameters**

        ipgroup: *list*
            Paired list of ``[[m/z values],[intensity valules]]`` to be grouped.


        **Returns**

        weighted average: *list*
            Returns the weighted average location along the x axis, and the
            sum of the y values in ``[x,y]`` format.

        """
        if sum(ipgroup[1]) == 0:  # catch for no intensity
            return sum(ipgroup[0]) / len(ipgroup[0]), 0.
        s = 0
        for ind, val in enumerate(ipgroup[0]):  # sum mz*int pairs
            s += val * ipgroup[1][ind]
        return s / sum(ipgroup[1]), sum(ipgroup[1])  # return weighted m/z, summed intensity


if __name__ == '__main__':  # for testing and troubleshooting
    # st.printstart()
    mol = Molecule(
        # 'C3900H4902N1500O2401P400',
        # 'W100',
        'L2PdAr+I',
        # charge= 2, # specify charge (if not specified in formula)
        # res=1050000, # specify spectrometer resolution (default 5000)
        verbose=True,
        # decpl=10,
        # dropmethod='threshold',
        # threshold=0.00001,
        # ipmethod='hybrid',
        # ipmethod='combinatorics',
        # keepall=True,
    )
    # mol.print_details()
    # st.printelapsed()
    # st.printprofiles()
