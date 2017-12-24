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

IGNORE:
CHANGELOG
---2.5---
- added the ability to not provide start and end points for an unfilled spectrum
---2.6
IGNORE
"""

from PythoMS._classes._ScriptTime import ScriptTime

st = ScriptTime(profile=True)


class Spectrum(object):
    def __init__(self, decpl, **kwargs):
        """
        A class for combining and otherwise manipulating spectra with non-equal dimensions.
        The object will track x values to a specified decimal place and can efficiently
        add a new value to a growing list of values.
        e.g. adding two spectra together that do not have an intensity value for every
        x value (common for mass spectra).

        **Parameters**

        decpl: *integer*
            The decimal places to track the x values two. e.g. a value of 3 will track
            x values to the nearest 0.001.


        **\*\*kwargs**

        start: 50.
            The minimum x value to track. Attempts to add an x value less than this
            will be ignored by the object. Options: float.

        end: 2000.
            The maximum x value to track. Attempts to add an x value greater than
            this will be ignored by the object. Options: float.

        specin: None
            An spectrum to be added to the object on initialization. The format should be
            ``[[x values],[y values]]``.
            Options: None, list

        empty: False
            Whether the spectrum object should be filled or empty. Options: bool.
            An empty spectrum will have no x or y values on initialization, and will
            add values with each call of ``Spectrum.addvalue()``.
            A filled spectrum will generate an x list from *start* to *end* with spacing
            10^-*decpl* and a y list of equal length filled with the value specified by
            *filler*.
            If the number of items to be contained in the spectrum is substantially less
            than ``(end-start)*10^decpl`` it can be more efficient to set this to True.
            If not, then set this to False to improve computation times.

        filler: None
            The y value to use if there is no y value. This can affect the functionality
            of some of the functions in this class. If ``Spectrum.addelement()`` is to be
            used (e.g. by the Molecule class), filler must be ``0.``.

        reusable: False
            If the same y list is to be reset and reused multiple times, set this to True.
            The method resety() can then be called to reset the y list without having to
            recreate the Spectrum instance.


        **Examples**

        ::

            >>> spec = Spectrum(3)
            >>> spec.addvalue(55.67839,100)
            >>> spec.trim()
            [[55.678], [100]]

            >>> spec.addvalue(55.67799,100)
            >>> spec.trim()
            [[55.678], [200]]


        """
        self.kw = {  # default keyword arguments
            'start': 50.,  # start m/z for the object (float)
            'end': 2000.,  # end m/z for the object (float)
            'specin': None,  # supplied spectrum on initialization
            'empty': False,  # whether the spectrum object should be filled or empty
            'filler': None,  # thing to fill the spectrum with
            'reusable': False,  # whether the same y list is to be reused
            # 'profile': False, # profile the functions
        }
        if set(kwargs.keys()) - set(self.kw.keys()):  # check for invalid keyword arguments
            string = ''
            for i in set(kwargs.keys()) - set(self.kw.keys()):
                string += ' %s' % str(i)
            raise KeyError('Unsupported keyword argument(s): %s' % string)
        self.kw.update(kwargs)  # update defaules with provided keyword arguments

        self.decpl = decpl
        self.sp = __import__('scipy')
        if self.kw['start'] is None:  # if no start is specified
            if self.kw['empty'] is False:
                raise ValueError('A start value must be specified for a filled Spectrum object')
            self.start = -self.sp.inf
        else:
            self.start = round(self.kw['start'], self.decpl)
        if self.kw['end'] is None:  # if no end is specified
            if self.kw['empty'] is False:
                raise ValueError('An end value must be specified for a filled Spectrum object')
            self.end = self.sp.inf
        else:
            self.end = round(self.kw['end'], self.decpl)
        self.filler = self.kw['filler']
        if self.kw['empty'] is True:
            from bisect import bisect_left
            self.bl = bisect_left
        # if self.kw['profile'] is True:
        #    from _ScriptTime import ScriptTime
        #    self.st = ScriptTime(profile=True)
        #    self.st.printstart()
        specin = self.kw.pop('specin')
        if self.kw['reusable'] is True:
            self.x, self.yimm = self.fullspeclist(self.start,
                                                  self.end)  # m/z and intensity lists (these are intended to remain immutable)
            self.y = list(self.yimm)  # create the list that will be actively modified
        else:
            self.x, self.y = self.fullspeclist(self.start, self.end)  # m/z and intensity lists

        if specin is not None:
            self.addspectrum(specin[0], specin[1])

    def __str__(self):
        return 'Full spectrum list from {} to {} keeping {} decimal places'.format(self.start, self.end, self.decpl)

    def __repr__(self):
        return '{}(decpl={},start={},end={})'.format(self.__class__.__name__, self.decpl, self.start, self.end)

    def __len__(self):
        return len(self.x)

    def __getitem__(self, ind):
        """
        if supplied index is an integer, return the x and y value of that index in the list
        if a float, return the intensity of that x value
        """
        if type(ind) is int:
            return [self.x[ind], self.y[ind]]
        elif type(ind) is float:  # returns the intensity value of the specified m/z
            if ind < self.start or ind > self.end:
                raise IndexError(
                    'The supplied float %f is outside of the m/z range of this Spectrum instance (%.3f -%.3f)' % (
                    ind, self.start, self.end))
            return self.y[self.index(ind)]

    def __add__(self, x):
        """
        Since addition to this class requires generating a complete copy of the class then addition,
        using the built in addition methods is recommended
        e.g. to add a single value use .addvalue()
        to add a spectrum use .addspectrum()
        """
        if isinstance(x, self.__class__) is True:  # if it is another Spectrum instance
            if x.decpl != self.decpl:
                raise ValueError(
                    'The decimal places of the two spectra to be added are not equal. Addition is not supported')
            newstart = min(self.start, x.start)  # find new start m/z
            newend = max(self.end, x.end)  # find new end m/z
            tempspec = Spectrum(self.decpl, start=newstart, end=newend, specin=[self.x, self.y], empty=self.kw['empty'],
                                filler=self.filler)  # temporary instance with self as specin
            tempspec.addspectrum(x.x, x.y)  # add incoming spectrum
            return tempspec
        elif type(x) is int:  # add this integer to every m/z
            tempspec = Spectrum(self.decpl, **self.kw)
            y = self.sp.asarray(self.y)
            y += x
            tempspec.y = y.tolist()
            return tempspec
        elif len(x) == 2 and len(x[0]) == len(x[1]):  # if it is a list of paired lists (another spectrum)
            tempspec = Spectrum(self.decpl, **self.kw)
            tempspec.y = list(self.y)
            tempspec.addspectrum(x[0], x[1])
            return tempspec
        else:
            return 'Addition of %s to the Spectrum class is unsupported' % str(x)

    def __sub__(self, x):
        if isinstance(x, self.__class__) is True:  # if it is another Spectrum instance
            if x.decpl != self.decpl:
                raise ValueError(
                    'The decimal places of the two spectra to be added are not equal. Subtraction is not supported')
            newstart = min(self.start, x.start)  # find new start m/z
            newend = max(self.end, x.end)  # find new end m/z
            tempspec = Spectrum(self.decpl, start=newstart, end=newend, specin=[self.x, self.y], empty=self.kw['empty'],
                                filler=self.filler)  # temporary instance
            tempspec.addspectrum(x.x, x.y, True)  # subtract incoming spectrum
            return tempspec
        elif type(x) is int:  # add this integer to every m/z
            tempspec = Spectrum(self.decpl, **self.kw)
            y = self.sp.asarray(self.y)
            y -= x
            tempspec.y = y.tolist()
            return tempspec
        elif len(x) == 2 and len(x[0]) == len(x[1]):  # if it is a list of paired lists (another spectrum)
            tempspec = Spectrum(self.decpl, **self.kw)
            tempspec.y = list(self.y)
            tempspec.addspectrum(x[0], x[1], True)  # subtract the incoming spectrum
            return tempspec
        else:
            return 'Subtraction of %s from the Spectrum class is unsupported' % str(x)

    def __mul__(self, x):
        raise AttributeError('Multiplication of the Spectrum class is unsupported')

    def __div__(self, x):
        raise AttributeError('Division of the Spectrum class is unsupported')

    def __pow__(self, x):
        raise AttributeError('Raising a Spectrum instance to a power is unsupported.\nAlso... really?!')

    @st.profilefn
    def addelement(self, masses, abunds):
        """
        Adds the masses and abundances of an element to the current spectrum object.
        This is more efficient than creating a new spectrum object every time an
        element is added.

        **Parameters**

        masses: *list*
            List of masses. e.g. for carbon
            ``masses = [12.0,13.0033548378]``

        abunds: *list*
            List of abundances, paired with *masses*. e.g. for carbon
            ``abunds = [0.9893, 0.0107]``

        **Note**

        This function will encounter an error if there are None values in the y list.
        If you intend to use this function, set the *filler* keyword argument to be
        some value that is not None (i.e. ``0.``).

        """
        if len(masses) != len(abunds):
            raise ValueError(
                'The dimensions of the supplied lists are not equal (length masses: %d, length abundances: %d)' % (
                len(masses), len(abunds)))
        # xboxed = [[0]]
        # yboxed = [[1]]
        xboxed = []
        yboxed = []
        for ind, val in enumerate(masses):  # box up the lists
            xboxed.append([val])
            yboxed.append([abunds[ind]])
        newx = self.sp.asarray(self.x) + self.sp.asarray(xboxed)  # matrix of new x values
        newy = self.sp.asarray(self.y) * self.sp.asarray(yboxed)  # matrix of new y values

        ## does not call for a new disposable object
        """
        extending the spectrum object then dropping is very fast,
        but requires subtracting the original spectrum before dropping
        slicing, deepcopy, list(), building the subtraction into the boxxed lists, and switching the array calls
        are all slower than generating a temporary Spectrum object
        """
        # self.newend(max(masses) + max(self.x)) # define new end point for Spectrum object
        # for i in range(newx.shape[0]): # add calculated masses and intensities to object
        #    if i == 0:
        #        self.addspectrum(newx[i],newy[i],True) # subtract original spectrum
        #        continue
        #    self.addspectrum(newx[i],newy[i])
        ##self.addspectrum(oldx,oldy,True) # subtract old spectrum
        # self.dropbelow(min(masses) + min(self.x)) # drop values below new start point

        self.kw['start'] = min(masses) + self.start - 10 ** -self.decpl  # bounds for new Spectrum object
        self.start = self.kw['start']
        self.kw['end'] = max(masses) + self.end + 10 ** -self.decpl
        self.end = self.kw['end']
        tempspec = Spectrum(self.decpl, **self.kw)
        for i in range(newx.shape[0]):
            tempspec.addspectrum(newx[i], newy[i])
        self.x = tempspec.x  # redefine the x and y lists
        self.y = tempspec.y

    # @st.profilefn
    def addvalue(self, xval, yval, subtract=False):
        """
        Adds an intensity value to the x value specified.

        **Parameters**

        xval: *float*
            The x value of the intensity. This value will be rounded to the
            decimal place specified on calling this class.

        yval: *float*
            The y value to add to the y list.

        subtract: *bool*
            Make this True if you wish to subtract the y value from the
            current y list.


        **Returns**

        return item: *type*
            description


        **Examples**

        ::

            >>> spec = Spectrum(3)
            >>> spec.trim()
            [[], []]
            >>> spec.addvalue(673.9082342357,100)
            >>> spec.trim()
            [[673.908], [100]]
            >>> spec.addvalue(1523.25375621,200)
            >>> spec.addvalue(50.89123,300)
            >>> spec.trim()
            [[50.891, 673.908, 1523.254], [300, 100, 200]]


        **Note**

        If the x value is not within the x bounds specified by the keyword
        arguments *start* and *end*, the supplied y value will not be added
        to the current spectrum object.

        """
        if yval is not None:  # if handed an actual value
            try:  # try indexing
                index = self.index(xval)
                if subtract is True:  # set sign based on input
                    sign = -1
                else:
                    sign = 1
                if self.kw['empty'] is False:  # if x list filled
                    try:
                        self.y[index] += yval * sign  # try to add value
                    except TypeError:
                        self.y[index] = yval * sign  # if None, then set to value
                else:
                    if len(self.x) == 0 and index == 0:
                        self.x.insert(index, round(xval, self.decpl))
                        self.y.insert(index, yval * sign)
                    elif index == len(self.x):  # if at end of list
                        self.x.append(round(xval, self.decpl))
                        self.y.append(yval)
                    elif self.x[index] != round(xval, self.decpl):  # if the index does not equal the value
                        self.x.insert(index, round(xval, self.decpl))  # insert x value at specified index
                        self.y.insert(index, yval * sign)  # insert the y value
                    else:
                        try:  # otherwise add
                            self.y[index] += yval * sign
                        except TypeError:  # or set to value if None
                            self.y[index] = yval * sign
            except ValueError:  # if index is not in spectrum
                pass  # do nothing (the value will not be added to the spectrum)

    @st.profilefn
    def addspectrum(self, x, y, subtract=False):
        """
        Adds an entire x and y list to the spectrum object.
        This avoids having to call ``Spectrum.addvalue()`` in loop form
        for a list of values.

        **Parameters**

        x: *list*
            List of x values. These do not need to be sorted, but are assumed
            to be paired with the supplied *y* list.

        y: *list*
            List of y values, paired with *x*.

        subtract: *bool*
            Whether or not to subtract the y intensities from the current
            spectrum object.


        """
        if len(x) != len(y):
            raise ValueError('The addspectrum() method only supports two lists of the same dimension')
        for ind, mz in enumerate(x):
            if y[ind] != self.filler:  # drops filler values at this point
                self.addvalue(mz, y[ind], subtract)

    def checknone(self):
        """counts the number of not-None values in the current y list (for debugging)"""
        count = 0
        for val in self.y:
            if val is not None:
                count += 1
        return count

    def countnone(self):
        """counts the number of None values in the current y list (for debugging)"""
        count = 0
        for val in self.y:
            if val is None:
                count += 1
        return count

    # @st.profilefn
    def consolidate(self, threshold, within, method='abs'):
        """
        A method of reducing the number of values in the spectrum object by
        consolidating y values below the specified threshold with nearby values.
        The method of combination is a weighted average. The intensities of
        adjacent values are combined until the threshold is passed or until no
        adjacent values within the specified x delta can be found.

        **Parameters**

        threshold: *float*
            The threshold value, below which the value will be consolidated
            into adjacent peaks.

        within: *float*
            The x delta to look within when consolidating peaks.

        method: 'abs' or 'rel'
            Whether to use an absolute or relative threshold.


        """

        def adjacent(index):
            """
            locates the index of the closest x value to the provided index
            (only returns an index if there is a value within the given delta
            """
            i = None
            if index != 0 and self.x[index] - self.x[index - 1] <= within:  # if previous index is nearer than within
                if self.y[index - 1] != 0:  # if the peak has intensity
                    i = index - 1
                    delta = self.x[index] - self.x[index - 1]
            if index != len(self.x) - 1 and self.x[index + 1] - self.x[
                index] <= within:  # if next index is nearer than within
                if self.y[index + 1] != 0:  # if the peak has intensity
                    if i is not None:  # if i is already defined
                        if self.x[index + 1] - self.x[index] <= delta:  # if the greater than is closer
                            i = index + 1
                    else:
                        i = index + 1
            return i

        def weightedaverage(ipgroup):
            """
            Determines the weighted average of a group of masses and abundances
            """
            s = 0
            for ind, val in enumerate(ipgroup[0]):  # sum mz*int pairs
                s += val * ipgroup[1][ind]
            return s / sum(ipgroup[1]), sum(ipgroup[1])  # return weighted m/z, summed intensity

        for ind in range(len(self.y)):
            if self.y[ind] < threshold and self.y[ind] != 0.:
                cur = ind
                closest = adjacent(cur)  # looks for adjacent peaks within the delta
                while self.y[cur] < threshold and closest is not None:
                    wx, wy = weightedaverage([[self.x[cur], self.x[closest]], [self.x[cur], self.y[closest]]])
                    self.addvalue(self.x[cur], self.y[cur], True)  # subtract the current value
                    self.addvalue(self.x[closest], self.y[closest], True)  # subtract the adjacent value
                    self.addvalue(wx, wy)  # add the weighted average to the spectrum
                    cur = self.index(wx)  # set current index to that of the new value
                    # cur = closest # set current index to the one being tested
                    closest = adjacent(cur)
        self.threshold(threshold, method)  # drop any peaks that could not be combined

    def cp(self):
        """returns a list (clone) of the spectrum"""
        return [list(self.x), list(self.y)]

    def dropabove(self, value):
        """drops all values at higher x values than the specified value"""
        index = self.index(value)  # find index
        self.newend(value)  # redefine end value
        self.x = self.x[:index]  # trim lists
        self.y = self.y[:index]

    # @st.profilefn
    def dropbelow(self, value):
        """drops all values at lower x values than the specified value"""
        index = self.index(value)  # find index
        self.newstart(value)  # redefine start value
        del self.x[:index]
        del self.y[:index]

    def fullspeclist(self, start, end):
        """
        Generates two paired lists (one m/z, one None) from start to end with a
        specified number of decimal places
        """
        """
        Generates two paired lists (one x, one y). 
        
        **Parameters**
        
        start: *float*
            The start value for the x list. 
        
        end: *float*
            The end value for the x list. 
        
        
        **Returns**
        
        x and y lists: *list*
            A list of x values and a list of y values (specified by the *filler* 
            keyword argument) of the same length. 
        
        
        **Notes**
        
        The maximum x value will be larger than end by 10^-decpl to include the 
        actual end value in the x list. 
        
        """
        if self.kw['empty'] is False:
            x = self.sp.arange(start, end + 10 ** -self.decpl,
                               10 ** -self.decpl).tolist()  # generate x values using arange and convert to list
            y = [self.filler] * len(x)  # generate y list of equal length
        else:
            x = []
            y = []
            # x = [start,end]
            # y = [self.filler,self.filler]
        return x, y

    def fillzeros(self, value=0.):
        """Replaces any None values in the y list with the specified value."""
        for ind, inten in enumerate(self.y):
            if inten is None:
                self.y[ind] = value
        return self.y

    # @st.profilefn
    def index(self, xval):
        """
        Locates the index of the specified x value in the object's x list.

        **Parameters**

        xval: *float*
            The x value to locate in the list.


        **Returns**

        index: *integer*
            The integer index for the x value in the x list.


        **Notes**

        If the *empty* keyword argument is True, the index will be located
        using the bisect module. If the x value is not in the current x
        list, the appropriate insertion index is returned.
        If the *empty* keyword argument is False, the index will be calculated
        based on the *start* and *decpl* values of the object (this is more
        computationally efficient than bisection).
        """
        if xval > self.end or xval < self.start:
            raise ValueError(
                'the m/z value ({}) is outside of the m/z range of this spectrum ({}-{})'.format(xval, self.start,
                                                                                                 self.end))
        if self.kw['empty'] is True:  # if spectrum is unfilled, searching is required
            return self.bl(self.x, round(xval, self.decpl))
        if self.kw['empty'] is False:  # otherwise, calculation of the index is more efficient
            return int(round((xval - self.start) * (10 ** self.decpl)))  # rounds after multiplication

    def keeptopn(self, n=5000):
        """
        Keeps the top n peaks and sets the intensity of those below that
        value to be zero.

        **Parameters**

        n: *integer*
            The number of values to keep in the list.


        **Notes**

        If there is more than one y value equal to the nth lowest value,
        then all of those values will be retained.

        """
        sort = sorted(self.y, reverse=True)  # sorts y values
        try:
            critical = sort[n]  # set the critical value to the nth value
            self.threshold(critical, 'abs')  # apply the absolute threshold
        except IndexError:
            pass

    def max(self):
        """
        Locates the maximum intensity in the y list and returns the
        x and y values of that point

        **Returns**

        values: *list*
            Returns the x and y values of the maximum y value.

        """
        locs = self.sp.where(self.sp.asarray(self.y) == max(self.y))  # locate index of maximum
        if len(locs[0]) > 1:  # if there is more than one value equal to the maximum
            out = []
            for i in locs[0]:
                out.append([self.x[i], self.y[i]])
            return out
        return self.x[locs[0][0]], self.y[locs[0][0]]

    def newend(self, value):
        """
        Defines a new end value for the object.

        **Parameters**

        value: *float*
            The new end x value for the object.

        """
        if self.ks['empty'] is True:
            raise ValueError('Modifying the end value in this way is not supported if empty is True')
        self.end = value + 10 ** -self.decpl
        self.kw['end'] = self.end

    def newstart(self, value):
        """
        Defines a new start value for the object.

        **Parameters**

        value: *float*
            The new start x value for the object.

        """
        if self.ks['empty'] is True:
            raise ValueError('Modifying the start value in this way is not supported if empty is True')
        self.start = value - 10 ** -self.decpl
        self.kw['start'] = self.start

    # @st.profilefn
    def normalize(self, top=100.):
        """
        Normalizes the y values to the specified value.

        **Parameters**

        top: *float*
            The new value for the maximum y value.
        """
        m = max(self.y)
        self.y = self.sp.asarray(self.y)
        self.y /= m
        self.y *= top
        self.y = self.y.tolist()
        # for ind,inten in enumerate(self.y):
        #    if inten is not None:
        #        self.y[ind] = inten/m*top

    def resety(self):
        """
        Resets the y list. The reusable kwarg must be true to use this method.
        """
        self.y = list(self.yimm)

    def shiftx(self, value):
        """
        Offsets all x values by a the specified value.

        **Parameters**

        value: *float*
            The amount to offset the x values by.
        """
        for ind, val in enumerate(self.x):
            self.x[ind] += value
        if self.start != -self.sp.inf:
            self.start += value
        if self.end != self.sp.inf:
            self.end += value

    def sum(self):
        """
        Calculates and returns the sum of all y values in the object.

        **Returns**

        sum: *float*
            The sum of all y values in the y list.
        """
        out = 0
        for val in self.y:
            if val is not None:
                out += val
        return out

    # @st.profilefn
    def threshold(self, thresh, method='abs'):
        """
        Removes all y values below the specified threshold value.

        **Parameters**

        thresh: *float*
            The threshold y value to drop below.

        method: 'abs' or 'rel'
            Whether the specifed *thresh* value is absolute or
            relative to the maximum y value.

        """
        if method == 'rel':  # if relative, calculate relative threshold
            thresh *= max(self.y)
        if self.kw['empty'] is True:  # removes values from the list
            x = []
            y = []
            for ind, inten in enumerate(self.y):
                if inten > thresh:
                    x.append(self.x[ind])
                    y.append(inten)
            self.x = x
            self.y = y
        if self.kw['empty'] is False:
            for ind, inten in enumerate(self.y):
                if inten < thresh:
                    self.y[ind] = self.filler  # sets the value to the filler value

    def trim(self, zeros=False, xbounds=None):
        """
        Trims x and y pairs that have None intensity and returns the trimmed list.
        This is the most efficient way of converting a Spectrum object to an x and y list.

        **Parameters**

        zeros: *bool*
            Specifies whether there should be zeros at the start and end values.
            This can be used to generate continuum spectra across the range [start,end].
            If there are non-zero intensity values at the start or end point, they will
            not be affected.

        xbounds: None or *list*
            This can specify a subsection of the x and y spectra to trim to. None will
            return the entire contents of the Spectrum object, and specifying
            ``[x1,x2]]`` will return the x and y lists between *x1* and *x2*.

        """
        if xbounds is None:
            xbounds = [self.start, self.end]
        elif xbounds[0] is None:
            xbounds[0] = self.start
        elif xbounds[1] is None:
            xbounds[1] = self.end
        xout = []
        yout = []
        for ind, inten in enumerate(self.y):
            if xbounds[0] <= self.x[ind] <= xbounds[1]:  # if within the x bounds
                if inten is not None:
                    xout.append(round(self.x[ind], self.decpl))  # rounded to avoid array floating point weirdness
                    yout.append(inten)
                # elif zeros is True: # if zeros at the edges of spectrum are desired
                #    if self.x[ind] == xbounds[0] or self.x[ind] == xbounds[1]: # at the edges of the output spectrum
                #        xout.append(self.x[ind])
                #        yout.append(0)
        if zeros is True:
            if xout[0] != self.start:
                xout.insert(0, self.start)
                yout.insert(0, 0.)
            if xout[-1] != self.end:
                xout.append(self.end)
                yout.append(0.)
        return [xout, yout]


def checkindexing(n=1000, dec=3):
    """validates the indexing function of the Spectrum object"""
    from random import random
    mismatch = 0
    mml = []
    for i in range(1000):
        num = random()
        mz = num * 2000.
        try:
            index = spec.index(mz)
        except ValueError:
            continue
        if round(mz, dec) != round(spec.x[index], dec):
            mismatch += 1
            mml.append([mz, round(mz, dec), round(spec.x[index], dec)])
    return mismatch, mml


if __name__ == '__main__':
    spec = Spectrum(3)
    # spec = Spectrum(4,start=12.0,end=13.0033548378,specin=[[12.0,13.0033548378],[0.9893, 0.0107]],empty=True,filler=0.)
    # masses = [12.0,13.0033548378]
    # abunds = [0.9893, 0.0107]

    # spec = Spectrum(10,start=1,end=3,specin=[[1.00782503207,2.0141017778],[0.999885,0.000115]],empty=True,filler=0.)
    # masses = [1.00782503207,2.0141017778]
    # abunds = [0.999885,0.000115]

    # dec = 4
    # spec = Spectrum(dec,start=1,end=3,specin=[[1.00782503207,2.0141017778],[0.999885,0.000115]],empty=True,filler=0.)
    # masses = [15.99491461956,16.9991317,17.999161]
    # abunds = [0.99757,0.00038,0.00205]
    #
    #
    # print spec
    #
    # thresh = 0.01
    # cons = 3*10**-dec
    #
    # import sys
    # for i in range(3900):
    #    sys.stdout.write('\rcarbons: %d' %(i+1))
    #    spec.addelement([12.0,13.0033548378],[0.9893, 0.0107])
    #    spec.normalize(100.)
    #    spec.consolidate(thresh,cons)
    # sys.stdout.write('\n')
    # print 'length of x:', len(spec.x)
    # for i in range(2401):
    #    sys.stdout.write('\roxygens: %d' %(i+1))
    #    spec.addelement(masses,abunds)
    #    spec.normalize(100.)
    #    spec.consolidate(thresh,cons)
    # sys.stdout.write('\n')
    # print 'length of x:',len(spec.x)
    #
    # st.printelapsed()
    # st.printprofiles()
    # mismatch,mml = checkindexing(1000,3)
    # print mismatch
