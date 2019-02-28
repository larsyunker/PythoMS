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
- added applycharge function to apply the charge to a mass list
IGNORE
"""
import numpy as np
from random import random
from bisect import bisect_left as bl


def weighted_average(xvals, yvals):
    """
    Determines the weighted average of a group of masses and abundances

    :param list xvals: x values
    :param list yvals: y values
    :return: weighted average, summed intensity
    :rtype: tuple of float
    """
    if sum(yvals) == 0:  # catch for no intensity
        return sum(xvals) / len(xvals), 0.
    return (
        sum([x * y for x, y in zip(xvals, yvals)]) / sum(yvals),  # weighted m/z
        sum(yvals)  # summed intensity
    )


def full_spectrum_list(start, end, decpl, filler=None):
    """
    Generates two paired lists (one m/z, one None) from start to end with a specified number of decimal places.

    :param float start: The start value for the x list.
    :param float end: The end value for the x list.
    :param int decpl: The decimal places to use for the generated list.
    :param filler: The filler value for the y list
    :return: A list of x values and a list of y values (specified by the ``filler`` keyword argument) of the
        same length.
    :rtype: tuple of lists

    **Notes**

    The maximum x value will be larger than end by 10^-``decpl`` to include the
    actual end value in the x list.

    """
    x = np.arange(  # generate x values
        start,
        end + 10 ** -decpl,
        10 ** -decpl
    )
    return (
        x.tolist(),  # convert to list
        [filler] * len(x),  # generate y list of equal length
    )


class Spectrum(object):
    _start = -np.inf
    _end = np.inf
    _charge = 1

    def __init__(self,
                 decpl,
                 start=50.,
                 end=2000.,
                 empty=False,
                 filler=None,
                 specin=None,
                 ):
        """
        A class for subtracting, combining, adding-to, and otherwise manipulating spectra with non-equal dimensions.
        The object will track *x* values to a specified decimal place and can efficiently add a new value to a growing
        list of values. e.g. adding two spectra together that do not have an intensity value for every *x* value
        (a common operation for combining mass spectra). On initialization, specify the number of decimal places to
        track using the ``decpl`` argument. Other behaviour of the class can be tweaked with the keyword arguments.

        :param int decpl: The decimal places to track the *x* values two. e.g. a value of 3 will track
            *x* values to the nearest 0.001.
        :param float,None start: The minimum *x* value to track. Attempts to add an *x* value less than this
            will be ignored by the instance.
        :param float,None end: The maximum *x* value to track. Attempts to add an *x* value greater than this will be
            ignored by the instance.
        :param list specin: An spectrum to be added to the object on initialization. The format should be
            ``[[x values],[y values]]``.
        :param bool empty: Whether the spectrum object should be filled or empty. An empty spectrum will have no *x*
            or *y* values on initialization, and will add values with each call of ``Spectrum.addvalue()``. A filled
            spectrum will generate an *x* list from *start* to *end* with spacing 10^-``decpl`` and a *y* list of equal
            length filled with the value specified by the ``filler`` kwarg. If the number of items to be contained in
            the spectrum is substantially less than ``(end-start)*10^decpl`` it can be more efficient to set this to
            ``True``. If not, then set this to False to reduce computational overhead.
        :param filler: The y value to use if there is no y value. This can affect the functionality of some of the
            functions in this class. If ``Spectrum.addelement()`` is to be used (e.g. by the Molecule class),
            filler must be ``0.``.

        **Basic Examples**

        Specify the number of decimal places to track on initialization.

        >>> spec = Spectrum(3)

        *x*, *y* pairs may be added using the ``add_value`` method

        >>> spec.add_value(55.67839, 100)


        When the spectrum has been manipulated to the user's satisfaction, it may be easily converted to
        ``[[x values], [y values]`` format using the ``trim()`` method.

        >>> spec.trim()
        [[55.678], [100]]

        The incoming x value will be compared to the current x list for equivalent x values. If a matching x value is
        found, the y value is added to the existing value.

        >>> spec.add_value(55.67799, 100)  # equivalent to 55.678
        >>> spec.trim()
        [[55.678], [200]]
        >>> spec.add_value(55.67744, 99)  # equivalent to 55.677
        >>> spec.trim()
        [[55.677, 55.678], [99, 200]]


        **y-value manipulation**

        The y values may be manipulated in a variety of ways.

        - The ``normalize()`` method will normalize the y values in the instance to the specified value.
        - The ``threshold()`` method will drop y values below a certain value (either relative or absolute).
        - The ``keep_top_n()`` method keeps the top n peaks.
        - The ``consolidate()`` method groups values together using a weighted average algorithm to keep the lowest
            y value above a given threshold but still retain the information in the spectrum.

        **spectrum constraint methods**

        - Values below a certain x value may be dropped by calling the ``drop_below()`` method.
        - Values above a certain x value may be dropped by calling the ``drop_above()`` method.
        """
        self.x = []
        self.y = []
        self.decpl = decpl
        self.empty = empty
        self.filler = filler

        if empty is False and any([val is None for val in [start, end]]):
            raise ValueError(f'A start and end value must be specified for a filled '
                             f'{self.__class__.__name__} instance. ')

        # set start and end values for the spectrum
        if start is not None:
            self._start = start
        if end is not None:
            self._end = end

        if self.empty is False:
            self.x, self.y = full_spectrum_list(  # m/z and intensity lists
                self.start,
                self.end,
                decpl=decpl,
                filler=filler,
            )

        if specin is not None:
            self.add_spectrum(specin[0], specin[1])

    def __str__(self):
        return f'Full spectrum from {self.start} to {self.end} keeping {self.decpl} decimal places'

    def __repr__(self):
        return f'{self.__class__.__name__}({self.start}, {self.end}, {self.decpl})'

    def __getinitargs__(self):
        return (
            self.decpl,
            self.start,
            self.end,
            self.empty,
            self.filler,
            [self.x, self.y]
        )

    def __reduce__(self):
        return (
            self.__class__,
            self.__getinitargs__()
        )

    def __copy__(self):
        return Spectrum(
            *self.__getinitargs__()
        )

    def __deepcopy__(self, memodict={}):
        return self.__copy__()

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
        kwargs = {
            'empty': self.empty,
            'filler': self.filler,
        }
        if isinstance(x, self.__class__) is True:  # if it is another Spectrum instance
            if x.decpl != self.decpl:
                raise ValueError(
                    'The decimal places of the two spectra to be added are not equal. Addition is not supported')
            newstart = min(self.start, x.start)  # find new start m/z
            newend = max(self.end, x.end)  # find new end m/z
            tempspec = Spectrum(  # temporary instance with self as specin
                self.decpl,
                start=newstart,
                end=newend,
                specin=[self.x, self.y],
                **kwargs,
            )
            tempspec.add_spectrum(x.x, x.y)  # add incoming spectrum
            return tempspec
        elif type(x) is int:  # add this integer to every m/z
            tempspec = Spectrum(
                self.decpl,
                start=self.start,
                end=self.end,
                **kwargs
            )
            y = np.asarray(self.y)
            y += x
            tempspec.y = y.tolist()
            return tempspec
        elif len(x) == 2 and len(x[0]) == len(x[1]):  # if it is a list of paired lists (another spectrum)
            tempspec = Spectrum(
                self.decpl,
                start=self.start,
                end=self.end,
                **kwargs,
            )
            tempspec.y = list(self.y)
            tempspec.add_spectrum(x[0], x[1])
            return tempspec
        else:
            return 'Addition of %s to the Spectrum class is unsupported' % str(x)

    def __sub__(self, x):
        kwargs = {
            'empty': self.empty,
            'filler': self.filler,
        }
        if isinstance(x, self.__class__) is True:  # if it is another Spectrum instance
            if x.decpl != self.decpl:
                raise ValueError(
                    'The decimal places of the two spectra to be added are not equal. Subtraction is not supported')
            newstart = min(self.start, x.start)  # find new start m/z
            newend = max(self.end, x.end)  # find new end m/z
            tempspec = Spectrum(
                self.decpl,
                start=newstart,
                end=newend,
                specin=[self.x, self.y],
                empty=self.empty,
                filler=self.filler
            )  # temporary instance
            tempspec.add_spectrum(x.x, x.y, True)  # subtract incoming spectrum
            return tempspec
        elif type(x) is int:  # add this integer to every m/z
            tempspec = Spectrum(
                self.decpl,
                start=self.start,
                end=self.end,
                **kwargs
            )
            y = np.asarray(self.y)
            y -= x
            tempspec.y = y.tolist()
            return tempspec
        elif len(x) == 2 and len(x[0]) == len(x[1]):  # if it is a list of paired lists (another spectrum)
            tempspec = Spectrum(
                self.decpl,
                start=self.start,
                end=self.end,
                **kwargs
            )
            tempspec.y = list(self.y)
            tempspec.add_spectrum(x[0], x[1], True)  # subtract the incoming spectrum
            return tempspec
        else:
            return 'Subtraction of %s from the Spectrum class is unsupported' % str(x)

    def __mul__(self, x):
        raise AttributeError('Multiplication of the Spectrum class is unsupported')

    def __truediv__(self, x):
        raise AttributeError('Division of the Spectrum class is unsupported')

    def __pow__(self, x):
        raise AttributeError('Raising a Spectrum instance to a power is unsupported.\nAlso... really?!')

    @property
    def start(self):
        """The start value for the spectrum object"""
        return self._start

    @start.setter
    def start(self, value):
        if value is None:
            value = -np.inf
        value = round(value, self.decpl)
        if value > self._start:  # if trimming is required
            index = self.index(value)  # find index
            del self.x[:index]  # trim spectra
            del self.y[:index]
        self._start = value

    @start.deleter
    def start(self):
        self._start = -np.inf

    @property
    def end(self):
        return self._end

    @end.setter
    def end(self, value):
        if value is None:
            value = np.inf
        value = round(value, self.decpl)
        if value < self._end:
            index = self.index(value)  # find index
            self.x = self.x[:index]  # trim lists
            self.y = self.y[:index]
        self._end = value

    @end.deleter
    def end(self):
        self._end = np.inf

    @property
    def charge(self):
        """Charge for the spectrum (in mass spectrometry, the x values are mass over charge)"""
        return self._charge

    @charge.setter
    def charge(self, charge):
        if charge == self._charge:  # if already set, ignore
            return
        try:  # if numpy array, cheat
            self.x /= charge
        except TypeError:  # otherwise iterate over list
            for ind, val in enumerate(self.x):
                self.x[ind] = val / (charge / self._charge)
        # set new bounds
        self.start /= charge
        self.end /= charge
        self._charge = charge

    @charge.deleter
    def charge(self):
        setattr(self, 'charge', 1)

    def add_element(self, masses, abunds):
        """
        Adds the masses and abundances of an element to the current spectrum object.
        This is more efficient than creating a new spectrum object every time an
        element is added.

        :param list masses: List of masses (*x* values).
        :param list abunds: abundances (*y* values, paired with ``masses``)

        For example, to add a single atom of carbon to the ``Spectrum`` object

            >>> Spectrum.add_element(
                [12.0, 13.0033548378],
                [0.9893, 0.0107]
            )

        **Note**

        This function will encounter an error if there are None values in the y list.
        If you intend to use this function, set the *filler* keyword argument to be
        some value that is not None (i.e. ``0.``).

        """
        if len(masses) != len(abunds):
            raise ValueError(
                f'The dimensions of the supplied lists are not equal ({len(masses)} != {len(abunds)})')
        if self.filler is None and self.count_none() > 0:
            raise ValueError('add_element cannot operate on a y list populated with None values')

        # create matricies of new x and y values
        newx = np.asarray(self.x) + np.asarray(
            [[val] for val in masses]  # values must be boxed for appropriate combination
        )
        newy = np.asarray(self.y) * np.asarray(
            [[val] for val in abunds]
        )

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
        # self.addspectrum(oldx,oldy,True) # subtract old spectrum
        # self.dropbelow(min(masses) + min(self.x)) # drop values below new start point

        tempspec = Spectrum(
            self.decpl,
            start=min(masses) + self.start - 10 ** -self.decpl,
            end=max(masses) + self.end + 10 ** -self.decpl,
            empty=self.empty,
            filler=self.filler,
        )
        for x, y in zip(newx, newy):
            tempspec.add_spectrum(x, y)

        # for i in range(newx.shape[0]):
        #     tempspec.addspectrum(newx[i], newy[i])
        self.x = tempspec.x  # redefine the x and y lists
        self.y = tempspec.y
        self._start = min(masses) + self.start - 10 ** -self.decpl
        self._end = max(masses) + self.end + 10 ** -self.decpl

    def add_value(self, xval, yval, subtract=False):
        """
        Adds an intensity value to the x value specified.

        :param float xval: The *x* value. This value will be rounded to the decimal place specified on calling this class.
        :param float yval:  The *y* value to add to the *y* list.
        :param bool subtract: Make this ``True`` if you wish to subtract the *y* value from the current *y* value at
            the specified x.

        **Examples**

            >>> spec = Spectrum(3)
            >>> spec.trim()
            [[], []]

            >>> spec.add_value(673.9082342357,100)
            >>> spec.trim()
            [[673.908], [100]]

            >>> spec.add_value(1523.25375621,200)
            >>> spec.add_value(50.89123,300)
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
                if self.empty is False:  # if x list filled
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

    def add_spectrum(self, x, y, subtract=False):
        """
        Adds an entire x and y list to the spectrum object.
        This avoids having to call ``Spectrum.addvalue()`` in loop form
        for a list of values.

        :param list x: List of x values. These may be unsorted, but are assumed to be paired with the supplied *y* list.
        :param list y: List of y values, paired with *x*.
        :param bool subtract: Whether or not to subtract the y intensities from the current Spectrum object.
        """
        if len(x) != len(y):
            raise ValueError('The add_spectrum() method only supports two lists of the same dimension')
        for ind, mz in enumerate(x):
            if y[ind] != self.filler:  # drops filler values at this point
                self.add_value(mz, y[ind], subtract)

    def check_none(self):
        """counts the number of not-None values in the current *y* list (for debugging)"""
        return len(self.y) - self.count_none()

    def count_none(self):
        """counts the number of None values in the current *y* list (for debugging)"""
        return self.y.count(None)

    def consolidate(self, threshold, within, method='abs'):
        """
        A method of reducing the number of values in the spectrum object by consolidating y values below the specified
        threshold with nearby values. The method of combination is a weighted average. The intensities of adjacent
        values are combined until the threshold is passed or until no adjacent values within the specified x delta
        can be found.

        :param float threshold: The threshold value, below which the value will be consolidated into adjacent peaks.
        :param float within: The x delta to look within when consolidating peaks.
        :param 'abs' or 'rel' method: Whether to use an absolute or relative threshold.
        :return:
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
                        if self.x[index + 1] - self.x[index] <= within:  # if the greater than is closer
                            i = index + 1
                    else:
                        i = index + 1
            return i

        for ind in range(len(self.y)):
            if self.y[ind] < threshold and self.y[ind] != 0.:
                cur = ind
                closest = adjacent(cur)  # looks for adjacent peaks within the delta
                while self.y[cur] < threshold and closest is not None:
                    self.add_value(  # subtract the current value
                        self.x[cur],
                        self.y[cur],
                        True
                    )
                    self.add_value(  # subtract the adjacent value
                        self.x[closest],
                        self.y[closest],
                        True
                    )
                    wx, wy = weighted_average(  # weighted average of removed values
                        [self.x[cur], self.x[closest]],
                        [self.x[cur], self.y[closest]]
                    )
                    self.add_value(wx, wy)  # add the weighted average to the spectrum
                    cur = self.index(wx)  # set current index to that of the new value
                    # cur = closest # set current index to the one being tested
                    closest = adjacent(cur)
        self.threshold(threshold, method)  # drop any peaks that could not be combined

    def cp(self):
        """returns a list (clone) of the spectrum"""
        return [list(self.x), list(self.y)]

    def fill_with_zeros(self, value=0.):
        """
        Replaces any ``None`` values in the *y* list with the specified value.

        :param value: value to replace ``None`` with.
        :return: new y list
        :rtype: list
        """
        for ind, inten in enumerate(self.y):
            if inten is None:
                self.y[ind] = value
        return self.y

    def index(self, xval):
        """
        Locates the index of the specified x value in the object's x list.

        :param float xval: The x value to locate in the list.
        :return: The integer index for the x value in the x list.
        :rtype: int

        **Notes**

        If the *empty* keyword argument is True, the index will be located using the bisect module. If the x value is
        not in the current x list, the appropriate insertion index is returned. If the *empty* keyword argument is
        False, the index will be calculated based on the ``start`` and ``decpl`` values of the object (this is more
        computationally efficient than bisection).
        """
        if xval > self.end or xval < self.start:
            raise ValueError(
                f'The x value {xval} is outside of the x-list range of this {self.__class__.__name__} instance '
                f'({self.start}, {self.end}).'
            )
        if self.empty is True:  # if spectrum is unfilled, searching is required
            return bl(self.x, round(xval, self.decpl))
        else:  # otherwise, calculation of the index is more efficient
            return int(
                round(  # round after multiplication
                    (xval - self.start) * (10 ** self.decpl)  # calculate index location
                )
            )

    def nearest_x_index(self, xval):
        """
        Finds the index of the closest x value to the one provided. This method differs from `index()` in that this
        finds the closest value and index finds the insertion point to maintain an ordered list.

        :param xval: x value to find
        :return: index of nearest value
        """
        if xval < self.start:
            return 0
        if xval > self.end:
            return len(self.x) - 1
        if self.empty is False:
            return self.index(xval)
        index = self.index(xval)
        if index == len(self.x):
            return len(self.x) - 1
        potentials = [index]
        if index + 1 != len(self.x):
            potentials.append(index + 1)
        if index != 0:
            potentials.append(index - 1)
        return min(
            potentials,
            key=lambda x: abs(xval - self.x[x])
        )

    def keep_top_n(self, n=5000):
        """
        Keeps the top n peaks and sets the intensity of those below that value to be zero.

        :param int n: The number of values to keep in the list.

        **Notes**

        If there is more than one y value equal to the nth lowest value, then all of those values will be retained.

        """
        if n > len(self.x):  # do nothing if number is longer than the number of values in the Spectrum object
            return
        self.threshold(  # use the threhold method
            sorted(  # sort the y list in reverse order
                self.y,
                reverse=True,
            )[n],  # use the value at the nth index as the threshold value
            'abs',  # apply the absolute threshold
        )

    def max(self):
        """
        Locates the maximum intensity in the y list and returns the x and y values of that point.

        :return: Returns the x and y values of the maximum y value.
        :rtype: tuple of float
        """
        locs = np.where(  # locate index of maximum values
            np.asarray(self.y) == max(self.y)
        )
        if len(locs[0]) > 1:  # if there is more than one value equal to the maximum
            out = []
            for i in locs[0]:
                out.append([self.x[i], self.y[i]])
            return out
        return self.x[locs[0][0]], self.y[locs[0][0]]

    def normalize(self, new_top=100.):
        """
        Normalizes the y values to the specified value.

        :param float new_top: The new value for the maximum y value.
        """
        # old numpy way (probably less efficient)
        # m = max(self.y)  # current maximum y value
        # self.y = np.asarray(self.y)
        # self.y /= m
        # self.y *= new_top
        # self.y = self.y.tolist()
        scalar = new_top / max(self.y)  # calculate the appropriate scalar
        for ind, inten in enumerate(self.y):
            if inten is not None:
                self.y[ind] *= scalar

    def shift_x(self, value):
        """
        Offsets all *x* values by the specified value.

        :param float value: The amount to offset the x values by.
        """
        for ind, val in enumerate(self.x):
            self.x[ind] += value
        if self.start != -np.inf:
            self.start += value
        if self.end != np.inf:
            self.end += value

    def sum(self):
        """
        Calculates and returns the sum of all y values in the object.

        :return: The sum of all y values in the y list.
        :rtype: float
        """
        return sum(
            [y for y in self.y if y is not None]
        )

    def reset_y(self):
        """
        Resets the y values in the Spectrum object. This allows reuse of the same Spectrum object without regenerating.
        """
        self.y = [self.filler for val in self.y]

    def threshold(self, thresh, method='abs'):
        """
        Removes all y values below the specified threshold value.

        :param float thresh: The threshold y value to drop below.
        :param 'abs' or 'rel' method: Whether the specifed *thresh* value is absolute or relative to the maximum y
            value.
        """
        if method == 'rel':  # if relative, calculate relative threshold
            thresh *= max(self.y)
        if self.empty is True:  # removes values from the list
            x = []
            y = []
            for ind, inten in enumerate(self.y):
                if inten >= thresh:
                    x.append(self.x[ind])
                    y.append(inten)
            self.x = x
            self.y = y
        else:
            for ind, inten in enumerate(self.y):
                if inten < thresh:
                    self.y[ind] = self.filler  # sets the value to the filler value

    def trim(self, zeros=False, xbounds=None):
        """
        Trims x and y pairs that have None intensity and returns the trimmed list.
        This is the most efficient way of converting a Spectrum object to an x and y list.

        :param bool zeros: Specifies whether there should be zeros at the start and end values. This can be used to
            generate continuum spectra across the range [start,end]. If there are non-zero intensity values at the
            start or end point, they will not be affected.
        :param list xbounds: This can specify a subsection of the x and y spectra to trim to. None will return the entire
            contents of the Spectrum object, and specifying ``[x1,x2]]`` will return the x and y lists between
            *x1* and *x2*.
        :return: trimmed spectrum in the form ``[[x values], [y values]]``
        :rtype: list of lists
        """
        # retrieve boundaries
        if xbounds is None:
            xbounds = [self.start, self.end]
        elif xbounds[0] is None:
            xbounds[0] = self.start
        elif xbounds[1] is None:
            xbounds[1] = self.end

        xout = []
        yout = []

        for ind in range(self.index(xbounds[0]), self.index(xbounds[1])):  # iterate over slice
            if self.y[ind] is not self.filler:
                xout.append(round(self.x[ind], self.decpl))  # rounded to avoid array floating point weirdness
                yout.append(self.y[ind])

        if zeros is True:  # if zeros was specified, check for and insert values as necessary
            if len(xout) == 0:  # if there is no intensity in the spectrum
                xout = [float(self.start), float(self.end)]
                yout = [0., 0.]
            if xout[0] != self.start:
                xout.insert(0, self.start)
                yout.insert(0, 0.)
            if xout[-1] != self.end:
                xout.append(self.end)
                yout.append(0.)
        return [xout, yout]


def check_indexing(n=1000, dec=3):
    """
    Validates the indexing functionality of the Spectrum class

    :param n: number of iterations
    :param dec: decimal place for the Spectrum object
    :return: number of mismatches, details
    """
    spec = Spectrum(3)
    mismatch = 0
    mml = []
    for i in range(n):
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
    pass
    # spec = Spectrum(3)
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
