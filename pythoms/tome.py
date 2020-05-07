"""
Tome v02 A compilation of all Lars' python scripts as callable functions

IGNORE:
functions:
    autoresolution (estimates the resolution of a spectrum)
    bindata (bins a list of values)
    binnspectra (bins n mass spectra together into a single mass spectrum)
    bincidspectra (bins mass spectra together based on their collision voltage)
    filepresent (checks for a file or directory in the current working directory)
    find_all (finds all locations of files of a given name in a provided directory)
    linmag (generates a list of values which is linear in magnification)
    linramp (generates a list of values which is linear from start to finish)
    locateinlist (locates a value or the closest value to it in a sorted list)
    lyround (rounds a number given a particular base number)
    mag (calculates and returns the magnification of a given y value relative to the start)
    normalize (normalizes a list to a given value)
    plotms (plots a mass spectrum)
    sigmafwhm (cacluates sigma and fwhm from a resolution and a mass)
    strtolist (converts a string to a list)  
    version_input (uses the appropriate user input function depending on the python version)      

changelog:
    created mzML class and moved many functions to work within that class (removed several functions from Tome)
    added strtolist
    moved classes to separate files
    fullspeclist has been moved to _Spectrum class (there were issues with mutation of the original)
    calcindex has also been moved to _Spectrum class (it is used solely in that class)
    moved colours to _Colour class
    removed automz (now handled in the Molecule class)
    created bincidspectra to bin spectra with the same cid together
    removed loadwb, openpyxlcheck, pullparams (now included in XLSX class)
    generalized filepresent
    removed pwconvert (now included in mzML class)
    completely rewrote resolution
    rewrote resolution again to check multiple portions of the spectrum
    significant change to plotms
    moved alpha to XLSX class
    ---v02---
IGNORE
"""
import os
import sys
import scipy as sci
import numpy as np
from .spectrum import Spectrum
from bisect import bisect_left, bisect_right
from .colour import Colour
from .molecule import IPMolecule
import pylab as pl

# ----------------------------------------------------------
# -------------------FUNCTION DEFINITIONS-------------------
# ----------------------------------------------------------


def resolution(x, y, index=None, threshold=5):
    """
    Finds the resolution and full width at half max of a spectrum

    :param x: list of mz values
    :param y: corresponding list of intensity values
    :param index: index of maximum intensity (optional; used if the resolution of a specific peak is desired)
    :param threshold: signal to noise threshold required to output a resolution
    :return: resolution
    """
    y = sci.asarray(y)  # convert to array for efficiency
    if index is None:  # find index and value of maximum
        maxy = max(y)
        index = sci.where(y == maxy)[0][0]
    else:
        maxy = y[index]
    # if intensity to average is below this threshold (rough estimate of signal/noise)
    if maxy / (sum(y) / len(y)) < threshold:
        return None
    halfmax = maxy / 2
    indleft = int(index) - 1  # generate index counters for left and right walking
    indright = int(index) + 1
    while y[indleft] > halfmax:  # while intensity is still above halfmax
        indleft -= 1
    while y[indright] > halfmax:
        indright += 1
    return x[index] / (x[indright] - x[indleft])  # return resolution (mz over full width at half max)


def autoresolution(x, y, n=10, v=True):
    """
    Attempts to determine the resolution of a provided spectrum by finding n pseudo-random
    samples, then finding a peak in each of those samples to determine the resolution.


    **Parameters**

    x: *list*
        List of x values (1D list)

    y: *list*
        List of y values (1D list, must be the same length as *x*)

    n: *int*, optional
        Number of sections to check in the supplied spectrum

    v: *Bool*, optional
        Verbose toggle

    
    **Returns**
    
    resolution: *float*
        The average resolution value determined by the function

    """
    if len(x) == 0 or len(y) == 0:
        raise ValueError('the function has been handed an empty list')

    if v is True:
        sys.stdout.write('\rEstimating resolution of the spectrum')

    # find some peaks in the spectrum
    split = int(len(y) / n)
    start = 0
    end = start + split
    splity = []
    for i in range(n):
        splity.append(sci.asarray(y[start:end]))
        start += split
        end += split
    inds = []
    for ind, section in enumerate(splity):
        maxy = max(section)
        if maxy == max(section[1:-1]):  # if max is not at the edge of the spectrum
            inds.append(sci.where(section == maxy)[0][0] + split * ind)

    res = []
    for ind in inds:  # for each of those peaks
        res.append(resolution(x, y, ind))
    res = [y for y in res if y is not None]  # removes None values (below S/N)
    res = sum(res) / len(res)  # calculate average
    if v is True:
        sys.stdout.write(': %.1f\n' % res)
    return res  # return average


def bindata(n, lst, v=1):
    """
    Bins a list of values into bins of size *n*.

    **Parameters**
    
    n: *int*
        Number of values to bin together. e.g. ``n = 4`` would bin the first four values into a single value, then the next 4, etc.

    lst: *list*
        List of values to bin.

    v: *int* or *float*, optional
        Bin scalar. The calculated bin values will be divided by this value. e.g. if ``n = v`` the output values will be an average of each bin.

    **Returns**

    binned list: *list*
        A list of binned values.


    **Notes**

    - If the list is not divisible by `n`, the final bin will not be included in the output list. (The last values will be discarded.)

    """
    out = []
    delta = 0
    ttemp = 0
    for ind, val in enumerate(lst):
        delta += 1
        ttemp += val  # add current value to growing sum
        if delta == n:  # critical number is reached
            out.append(ttemp / float(v))  # append sum to list
            delta = 0  # reset critical count and sum
            ttemp = 0
    return out


def binnspectra(lst, n, dec=3, start=50., end=2000.):
    """
    Sums n spectra together.

    **Parameters**

    lst: *list*
        A list of paired lists of the form ``[ [[x1,x2,...,xn],[y1,y2,...,yn]] , [[],[]] ,...]``
        where each index of the parent list is one paired spectrum of x and y values.
        The x values of one index do not have to be the same. The spectra will be combined based on the x value rounded to the nearest 10^-`dec`.

    n: *int*
        The number of adjacent spectra to bin together. e.g. ``n = 4`` would bin the first four spectra into a single spectrum, then the next 4, etc.
    
    dec: *int*
        The decimal place to track the x values to. e.g. ``dec = 3`` would track x values to the nearest 0.001 (10^-3)

    start: *float*, optional
        The minimum x value to track in the summed spectra.

    end: *float*, optional
        The maximum x value to track in the summed spectra.

    **Returns**

    binned spectrum list: *list*
        A list of paired lists (similar to *lst*) where each index is a binned spectrum.
        If there is only one item in the binned spectra, this returns a single paired list
        of the form ``[[x values],[y values]]``.
    """
    out = []
    delta = 0
    spec = Spectrum(
        dec,
        start=start - 1,
        end=end + 1,
    )
    for ind, (x, y) in enumerate(lst):  # for each timepoint
        delta += 1
        sys.stdout.write('\rBinning spectrum #%i/%i  %.1f%%' % (ind + 1, len(lst), float(ind) / float(len(lst)) * 100.))
        spec.add_spectrum(x, y)  # add spectrum
        if delta == n:  # critical number is reached
            out.append(spec.trim(zeros=True))  # append list
            spec.reset_y()  # reset y list in object
            delta = 0  # reset critical sum
    sys.stdout.write(' DONE\n')
    if len(out) == 1:  # if there is only one item
        return out[0]
    return out


def bincidspectra(speclist, celist, dec=3, startmz=50., endmz=2000., threshold=0, fillzeros=False):
    """
    Bins mass spectra together based on the collision voltage of associated with each spectrum.

    **Parameters**

    speclist: *list*
        A list of lists of the form ``[ [[x1,x2,...,xn],[y1,y2,...,yn]] , [[],[]] ,...]``
        where each index of the parent list is one paired spectrum of x and y values.
        The x values of one index do not have to be the same. The spectra will be combined based on the x value rounded to the nearest 10^-`dec`.

    celist: *list*
        A list of collision energy values, where each index corresponds to the spectrum at that index of *speclist*. This list must be the same length as *speclist*.

    dec: *int*
        The decimal place to track the x values to. e.g. ``dec = 3`` would track x values to the nearest 0.001 (10^-3)
    
    startmz: *float*, optional
        The minimum mass to charge value to track in the summed spectra.

    end: *float*, optional
        The maximum mass to charge value to track in the summed spectra.

    threshold: *float*, optional
        The minimum y value intensity to track.

    fillzeros: *bool*, optional
        Whether to fill the resulting spectra with 0. for every value of x that does not have intensity.

    **Returns**

    specout: *list*
        A list of paired lists (similar to *speclst*) where each index is a binned spectrum.

    cv: *list*
        A sorted list of collision voltages with each index corresponding to that index in *specout*.

    """
    binned = {}

    for ind, ce in enumerate(celist):
        sys.stdout.write('\rBinning spectrum by CID value #%i/%i  %.1f%%' % (
        ind + 1, len(celist), float(ind + 1) / float(len(celist)) * 100.))
        if ce not in binned:  # generate key and spectrum object if not present
            binned[ce] = Spectrum(dec, start=startmz, end=endmz)
        else:  # otherwise add spectrum
            binned[ce].add_spectrum(speclist[ind][0], speclist[ind][1])

    if threshold > 0 or fillzeros is True:  # if manipulation is called for
        for vol in binned:  # for each voltage
            sys.stdout.write('\rZero filling spectrum for %s eV' % str(vol))
            if threshold > 0:
                binned[vol].threshold(threshold)  # apply threshold
            if fillzeros is True:
                binned[vol].fill_with_zeros()  # fill with zeros
        sys.stdout.write(' DONE\n')

    cv = []  # list for collision voltages
    specout = []  # list for spectra
    for vol, spec in sorted(binned.items()):
        sys.stdout.write('\rTrimming spectrum for %s eV' % str(vol))
        cv.append(vol)  # append voltage to list
        specout.append(spec.trim())  # append trimmed spectrum to list
    sys.stdout.write(' DONE\n')
    sys.stdout.flush()
    return specout, cv


def find_all(fname, path):
    """
    Finds all files matching a specified name within the directory specified.

    **Parameters**

    fname: *string*
        The name of the file to be located

    path: *string*
        The absolute directory path to search.

    
    **Returns**

    list of locations: *list*
        A list of all possible paths matching the filename in the specified directory.

    """
    locations = []
    for root, dirs, files in os.walk(path):
        if fname in files:
            locations.append(os.path.join(root, fname))
    return locations


def linmag(vali, magstart, magend, dur):
    """
    Generates a ramp of values that is linear in magnification.

    **Parameters**

    vali: *float*
        The initial y value at the start of the ramp.

    magstart: *float*
        The magnification at the start of the ramp.

    magend: *float*
        The magnification at the end of the ramp.

    dur: *int*
        The desired number of steps to get from *magstart* to *magend*.


    **Returns**

    list of magnifications: *list*
        A list of magnifications corresponding to the ramp.

    """
    out = []
    for i in range(dur):
        out.append(float(vali) / ((magend - magstart) / dur * i + magstart))
    return out


def linramp(valstart, valend, dur):
    """
    Generates a linear ramp of values.

    **Parameters**

    valstart: *float*
        The value at the start of the ramp.

    valend: *float*
        The value at the end of the ramp.

    dur: *int*
        The number of steps in the ramp.


    **Returns**

    List of ramped values: *list*

    """
    out = []
    for i in range(int(dur)):
        out.append(((float(valend - valstart)) / (float(dur))) * i + valstart)
    return out


def locate_in_list(lst, value, bias='closest', within=0.1):
    """
    Finds index in a sorted list of the value closest to a given value

    If two numbers are equally close, return the smallest number.
    roughly based on http://stackoverflow.com/questions/12141150/from-list-of-integers-get-number-closest-to-a-given-value

    :param lst: list of values to search
    :param value: value number to find
    :param bias: 'lesser' will return the position of the value just less than the provided value.
        'greater' will return the position of the value just greater than the provided value.
        'closest' will return the index of the nearest value to the one provided
    :param within: If the bias is closest, the position will only be returned if the position is this value away from
        the actual value
    :return: index of the position
    :rtype: int
    """
    pos = bisect_left(lst, value)
    if pos == 0:  # if at start of list
        return pos
    elif pos == len(lst):  # if insertion is beyond index range
        return pos - 1
    if lst[pos] == value:  # if an exact match is found
        return pos
    if bias == 'greater':  # return value greater than the value (bisect_left has an inherent bias to the right)
        return pos
    if bias == 'lesser':  # return value lesser than the provided
        return pos - 1
    if bias == 'closest':  # check differences between index and index-1 and actual value, return closest
        adjval = abs(lst[pos - 1] - value)
        curval = abs(lst[pos] - value)
        if adjval > within and curval > within:  # if the value is outside of the lookwithin bounds
            return None
        if adjval < curval:  # if the lesser value is closer
            return pos - 1
        if adjval == curval:  # if values are equidistant
            return pos - 1
        else:
            return pos


def lyround(x, basen):
    """
    Rounds the specified number using a specific base

    **Parameters**

    x: *float*
        The value to be rounded

    basen: *int*
        The number base to use for rounding


    **Returns**

    value: *float*
        The rounded value.

    **Notes**

    This function is based on http://stackoverflow.com/questions/2272149/round-to-5-or-other-number-in-python
    """
    base = basen ** (int(len(str(int(x)))) - 1)
    return int(base * round(float(x) / base))


def mag(initial, current):
    """
    Calculates the magnification of a specified value

    **Parameters**

    intial: *float*
        initial value (magnificiation of 1)

    current: *float*
        current value


    **Returns**

    magnification: *float*
        the magnification of the current value
    """
    return float(initial) / float(current)


def normalize(lst, maxval=1.):
    """
    Normalizes a list of values with a specified value.

    **Parameters**

    lst: *list*
        List of values to be normalized

    maxval: *float*, optional
        The maximum value that the list will have after normalization.


    **Returns**

    normalized list: *list*
        A list of values normalized to the specified value.

    """
    listmax = max(lst)
    for ind, val in enumerate(lst):
        lst[ind] = float(val) / float(listmax) * maxval
    return lst


def localmax(x: list, y: list, xval: float, lookwithin: float = 1.):
    """
    Finds the local maximum within +/- lookwithin of the xval

    :param x: x list
    :param y: y list
    :param xval:
    :param lookwithin:
    :return: maximum y value
    :rtype: float
    """
    l = bisect_left(x, xval - lookwithin)
    r = bisect_right(x, xval + lookwithin)
    return max(y[l:r])


def trimspectrum(x: list, y: list, left: float, right: float, outside: bool = False):
    """
    Trims a spectrum to the left and right bounds specified

    :param x: x value list
    :param y: y value list
    :param left: left trim value
    :param right: right trim value
    :param outside: Whether to include the next point outside of the trimmed spectrum. This provides continuity if the
        spectrum is to be used for image generation.
    :return: new spectrum
    :rtype: tuple of list
    """
    # find indicies
    l = locate_in_list(x, left, 'greater')
    r = locate_in_list(x, right, 'lesser')
    if outside is True:
        l -= 1
        r += 1
    return x[l:r + 1], y[l:r + 1]  # trim spectrum


def estimated_exact_mass(
        x: list,
        y: list,
        em: float,
        simmin: float,
        simmax: float,
        lookwithin: float = 1,
):
    """
    Estimates the exact mass of a peak in a spectrum within a provided simulated set of bounds

    :param x: x value list
    :param y: y value list
    :param em: estimated exact mass for the species
    :param simmin: minimum bound for the simulated isotope pattern
    :param simmax: maximum bounds for the simulated isotope pattern
    :param lookwithin: +/- bounds for the search
    :return: estimated exact mass
    :rtype: float
    """
    # narrow range to that of the isotope pattern
    l = bisect_left(x, simmin - lookwithin)
    r = bisect_right(x, simmax + lookwithin)
    locmax = max(y[l:r])  # find local max in that range
    for ind, val in enumerate(y):
        if val == locmax:  # if the y-value equals the local max
            if l <= ind <= r:  # and if the index is in the range (avoids false locations)
                return x[ind]
    difleft = abs(em - simmin)
    difright = abs(em - simmax)
    return '>%.1f' % max(difleft, difright)  # if no match is found, return maximum difference


# TODO change simdict to be nonmutable
def plot_mass_spectrum(
        realspec,
        simdict={},
        mz='auto',  # m/z bounds for the output spectrum
        outname='spectrum',  # name for the output file
        output='save',  # 'save' or 'show' the figure
        simtype='bar',  # simulation overlay type ('bar' or 'gaussian')
        spectype='continuum',  # spectrum type ('continuum' or 'centroid')
        maxy='max',  # max or value
        norm=True,  # True or False
        simnorm='spec',  # top, spec, or value
        xlabel=True,  # show x label
        ylabel=True,  # show y label
        xvalues=True,  # show x values
        yvalues=True,  # show y values
        showx=True,  # show x axis
        showy=True,  # how y axis
        offsetx=True,  # offset x axis (shows low intensity species better)
        fs=16,  # font size
        lw=1.5,  # line width for the plotted spectrum
        axwidth=1.5,  # axis width
        simlabels=False,  # show labels isotope for patterns
        bw='auto',  # bar width for isotope patterns (auto does 2*fwhm)
        specfont='Arial',  # the font for text in the plot
        size=[7.87, 4.87],  # size in inches for the figure
        dpiout=300,  # dpi for the output figure
        exten='png',  # extension for the output figure
        resolution=None,  # resolution to use for simulations (if not specified, automatically calculates)
        res_label=False,  # output the resolution of the spectrum
        delta=False,  # output the mass delta between the spectrum and the isotope patterns
        stats=False,  # output the goodness of match between the spectrum and the predicted isotope patterns,
        speccolour='k',  # colour for the spectrum to be plotted
        padding='auto',  # padding for the output plot
        verbose=True,  # verbose setting
        normwindow='fwhm',  # the width of the window to look for a maximal value around the expected exact mass for a peak
        annotations=None,  # annotations for the spectrum in dictionary form {'thing to print':[x,y],}
        normrel=100.,  # the maximum value for normalization
        ipmol_kwargs={},  # IPMolecule keyword arguments
        **kwargs
):
    """
    Plots and saves a publication quality mass spectrum with optional overlaid isotope patterns

    :param list realspec: A paired list of x and y values of the form ``[[x values],[y values]]``
    :param dict simdict: This can either be a molecular formula to predict the isotope pattern of (string),
        a list of formulae, or a dictionary of the form
        ``simdict = {'formula1':{'colour':<hex or name or RGB tuple>, 'alpha':float}, ...}``.
        If this is dictionary is left empty, no isotope patterns will be overlaid on the output
        spectrum.
    :param list mz: The *m/z* bounds for the output spectrum. Default: 'auto', but can be supplied
        with a tuple or list of length 2 of the form ``[x start, x end]``.
    :param str outname: Name of the file to be saved.
    :param str output: Save ('save') or show ('show') the figure.
    :param str simtype: The type for the isotope pattern simulation overlay. Options: 'bar' or 'gaussian'.
    :param str spectype: The type of spectrum being handed to the function. Options: 'continuum' or 'centroid'.
    :param float maxy: The maximum y value for the spectrum. Options: 'max' or specify a value
    :param bool norm: Normalize the spectrum. Options: bool
    :param str, float simnorm: Normalize the isotope pattern simulations to what value. Options: 'top', 'spec', or
        specify a value. Top will normalize the patterns to ``maxy``, and will only function if maxy is not 'max'.
        Spec will normalize the patterns to the maximum spectrum y value within the x bounds of the
        simulated pattern.
        Specifying a value will normalize all isotope patterns to that value.
    :param bool xlabel: Whether to show the label for the *m/z* axis.
    :param bool ylabel: Whether to show the y-axis label.
    :param bool xvalues: Whether to show the values of the x-axis.
    :param bool yvalues: Whether to show the values of the y-axis.
    :param bool showx: Whether to show the x-axis line.
    :param bool showy: Whether to show the y-axis line.
    :param bool offsetx: Whether to offset the x-axis slightly.
        Enabling this shows makes it easier to see low intensity peaks.
    :param int fs: Font size to use for labels.
    :param float lw: Line width for the plotted spectrum.
    :param float axwidth: Line width for the axes and tick marks. Default 1.5
    :param bool simlabels: Whether to show the names of the simulated isotope patterns.
        The names will be exactly as supplied in ``simdict``.
    :param float bw: The width of the bar in *m/z* for bar isotope patterns. Options: 'auto' or float.
        This only has an affect if *simtype* is 'bar'.
        Auto make the bars equal to 2 times the full width at half max of the peak they are simulating.
    :param str specfont: The font to use for text in the plot. The specified font must be accepted by matplotlib.
    :param list size: The size in inches for the output figure. This must be a list of length 2 of the form
        ``[width,height]``.
    :param int dpiout: The dots per inch for the output figure.
    :param str exten: The file extension for the output figure. Options: 'png', 'svg', or other supported by matplotlib.
    :param float resolution: Override the auto-resolution calculation with a specified instrument resolution
    :param bool res_label: Whether to output the resolution of the spectrum onto the figure.
    :param bool delta: Whether to calculate and output the mass delta between the exact mass predicted by the isotope
        pattern simulation and the location of the maximum intensity within the bounds specified by *normwindow*.
    :param bool stats: Whether to calculate and output the goodness of fit between the predicted isotope pattern and
        the supplied spectrum. This functionality is still a work in progress.
    :param speccolour: The colour for the real spectrum , # colour for the spectrum to be plotted
    :param list padding: This allows the user to specify the subplot padding of the output figure.
        Options: 'auto' or list of the form ``[left,right,bottom,top]`` scalars.
    :param bool verbose: Verbose option for the script. Options: bool.
    :param float normwindow: The *m/z* window width within with too look for a maximum intensity value.
        This will only have an effect if *delta* is ``True``.
        Options: 'fwhm' for full width at half max or float.
    :param dict annotations: Annotations for the spectrum in dictionary form: ``{'thing to print':[x,y],}``.
    :param normrel: The maximum value for normalization. This can be used to globally set the top value for normalizing
        simulated isotope patterns. This is used most often to show the lack of an isotope pattern in the shown area.
    :param ipmol_kwargs: Keyword arguments to use for IPMolecule calls. See IPMolecule for more details.
    :param kwargs: catch for unused kwargs
    """
    def checksimdict(dct):
        """
        checks the type of simdict, converting to dictionary if necessary
        also checks for alpha and colour keys and adds them if necessary (defaulting to key @ 0.5)
        """
        if type(dct) is not dict:
            if type(dct) is str:
                dct = {dct: {}}
            elif type(dct) is list or type(dct) is tuple:
                tdct = {}
                for i in dct:
                    tdct[i] = {}
                dct = tdct
        for species in dct:
            if 'colour' not in dct[species]:
                dct[species]['colour'] = 'k'
            if 'alpha' not in dct[species]:
                dct[species]['alpha'] = 0.5
        return dct

    if resolution is None:
        if spectype != 'centroid':
            resolution = autoresolution(realspec[0], realspec[1])  # calculate resolution
        else:
            resolution = 5000

    simdict = checksimdict(simdict)  # checks the simulation dictionary
    for species in simdict:  # generate Molecule object and set x and y lists
        simdict[species]['colour'] = Colour(simdict[species]['colour'])
        simdict[species]['mol'] = IPMolecule(
            species,
            resolution=resolution,
            **ipmol_kwargs,
        )
        # simdict[species]['mol'] = Molecule(species, res=res, dropmethod='threshold')
        if simtype == 'bar':
            simdict[species]['x'], simdict[species]['y'] = simdict[species]['mol'].barip
        if simtype == 'gaussian':
            simdict[species]['x'], simdict[species]['y'] = simdict[species]['mol'].gausip

    if mz == 'auto':  # automatically determine m/z range
        if verbose is True:
            sys.stdout.write('Automatically determining m/z window')
        mz = [10000000, 0]
        for species in simdict:
            simdict[species]['bounds'] = simdict[species]['mol'].bounds  # calculate bounds
            if simdict[species]['bounds'][0] < mz[0]:
                mz[0] = simdict[species]['bounds'][0] - 1
            if simdict[species]['bounds'][1] > mz[1]:
                mz[1] = simdict[species]['bounds'][1] + 1
        if mz == [10000000, 0]:
            mz = [min(realspec[0]), max(realspec[0])]
        if verbose is True:
            sys.stdout.write(': %i - %i\n' % (int(mz[0]), int(mz[1])))
            sys.stdout.flush()

    realspec[0], realspec[1] = trimspectrum(  # trim real spectrum for efficiency
        realspec[0],
        realspec[1],
        mz[0] - 1,
        mz[1] + 1
    )

    if len(realspec[0]) == 0:  # catch for empty spectrum post-trim (usually user error)
        raise ValueError(f'There are no spectral values in the specified m/z bounds ({mz[0]}-{mz[1]}). Common causes: '
                         f'no values in the loaded spectrum within the window of interest, an error in specifying the '
                         f'molecule(s) to simulate')

    if norm is True:  # normalize spectrum
        realspec[1] = normalize(realspec[1], normrel)

    for species in simdict:  # normalize simulations
        if simnorm == 'spec':  # normalize to maximum around exact mass
            if normwindow == 'fwhm':  # if default, look within the full width at half max
                window = simdict[species]['mol'].fwhm
            else:  # otherwise look within the specified value
                window = normwindow
            simdict[species]['y'] = normalize(
                simdict[species]['y'],
                localmax(
                    realspec[0],
                    realspec[1],
                    simdict[species]['mol'].estimated_exact_mass,
                    window
                )
            )
        elif simnorm == 'top':  # normalize to top of the y value
            if maxy == 'max':
                raise ValueError('Simulations con only be normalized to the top of the spectrum when the maxy setting '
                                 'is a specific value')
            simdict[species]['y'] = normalize(simdict[species]['y'], maxy)
        elif type(simnorm) is int or type(simnorm) is float:  # normalize to specified value
            simdict[species]['y'] = normalize(simdict[species]['y'], simnorm)
        if delta is True:
            if normwindow == 'fwhm':  # if default, look within the full width at half max
                window = simdict[species]['mol'].fwhm
            else:  # otherwise look within the specified value
                window = normwindow
            est = estimated_exact_mass(  # try to calculate exact mass
                realspec[0],
                realspec[1],
                simdict[species]['mol'].estimated_exact_mass,
                simmin=simdict[species]['mol'].estimated_exact_mass,
                simmax=simdict[species]['mol'].estimated_exact_mass,
                lookwithin=window,
                # min(simdict[species]['x']),
                # max(simdict[species]['x'])
            )
            if type(est) is float:
                simdict[species]['delta'] = '%.3f (%.1f ppm)' % (
                simdict[species]['mol'].estimated_exact_mass - est, simdict[species]['mol'].compare_exact_mass(est))
            else:
                simdict[species]['delta'] = est

    pl.clf()  # clear and close figure if open
    pl.close()
    fig = pl.figure(figsize=tuple(size))
    ax = fig.add_subplot(111)

    ax.spines["right"].set_visible(False)  # hide right and top spines
    ax.spines["top"].set_visible(False)

    if showx is False:
        ax.spines["bottom"].set_visible(False)  # hide bottom axis
    if showy is False:
        ax.spines["left"].set_visible(False)  # hide left axis

    for axis in ["top", "bottom", "left", "right"]:
        ax.spines[axis].set_linewidth(axwidth)

    if offsetx is True:  # offset x axis
        ax.spines["bottom"].set_position(('axes', -0.01))

    font = {'fontname': specfont, 'fontsize': fs}  # font parameters for axis/text labels
    tickfont = pl.matplotlib.font_manager.FontProperties(  # font parameters for axis ticks
        family=specfont,
        size=fs
    )

    ax.set_xlim(mz)  # set x bounds

    if maxy == 'max':  # set y bounds
        ax.set_ylim(0., max(realspec[1]))
        top = max(realspec[1])
    elif type(maxy) is int or type(maxy) is float:
        ax.set_ylim(0., maxy)
        top = maxy

    if simtype == 'bar':  # generates zeros for bottom of bars (assumes m/z spacing is equal between patterns)
        for species in simdict:
            simdict[species]['zero'] = []
            for i in simdict[species]['x']:
                simdict[species]['zero'].append(0.)
        for species in simdict:  # for each species
            for subsp in simdict:  # look at all the species
                if subsp != species:  # if it is not itself
                    ins = bisect_left(simdict[subsp]['x'], simdict[species]['x'][-1])  # look for insertion point
                    if 0 < ins < len(simdict[subsp]['x']):  # if species highest m/z is inside subsp list
                        for i in range(ins):  # add intensity of species to subsp zeros
                            # used -ins+i-1 to fix an error, with any luck this won't break it next time
                            simdict[subsp]['zero'][i] += simdict[species]['y'][-ins + i]
    # include resolution if specified (and spectrum is not centroid)
    if res_label is True and spectype != 'centroid':
        ax.text(
            mz[1],
            top * 0.95,
            f'resolution: {str(round(resolution))}',
            horizontalalignment='right',
            **font
        )

    for species in simdict:  # plot and label bars
        if simtype == 'bar':
            if bw == 'auto':
                bw = simdict[species]['mol'].fwhm * 2
            else:
                bw = bw
            ax.bar(
                simdict[species]['x'],
                simdict[species]['y'],
                bw,
                alpha=simdict[species]['alpha'],
                color=simdict[species]['colour'].mpl,
                linewidth=0,
                align='center',
                bottom=simdict[species]['zero']
            )
        elif simtype == 'gaussian':
            ax.plot(
                simdict[species]['x'],
                simdict[species]['y'],
                alpha=simdict[species]['alpha'],
                color=simdict[species]['colour'].mpl,
                linewidth=lw
            )
            ax.fill_between(
                simdict[species]['x'],
                0,
                simdict[species]['y'],
                alpha=simdict[species]['alpha'],
                color=simdict[species]['colour'].mpl,
                linewidth=0
            )
        # if any labels are to be shown
        if simlabels is True or stats is True or delta is True:
            string = ''
            bpi = simdict[species]['y'].index(max(simdict[species]['y']))  # index of base peak
            if simlabels is True:  # species name
                string += species
                if stats is True or delta is True:  # add return if SER or delta is called for
                    string += '\n'
            if stats is True:  # standard error of regression
                string += f'SER: {simdict[species]["mol"].compare(realspec)} '
            if delta is True:  # mass delta
                string += f'mass delta: {simdict[species]["delta"]}'
            ax.text(
                simdict[species]['x'][bpi],
                top * 1.01,
                string,
                color=simdict[species]['colour'].mpl,
                horizontalalignment='center',
                **font
            )

    if spectype == 'continuum':
        ax.plot(
            realspec[0],
            realspec[1],
            linewidth=lw,
            color=Colour(speccolour).mpl
        )
    elif spectype == 'centroid':
        dist = []
        for ind, val in enumerate(realspec[0]):  # find distance between all adjacent m/z values
            if ind == 0:
                continue
            dist.append(realspec[0][ind] - realspec[0][ind - 1])
        dist = sum(dist) / len(dist)  # average distance
        ax.bar(
            realspec[0],
            realspec[1],
            dist * 0.75,
            linewidth=0,
            color=Colour(speccolour).mpl,
            align='center',
            alpha=0.8
        )

    if annotations is not None:
        for label in annotations:
            ax.text(
                annotations[label][0],
                annotations[label][1],
                label,
                horizontalalignment='center',
                **font
            )

            # show or hide axis values/labels as specified
    if yvalues is False:  # y tick marks and values
        ax.tick_params(axis='y', labelleft='off', length=0)
    else:  # y value labels
        ax.tick_params(
            axis='y',
            length=axwidth * 3,
            width=axwidth,
            direction='out',
            right=False,
        )
        for label in ax.get_yticklabels():
            label.set_fontproperties(tickfont)
    if ylabel is True:  # y unit
        if top == 100:  # normalized
            ax.set_ylabel('relative intensity', **font)
        else:  # set to counts
            ax.set_ylabel('intensity (counts)', **font)

    if xvalues is False:  # x tick marks and values
        ax.tick_params(axis='x', labelbottom='off', length=0)
    else:  # x value labels
        ax.tick_params(
            axis='x',
            length=axwidth * 3,
            width=axwidth,
            direction='out',
            top=False,
        )
        for label in ax.get_xticklabels():
            label.set_fontproperties(tickfont)
    if xlabel is True:  # x unit
        ax.set_xlabel('m/z', style='italic', **font)

    pl.ticklabel_format(useOffset=False)  # don't use the stupid shorthand thing
    if padding == 'auto':
        pl.tight_layout(pad=0.5)  # adjust subplots
        if simlabels is True or stats is True or delta is True:
            pl.subplots_adjust(top=0.90)  # lower top if details are called for
    elif type(padding) is list and len(padding) == 4:
        pl.subplots_adjust(
            left=padding[0],
            right=padding[1],
            bottom=padding[2],
            top=padding[3]
        )

    if output == 'save':  # save figure
        outname = ''  # generate tag for filenaming
        for species in simdict:
            outname += ' ' + species
        outname = outname + outname + '.' + exten
        pl.savefig(
            outname,
            dpi=dpiout,
            format=exten,
            transparent=True
        )
        if verbose is True:
            sys.stdout.write('Saved figure as:\n"%s"\nin the working directory' % outname)

    elif output == 'show':  # show figure
        pl.show()


def plotuv(wavelengths, intensities, **kwargs):
    """
    Plots and saves a publication quality UV-Vis figure.

    **Parameters**

    wavelengths: *list*
        A list of wavelengths

    intensities: *list*
        A list of intensity values paired by index to *wavelengths*


    **Returns**

    return item: ``None``
        This function has no pythonic return.

    **\*\*kwargs**

    axwidth: 1.5
        Line width for the axes and tick marks. Options: float.

    colours:
        A list of colours to be used if the fuction is supplied with multiple traces.

    dpiout: 300
        The dots per inch for the output figure. Options: integer.

    exten: 'png'
        The file extension for the output figure. Options: 'png', 'svg', or other supported by matplotlib.

    fs: 16
        Font size to use for labels. Options: integer or float.

    legloc: 0
        The matplotlib legend location key.
        See http://matplotlib.org/api/legend_api.html for location codes.

    lw: 1.5
        Line width for the plotted spectrum. Options: float.

    outname: 'UV-Vis spectrum'
        Name for the output file. Options: string.

    output: 'save'
        Save ('save') or show ('show') the figure.

    padding: 'auto'
        This allows the user to specify the subplot padding of the output figure.
        Options: 'auto' or list of the form ``[left,right,bottom,top]`` scalars.

    size: [7.87,4.87]
        The size in inches for the output figure. This must be a list of length 2 of the form
        ``[width,height]``.

    specfont: 'Arial'
        The font to use for text in the plot. The specified font must be accepted by matplotlib.

    times: None
        A list of timepoints for each provided trace. These are used as labels in the legend.

    verbose: True
        Verbose option for the script. Options: bool.

    xrange: None
        The limits for the x axis. Options None or ``[x min,x max]``

    yrange: None
        The limits for the y axis. Options None or ``[y min,y max]``


    """
    settings = {  # default settings for the function
        'outname': 'UV-Vis spectrum',  # name for the output file
        'fs': 16,  # font size
        'lw': 1.5,  # line width for the plotted spectrum
        'axwidth': 1.5,  # axis width
        'size': [7.87, 4.87],  # size in inches for the figure
        'dpiout': 300,  # dpi for the output figure
        'exten': 'png',  # extension for the output figure
        'specfont': 'Arial',  # the font for text in the plot
        # colours to use for multiple traces in the same spectrum (feel free to specify your own)
        'colours': ['#a6cee3', '#1f78b4', '#b2df8a', '#33a02c', '#fb9a99', '#e31a1c', '#fdbf6f', '#ff7f00', '#cab2d6',
                    '#6a3d9a', '#ffff99', '#8dd3c7', '#ffffb3', '#bebada', '#fb8072', '#80b1d3', '#fdb462', '#b3de69',
                    '#fccde5', '#d9d9d9', '#bc80bd', '#ccebc5', ],
        'xrange': None,  # the limits for the x axis
        'yrange': None,  # the limits for the y axis
        'times': None,  # time points for each provided trace (for legend labels)
        'output': 'save',  # 'save' or 'show' the figure
        'padding': None,  # padding for the output plot
        'verbose': True,  # chatty
        'legloc': 0,  # legend location (see http://matplotlib.org/api/legend_api.html location codes)
    }
    if set(kwargs.keys()) - set(settings.keys()):  # check for invalid keyword arguments
        string = ''
        for i in set(kwargs.keys()) - set(settings.keys()):
            string += ' %s' % i
        raise KeyError('Unsupported keyword argument(s): %s' % string)

    settings.update(kwargs)  # update settings from keyword arguments

    pl.clf()  # clear and close figure if open
    pl.close()
    fig = pl.figure(figsize=tuple(settings['size']))
    ax = fig.add_subplot(111)

    ax.spines["right"].set_visible(False)  # hide right and top spines
    ax.spines["top"].set_visible(False)

    font = {'fontname': settings['specfont'], 'fontsize': settings['fs']}  # font parameters for axis/text labels
    tickfont = pl.matplotlib.font_manager.FontProperties(family=settings['specfont'],
                                                         size=settings['fs'])  # font parameters for axis ticks

    if type(intensities[0]) is float:  # if the function has only been handed a single spectrum
        intensities = [intensities]

    # determine and set limits for axes
    if settings['xrange'] is None:  # auto determine x limits
        settings['xrange'] = [min(wavelengths), max(wavelengths)]
    if settings['yrange'] is None:  # auto determine y limits
        settings['yrange'] = [0, 0]
        for spec in intensities:
            if max(spec) > settings['yrange'][1]:
                settings['yrange'][1] = max(spec)
    ax.set_xlim(settings['xrange'])  # set x bounds
    ax.set_ylim(settings['yrange'])  # set y bounds

    # apply font and tick parameters to axes
    ax.tick_params(axis='x', length=settings['axwidth'] * 3, width=settings['axwidth'], direction='out', top='off')
    for label in ax.get_xticklabels():
        label.set_fontproperties(tickfont)
    ax.tick_params(axis='y', length=settings['axwidth'] * 3, width=settings['axwidth'], direction='out', right='off')
    for label in ax.get_yticklabels():
        label.set_fontproperties(tickfont)
    for axis in ["top", "bottom", "left", "right"]:
        ax.spines[axis].set_linewidth(settings['axwidth'])

    if settings['times'] is not None:
        if len(settings['times']) != len(intensities):
            raise IndexError('The numer of times provided do not match the number of traces provided.')

    for ind, spec in enumerate(intensities):  # plot traces
        if settings['times'] is not None:
            string = 't = ' + str(round(settings['times'][ind], 1)) + 'm'
            ax.plot(wavelengths, spec, label=string, color=Colour(settings['colours'][ind]).mpl,
                    linewidth=settings['lw'])
        else:
            ax.plot(wavelengths, spec, color=Colour(settings['colours'][ind]).mpl, linewidth=settings['lw'])

    if settings['times'] is not None:
        ax.legend(loc=0, frameon=False)

    ax.set_xlabel('wavelength (nm)', **font)
    ax.set_ylabel('absorbance (a.u.)', **font)

    if settings['padding'] is None:
        pl.tight_layout(pad=0.5)  # adjust subplots
    elif type(settings['padding']) is list and len(settings['padding']) == 4:
        pl.subplots_adjust(left=settings['padding'][0], right=settings['padding'][1], bottom=settings['padding'][2],
                           top=settings['padding'][3])

    if settings['output'] == 'save':  # save figure
        outname = settings['outname'] + '.' + settings['exten']
        pl.savefig(outname, dpi=settings['dpiout'], format=settings['exten'], transparent=True)
        if settings['verbose'] is True:
            sys.stdout.write('Saved figure as:\n"%s"\nin the working directory' % outname)

    elif settings['output'] == 'show':  # show figure
        pl.show()


def sigmafwhm(res, x):
    """
    Calculates the full width at half max and standard deviation for a spectrum peak.

    **Parameters**

    res: *float*
        The resolution of the peak in question

    x: *float*
        The x value of the peak in question


    **Returns**

    fwhm: *float*
        The full width at half max of the peak.

    sigma: *float*
        The standard deviation of the peak.

    """
    fwhm = x / res
    sigma = fwhm / (2 * np.sqrt(2 * np.log(2)))  # based on the equation FWHM = 2*sqrt(2ln2)*sigma
    return fwhm, sigma


def strtolist(string):
    """
    Converts a string to a list with more flexibility than ``string.split()``
    by looking for both brackets of type ``(,),[,],{,}`` and commas.

    **Parameters**

    string: *string*
        The string to be split.


    **Returns**

    split list: *list*
        The split list


    **Examples**

    ::

        >>> strtolist('[(12.3,15,256.128)]')
            [12.3, 15, 256.128]

    """
    out = []
    temp = ''
    brackets = ['(', ')', '[', ']', '{', '}']
    for char in list(string):
        if char not in brackets and char != ',':
            temp += char
        if char == ',':
            try:
                out.append(int(temp))
            except ValueError:
                out.append(float(temp))
            temp = ''
    if len(temp) != 0:  # if there is a weird ending character
        try:
            out.append(int(temp))
        except ValueError:
            out.append(float(temp))
    return out
