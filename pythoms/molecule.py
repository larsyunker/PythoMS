"""
IGNORE:
Molecule class (previously "isotope pattern generator" and "MolecularFormula")

The output of this builder has been validated against values calculated by ChemCalc (www.chemcalc.org)
Negligable differences are attributed to different low value discarding techniques
(ChemCalc keeps the top 5000 peaks, this script drops values less than a threshold 5 orders of magnitude below the
maximum value)

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
import sys
import pickle
import os
import importlib.util
import numpy as np
from scipy import stats
from datetime import datetime
import sympy as sym
import pylab as pl
import copy
from tqdm import tqdm
from .scripttime import ScriptTime
from .spectrum import Spectrum, weighted_average
from . import mass_dictionaries  # import mass dictionaries
from itertools import combinations_with_replacement as cwr
from IsoSpecPy import IsoThreshold

# flag for reminding folk to cite people
_CITATION_REMINDER = False

# attempt to load abbreviation dictionary from current working directory
from .mass_abbreviations import abbrvs

try:
    abbrv_spec = importlib.util.spec_from_file_location(
        'user_abbrvs',
        os.path.join(
            os.getcwd(),
            'user_mass_abbreviations.py'
        )
    )
    abbrv_module = importlib.util.module_from_spec(abbrv_spec)
    abbrv_spec.loader.exec_module(abbrv_module)
    user_abbrvs = abbrv_module.user_abbrvs
    abbrvs.update(user_abbrvs)
except FileNotFoundError:  # if it can't find the file, continue with default abbreviations
    pass


"""Mass dictionary associated with the instance"""
MASS_KEY = 'crc_mass'
mass_dict = getattr(
    mass_dictionaries,
    MASS_KEY,
)


st = ScriptTime(profile=True)

# valid start and end brackets
OPENING_BRACKETS = ['(', '{', '[']  # opening brackets
CLOSING_BRACKETS = [')', '}', ']']  # closing brackets
SIGNS = ['+', '-']  # charge signs
VERBOSE = False  # toggle for verbose

# valid grouping methods
VALID_GROUP_METHODS = [
    'weighted',
    'centroid',
]
# valid isotope pattern generation methods
VALID_IPMETHODS = [
    'combinatorics',
    'multiplicative',
    'hybrid',
    'isospec',  # use isospecpy package
    # 'cuda',
]

# valid dropping methods
VALID_DROPMETHODS = [
    None,  # no dropping
    'threshold',  # drop values below threshold
    'npeaks',  # keep top n number of peaks
    # 'consolidate',  # consolidate intensities
]

# default threshold for low-intensity peak dropping
THRESHOLD = 0.01
# number of peaks to keep for low-intensity peak dropping
NPEAKS = 5000
# consolidation threshold for low-intensity peak combination
CONSOLIDATE = 3


def interpret(block: str):
    """
    Interprets an element block, breaking it into element and number of that element.

    :param block: string block describing an element
    :return: composition dictionary
    :rtype: dict
    """
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


def interpret_charge(string: str):
    """
    Interprets a charge string.

    :param string: string describing the charge (e.g. '2+')
    :return: charge, sign
    :rtype: tuple
    """
    value = ''
    sign = '+'  # default value for sign
    if type(string) is int:
        return string, sign
    for ind, val in enumerate(string):
        if val in SIGNS:  # if val sets mode
            sign = val
        else:  # number
            value += val
    if value == '':  # if no number was specified (e.g. "+")
        value = 1
    return int(value), sign


def string_to_isotope(string: str):
    """
    Attempts to interpret an undefined key as an isotope/element combination (e.g. "13C" becomes 'C', 13). Raises a
    ValueError if the string cannot be interpreted as such.

    :param string: string to interpret
    :return: element, isotope
    :rtype: (str, int)
    """
    iso = string[0]
    if iso.isdigit() is False:
        raise TypeError(f'The isotope "{string}" is not a valid format. Use isotope/element format e.g. "12C"')
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
            f'The string "{string}" could not be interpreted as an element, isotope combination, please check'
            f'your input')


unicode_subscripts = {  # subscripts values for unit representations
    0: f'\u2080',
    1: f'\u2081',
    2: f'\u2082',
    3: f'\u2083',
    4: f'\u2084',
    5: f'\u2085',
    6: f'\u2086',
    7: f'\u2087',
    8: f'\u2088',
    9: f'\u2089',
}
unicode_superscripts = {  # superscript values for unit representations
    0: f'\u2070',
    1: f'\u00b9',
    2: f'\u00b2',
    3: f'\u00b3',
    4: f'\u2074',
    5: f'\u2075',
    6: f'\u2076',
    7: f'\u2077',
    8: f'\u2078',
    9: f'\u2079',
}


def to_subscript(number):
    """
    Converts the value to subscript characters.

    :param int number: number to convert
    :return: subscript
    :rtype: str
    """
    return ''.join(
        [unicode_subscripts[int(val)] for val in str(abs(number))]
    )


def to_superscript(val):
    """
    Returns the integer value represented as a superscript string.

    :param int val: value to represent
    :return: superscript string
    :rtype: str
    """
    return ''.join(
        [unicode_superscripts[int(val)] for val in str(abs(val))]
    )


def check_in_mass_dict(comp: dict):
    """
    Checks for the presence of the dictionary keys in the mass dictionary. Raises a ValueError if the key is not found.

    :param comp: composition dictionary
    """
    for key in comp:
        if key not in mass_dict:
            ele, iso = string_to_isotope(key)
            if ele not in mass_dict:
                raise ValueError(f'The element {ele} is not defined in the mass dictionary. Please check your input.')
            elif iso not in mass_dict[ele]:
                raise ValueError(
                    f'The element "{ele}" does not have a defined isotope "{iso}" in the mass dictionary. '
                    f'Please check your input.'
                )


def element_intensity_list(element: str):
    """
    Returns the non-zero element intensity for the specified element.

    :param element: element key
    :return: mass, intensity lists
    :rtype: list
    """
    if element not in mass_dict:
        raise KeyError(f'The element {element} is not defined in the mass dictionary.')
    ele_dict = mass_dict[element]
    mass_out = []
    intensity_out = []
    for isotope in ele_dict:
        if isotope != 0 and ele_dict[isotope][1] != 0.:
            mass_out.append(ele_dict[isotope][0])
            intensity_out.append(ele_dict[isotope][1])
    return [mass_out, intensity_out]


def chew_formula(formula: str):
    """
    Iterates through provided formula, extracting blocks, interpreting the blocks,
    and returning the formula minus the blocks.

    :param formula: string formula
    :return: remaining formula, interpreted block
    :rtype: str, dict
    """
    if formula[0].isupper() is True:  # element is recognized by an uppercase letter
        block = formula[0]  # element block
        for loc in range(len(formula)):
            if loc == 0:
                continue
            if formula[loc].isupper() is True:  # if an uppercase character is encountered
                break
            elif formula[loc] in OPENING_BRACKETS:  # if a bracket is encountered
                break
            else:
                block += formula[loc]
        return formula[len(block):], interpret(block)  # return remaining formula and the interpreted block
    elif formula[0] in OPENING_BRACKETS:  # if a bracket is encountered, intialize bracket interpretation
        return bracket(formula)
    elif formula[0].isdigit() is True:  # either isotope or charge
        if any([sign in formula for sign in SIGNS]):  # if the block is a value-sign charge specification
            return '', {'charge': formula}
        for ind, val in enumerate(formula):
            if formula[ind].isalpha() is True:  # if isotope encountered, return that isotope with n=1
                return '', {formula: 1}
    elif formula[0] in SIGNS:  # charge specification
        return '', {'charge': formula}  # assign as charge for later interpretation
    else:
        raise ValueError(f'An uninterpretable formula chunck was encountered: {formula}')


def bracket(form):
    """finds the string block contained within a bracket and determines the formula within that bracket"""
    bracktype = OPENING_BRACKETS.index(form[0])  # sets bracket type (so close bracket can be identified)
    bnum = ''  # number of things indicated in the bracket
    block = ''  # element block
    nest = 1  # counter for nesting brackets
    for loc in range(len(form)):  # look for close bracket
        if loc == 0:
            continue
        elif form[loc] == OPENING_BRACKETS[bracktype]:  # if a nested bracket is encountered
            nest += 1
            block += form[loc]
        elif form[loc] == CLOSING_BRACKETS[bracktype]:  # if close bracket is encountered
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
            f'A close bracket was not encountered for the "{form[0]}" bracket in the formula segment "{form}". '
            f'Please check your input molecular formula.')

    lblock = len(block) + len(
        bnum) + 2  # length of the internal block + the length of the number + 2 for the brackets
    if bnum == '':  # if no number is specified
        bnum = 1
    else:
        bnum = int(bnum)
    outdict = {}
    while len(block) > 0:  # chew through bracket
        ftemp, tempdict = chew_formula(block)
        for key in tempdict:
            try:
                outdict[key] += tempdict[key] * bnum
            except KeyError:
                outdict[key] = tempdict[key] * bnum
        block = ftemp
    return form[lblock:], outdict  # returns remaining formula and composition of the block


def abbreviations(dct: dict):
    """
    Searches for abbreviations predefined in mass_abbreviations.py either in the pythoms package or in the current
    working directory. Any found abbreviations will be added to the current dictionary.

    :param dct: incoming dictionary
    :return: un-abbreviated dictionary
    :rtype: dict
    """
    comptemp = {}
    for key in dct:
        if key in abbrvs:  # if a common abbreviation is found in formula
            for subkey in abbrvs[key]:
                try:
                    comptemp[subkey] += abbrvs[key][subkey] * dct[key]
                except KeyError:
                    comptemp[subkey] = abbrvs[key][subkey] * dct[key]
        else:
            try:
                comptemp[key] += dct[key]
            except KeyError:
                comptemp[key] = dct[key]
    return comptemp


def composition_from_formula(formula):
    """
    Interprets a provided string as a molecular formula.
    Supports nested brackets, charges, and isotopes.

    :param formula:  A molecular formula. Charge may be specified in the formula, but care must be taken to specify
        the charge in sign-value format (e.g. '+2' if value-sign is specified, the script will attempt to interpret the
        key as an isotope).
    :return: A dictionary where each key is an element or isotope with its value
        being the number of each of the elements or isotopes. e.g. the
        molecule CH4 would have the composition ``comp = {'C':1, 'H':4}``
    :rtype: dict
    """
    comp = {}
    while len(formula) > 0:  # chew through formula
        ftemp, nomdict = chew_formula(formula)  # find the next block
        for ele in nomdict:
            try:
                comp[ele] += nomdict[ele]
            except KeyError:
                comp[ele] = nomdict[ele]
        formula = ftemp
    comp = abbreviations(comp)  # look for common abbreviations
    return comp


def standard_deviation(fwhm):
    """determines the standard deviation for a normal distribution with the full width at half max specified"""
    return fwhm / (2 * np.sqrt(2 * np.log(2)))  # based on the equation FWHM = 2*sqrt(2ln2)*sigma


def group_masses(ip, dm: float = 0.25):
    """
    Groups masses in an isotope pattern looking for differences in m/z greater than the specified delta.
    expects

    :param ip: a paired list of [[mz values],[intensity values]]
    :param dm: Delta for looking +/- within
    :return: blocks grouped by central mass
    :rtype: list
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
    this method results in substantially more error compared to the weighted_average method (~9 orders of
    magnitude for C61H51IP3Pd)
    """
    return sum(ipgroup[0]) / len(ipgroup[0]), sum(ipgroup[1]) / len(ipgroup[1])


def bar_isotope_pattern(
        rawip: list,
        delta: float = 0.5,
        method: str = 'weighted',
        verbose: bool = VERBOSE,
):
    """
    Converts a raw isotope pattern into a bar isotope pattern. This groups mass defects
    that are within a given difference from each other into a single *m/z* value and
    intensity.

    :param rawip: The raw isotope pattern with no grouping applied
    :param delta: The *m/z* difference to check around a peak when grouping it into a single *m/z* value.
        The script will look delta/2 from the peak being checked
    :param method: Method of combining (weighted or centroid). Weighted is recommended for accuracy
    :param verbose: chatty mode
    :return: bar isotope pattern in ``[[m/z values],[intensity values]]`` format.
    :rtype: list
    """
    if method not in VALID_GROUP_METHODS:
        raise ValueError(f'The grouping method {method} is invalid. Choose from {", ".join(VALID_GROUP_METHODS)}')
    if verbose is True:
        sys.stdout.write('Generating bar isotope pattern')
    if isinstance(rawip, Spectrum):  # if handed a Spectrum object, trim before continuing
        rawip = rawip.trim()
    groupedip = group_masses(rawip, delta / 2)
    out = [[], []]
    for group in groupedip:
        if method == 'weighted':
            x, y = weighted_average(*group)  # determine weighted mass and summed intensity
        elif method == 'centroid':
            x, y = centroid(group)
        out[0].append(x)
        out[1].append(y)
    maxint = max(out[1])
    for ind, val in enumerate(out[1]):
        out[0][ind] = out[0][ind]  # / abs(charge)
        out[1][ind] = val / maxint * 100.  # normalize to 100
    if verbose is True:
        sys.stdout.write(' DONE\n')
    return out


def normal_distribution(center, fwhm, height):
    """
    Generates a normal distribution about the center with the full width at half max specified. Y-values will be
    normalized to the height specified.

    :param center: center x value for the distribution
    :param fwhm: full-width-at-half-maximum
    :param height: maximum value for the y list
    :return: x values, y values
    :rtype: list
    """
    x = np.arange(
        center - fwhm * 2,
        center + fwhm * 2,
        10 ** -autodec(fwhm),
        dtype=np.float64,
    )
    y = stats.norm.pdf(  # generate normal distribution
        x,
        float(center),  # type-convert from sympy Float
        standard_deviation(fwhm),
    )
    y /= max(y)  # normalize
    y = y * height
    return [x.tolist(), y.tolist()]


def autodec(fwhm):
    """
    Automatically calculates the appropriate decimal place to track based on a full-width-at-half-maximum

    :param fwhm: full-width-at-half-maximum
    :return: decimal power
    :rtype: int
    """
    shift = fwhm
    n = 0
    while shift < 1.:
        n += 1
        shift = fwhm * 10 ** n
    return n + 1  # track 1 higher


def gaussian_isotope_pattern(
        barip: list,
        fwhm: float,
        verbose: bool = VERBOSE,
):
    """
    Simulates the isotope pattern that would be observed in a mass
    spectrometer with the resolution specified as a fwhm value.

    :param barip: The isotope pattern to be simulated. This can be either the bar isotope
        pattern or the raw isotope pattern (although this will be substantially
        slower for large molecules).
    :param fwhm: full-width-at-half-maximum
    :param verbose: chatty mode
    :return: The predicted gaussian isotope pattern in the form of a paired list
        ``[[m/z values],[intensity values]]``
    :rtype: list
    """
    spec = Spectrum(  # generate Spectrum object to encompass the entire region
        autodec(fwhm),
        start=min(barip[0]) - fwhm * 2,
        end=max(barip[0]) + fwhm * 2,
        empty=False,  # whether or not to use emptyspec
        filler=0.,  # fill with zeros, not None
    )
    for ind, val in enumerate(barip[0]):  # generate normal distributions for each peak
        # if verbose is True:
        #    sys.stdout.write('\rSumming m/z %.3f %d/%d' %(val,ind+1,len(self.barip[0])))
        nd = normal_distribution(val, fwhm, barip[1][ind])  # generate normal distribution for that peak
        spec.add_spectrum(nd[0], nd[1])  # add the generated spectrum to the total spectrum
    spec.normalize()  # normalize
    gausip = spec.trim()  # trim None values and output
    if verbose is True:
        sys.stdout.write(' DONE\n')
    return gausip


def isotope_pattern_hybrid(
        composition: dict,
        fwhm: float,
        decpl: int,
        verbose: bool = VERBOSE,
        dropmethod: str = None,
        threshold: float = THRESHOLD,
        npeaks: int = NPEAKS,
        consolidate: float = CONSOLIDATE,
        **kwargs,
):
    """
    A hybrid isotope pattern calculator which calculates the isotope pattern from each element, then multiplies the
    lists.

    :param composition: composition dictionary
    :param fwhm: full-width-at-half-maximum
    :param decpl: decimal places to track in the Spectrum object
    :param verbose: chatty mode
    :param dropmethod: optional method to use for low-intensity peak dropping or consolidation. Valid options are
        'threshold', 'npeaks', or 'consolidate'.
    :param threshold: if the dropmethod is set to 'threshold', any peaks below this threshold will be dropped.
    :param npeaks: if the dropmethod is set to 'npeaks', the top n peaks will be kept, with the rest being dropped.
    :param consolidate: if the dropmethod is set to 'consolidate', any peaks below the threshold will be consolidated
        into adjacent peaks using a weighted average. Any peaks that do not have a neighbour within 10^-`consolidate`
        will be dropped entirely.
    :return: isotope pattern as a Spectrum object
    :rtype: Spectrum
    """
    eleips = {}  # dictionary for storing the isotope patterns of each element
    for element, number in composition.items():
        eleips[element] = isotope_pattern_combinatoric(  # calculate the isotope pattern for each element
            {element: number},
            decpl=decpl,
            verbose=verbose,
        ).trim()  # trim the generated spectra to lists

    sortlist = []
    for element in eleips:
        sortlist.append((
            len(eleips[element][0]),
            element
        ))
    sortlist = sorted(sortlist)  # sorted list of elements based on the length of their isotope patterns
    sortlist.reverse()

    spec = None
    # todo convert to context tqdm (update string)
    for lenlist, element in tqdm(sortlist, desc='adding element to isotope pattern', disable=not verbose):
        if spec is None:
            spec = Spectrum(
                autodec(fwhm),  # decimal places
                start=None,  # minimum mass
                end=None,  # maximum mass
                empty=True,  # whether or not to use emptyspec
                filler=0.,  # fill with zeros, not None
                specin=eleips[element],  # supply masses and abundances as initialization spectrum
            )
            continue
        spec.add_element(eleips[element][0], eleips[element][1])
        spec.normalize(100.)  # normalize spectrum object
        if dropmethod == 'threshold':  # drop values below threshold
            spec.threshold(threshold)
        elif dropmethod == 'npeaks':  # keep top n number of peaks
            spec.keep_top_n(npeaks)
        elif dropmethod == 'consolidate':  # consolidate values being dropped
            spec.consolidate(
                threshold,
                3 * 10 ** -consolidate
            )
    return spec


class ReiterableCWR(object):
    def __init__(self, isos, number):
        """a reiterable version of combinations with replacements iterator"""
        self.isos = isos  # isotopes group
        self.number = number  # number of atoms of the element

    def __iter__(self):
        return cwr(self.isos, self.number)


@st.profilefn
def num_permu(lst, isos):
    """
    Calculates the number of unique permutations of the given set of isotopes for an element.
    The calculation is generated as a sympy function before evaluation. numpy factorial is limited in the size of
    factorials that are calculable, so sympy is required.

    :param lst: list of isotopes in the combination
    :param isos: possible isotopes for that element
    :return: number of occurrences of this list of isotopes
    :rtype: int
    """
    counts = [lst.count(x) for x in isos]  # counts the number of each isotope in the set
    num = sym.factorial(len(lst))  # numerator is the factorial of the length of the list
    denom = 1  # denominator is the product of the factorials of the counts of each isotope in the list
    for count in counts:
        denom *= sym.factorial(count)
    return int((num / denom).evalf())  # divide, evaluate, and return integer


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


def cpu_list_product(iterable):
    """returns the product of a list"""
    prod = 1
    for n in iterable:
        prod *= n
    return prod


@st.profilefn
def isotope_pattern_combinatoric(
        comp: dict,
        decpl: int,
        verbose: bool = VERBOSE,
        **kwargs,  # catch for extra keyword arguments
):
    """
    Calculates the raw isotope pattern of a given molecular formula with mass defects preserved.
    Uses a combinatorial method to generate isotope formulae

    :param comp: composition dictionary
    :param decpl: decimal places to track in the Spectrum object
    :param verbose: chatty mode
    :return: raw isotope pattern as a Spectrum object
    :rtype: Spectrum
    """
    speciso = False  # set state for specific isotope
    isos = {}  # isotopes dictionary
    isosets = {}  # set of isotopes for each element
    iterators = []  # list of iterators
    nk = []
    for element in comp:  # for each element
        if element in mass_dict:
            isosets[element] = []  # set of isotopes
            for isotope in mass_dict[element]:  # for each isotope of that element in the mass dictionary
                if isotope != 0 and mass_dict[element][isotope][1] != 0:  # of the intensity is nonzero
                    isosets[element].append(isotope)  # track set of isotopes
                    isos[isotope] = element  # create isotope,element association for reference
            iterators.append(
                ReiterableCWR(  # create iterator instance
                    isosets[element],
                    comp[element]
                )
            )
            if verbose is True:
                nk.append([  # track n and k for list length output
                    len(isosets[element]),
                    comp[element]
                ])
        else:  # if it's an isotope
            speciso = True

    spec = Spectrum(  # initiate spectrum object
        decpl,  # decimal places
        start=None,  # no minimum mass
        end=None,  # no maximum mass
        empty=True,  # whether or not to use emptyspec
        filler=0.,  # fill with zeros, not None
    )

    iterations = int(cpu_list_product([numberofcwr(n, k) for n, k in nk]))  # number of iterations

    for comb in tqdm(product(*iterators), desc='processing isotope combination', total=iterations, disable=not verbose):
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
                x += mass_dict[element][isotope][0]  # shift x
                y *= mass_dict[element][isotope][1]  # multiply intensity
        # add the x and y combination factored by the number of times that combination will occur
        spec.add_value(x, y * num)

    if speciso is True:  # if an isotope was specified
        for element in comp:
            if element not in mass_dict:  # if an isotope
                ele, iso = string_to_isotope(element)  # determine element and isotope
                spec.shift_x(mass_dict[ele][iso][0] * comp[element])  # shift the x values by the isotopic mass
    spec.normalize()  # normalize the spectrum object
    return spec


@st.profilefn
def isotope_pattern_multiplicative(
        comp: dict,
        decpl: int,
        verbose: bool = VERBOSE,
        dropmethod: str = None,
        threshold: float = THRESHOLD,
        npeaks: int = NPEAKS,
        consolidate: float = CONSOLIDATE,
        **kwargs,
):
    """
    Calculates the raw isotope pattern of a given molecular formula with mass defects preserved.

    :param comp: The molecular composition dictionary. See ``Molecule.composition`` for more details.
    :param decpl: The number of decimal places to track. This is normally controlled by the keyword
        arguments of the class, but can be specified if called separately.
    :param verbose: chatty mode
    :param dropmethod: optional method to use for low-intensity peak dropping or consolidation. Valid options are
        'threshold', 'npeaks', or 'consolidate'.
    :param threshold: if the dropmethod is set to 'threshold', any peaks below this threshold will be dropped.
    :param npeaks: if the dropmethod is set to 'npeaks', the top n peaks will be kept, with the rest being dropped.
    :param consolidate: if the dropmethod is set to 'consolidate', any peaks below the threshold will be consolidated
        into adjacent peaks using a weighted average. Any peaks that do not have a neighbour within 10^-`consolidate`
        will be dropped entirely.
    :return: Returns the isotope pattern with mass defects preserved (referred to as the 'raw'
        isotope pattern in this script).
    :rtype: Spectrum
    """
    spec = None  # initial state of spec
    if verbose is True:
        sys.stdout.write('Generating raw isotope pattern.\n')

    for key in comp:  # for each element
        if key in mass_dict:  # if not a single isotope
            masses = []  # list for masses of each isotope
            abunds = []  # list for abundances
            for mass in mass_dict[key]:
                if mass != 0:
                    if mass_dict[key][mass][1] > 0:  # if abundance is nonzero
                        masses.append(mass_dict[key][mass][0])
                        abunds.append(mass_dict[key][mass][1])
            msg = f'Processing element {key}'
            for n in tqdm(range(comp[key]), desc=msg, disable=not verbose):  # for n number of each element
                if spec is None:  # if spectrum object has not been defined
                    spec = Spectrum(
                        decpl,  # decimal places
                        start=min(masses) - 10 ** -decpl,  # minimum mass
                        end=max(masses) + 10 ** -decpl,  # maximum mass
                        specin=[masses, abunds],  # supply masses and abundances as initialization spectrum
                        empty=True,  # whether or not to use emptyspec
                        filler=0.,  # fill with zeros, not None
                    )
                    continue
                spec.add_element(masses, abunds)  # add the element to the spectrum object
                spec.normalize(100.)  # normalize spectrum
                if dropmethod == 'threshold':  # drop values below threshold
                    spec.threshold(threshold)
                elif dropmethod == 'npeaks':  # keep top n number of peaks
                    spec.keep_top_n(npeaks)
                elif dropmethod == 'consolidate':  # consolidate values being dropped
                    # todo figure out what's wrong here
                    raise NotImplementedError("There are bugs here, for the time being don't use the 'consolidate' "
                                              "dropmethod.")
                    spec.consolidate(
                        threshold,
                        3 * 10 ** -consolidate
                    )
        else:  # if specific isotope
            ele, iso = string_to_isotope(key)  # find element and isotope
            if spec is None:  # if spectrum object has not been defined
                spec = Spectrum(
                    decpl,  # decimal places
                    start=(mass_dict[ele][iso][0] * float(comp[key])) - 10 ** -decpl,  # minimum mass
                    end=(mass_dict[ele][iso][0] * float(comp[key])) + 10 ** -decpl,  # maximum mass
                    specin=[[mass_dict[ele][iso][0] * float(comp[key])], [1.]],
                    # supply masses and abundances as initialization spectrum
                    empty=True,  # whether or not to use emptyspec
                    filler=0.  # fill with zeros, not None
                )
                continue
            # todo add tqdm progress bar
            spec.shift_x(mass_dict[ele][iso][0])  # offset spectrum object by the mass of that
    spec.normalize()
    if verbose is True:
        sys.stdout.write('DONE\n')
    return spec


def isotope_pattern_isospec(
        comp: dict,
        decpl: int,
        verbose: bool = VERBOSE,
        threshold: float = THRESHOLD,
        **kwargs,
):
    """
    Generates a raw isotope pattern using the isospecpy package. http://matteolacki.github.io/IsoSpec/

    :param comp: composition dictionary
    :param decpl: decimal places to track while converting from isospec to Spectrum
    :param verbose: chatty mode
    :param threshold: threshold level (relative, seems slightly buggy)
    :param kwargs: catch for extra kwargs
    :return: Spectrum object
    """
    global _CITATION_REMINDER
    if _CITATION_REMINDER is False:  # remind the user on the first use
        print('IsoSpecPy package was used, please cite https://dx.doi.org/10.1021/acs.analchem.6b01459')
        _CITATION_REMINDER = True

    if any([key not in mass_dict for key in comp]):
        # todo see if there's a workaround for isotope specification
        raise KeyError(f'Isotope specification is not supported in IsoSpec calling. Please use a different isotope '
                       f'pattern generation method for isotopes. ')

    # todo see if there's a way to use IsoThresholdGenerator instead
    # use IsoSpec algorithm to generate configurations
    iso_spec = IsoThreshold(
        formula="".join(f'{ele}{num}' for ele, num in comp.items()),
        threshold=threshold * 0.1,
    )

    spec = Spectrum(
        decpl,  # decimal places
        start=min(iso_spec.masses) - 10 ** -decpl,  # minimum mass
        end=max(iso_spec.masses) + 10 ** -decpl,  # maximum mass
        empty=True,
        filler=0.  # fill with zeros, not None
    )
    # add values to Spectrum object
    for mass, abund in zip(iso_spec.masses, iso_spec.probs):
        spec.add_value(
            mass,
            abund
        )
    spec.normalize()  # normalize values to 100.
    return spec


def pattern_molecular_weight(mzs: list, intensities: list, charge: int = 1):
    """
    Calculates the molecular weight given by an isotope pattern.

    :param mzs: m/z (x) values for pattern
    :param intensities: intensity (y) values for the pattern
    :param charge: charge for the molecule
    :return: molecular weight
    :rtype: float
    """
    return sum([  # sum
        mz * intensity * charge  # of the product of the m/z, intensity, and charge
        for mz, intensity in zip(mzs, intensities)  # for all the values
    ]) / sum(intensities)  # divided by the sum of the intensities


def molecular_weight_error(calculated: float, expected: float):
    """
    Calculate the error between a calculated and expected molecular weight. This method may be used as a validation
    tool for calculated isotope patterns.

    :param calculated: calculated molecular weight (derived from an isotope pattern)
    :param expected: expected (true) molecular weight (derived from the molecular weights of the constituent elements)
    :return: Calculated error. Typically a difference of 3 parts per million (3*10^-6) is deemed an acceptable
        error.
    :rtype: float
    """
    return (calculated - expected) / expected


class Molecule(object):
    _comp = {}  # storage for composition of the molecule
    _mf = ''
    verbose = VERBOSE

    def __init__(self,
                 string: (str, dict),
                 charge=1,
                 mass_key='nist_mass',
                 verbose=False,
                 ):
        """
        Calculates many properties of a specified molecule.

        :param str, dict string: The molecule to interpret. A composition dictionary may also be specified here.
        :param int, str charge: the charge of the molecule (for mass spectrometric applications).
            This will affect any properties related to the mass to charge
            ratio. If the charge is specified in the input molecular formula, this will be
            overridden.
        :param str mass_key: The mass dictionary to use for calculations. Default is nist_mass, but additional mass
            dictionaries may be defined in the mass_dictionary file and retrieved using the dictionary name
            used to define them.
        :param bool verbose: Verbose output. Mostly useful when calculating for large molecules or while debugging.

        **Notes regarding string specification**

        - Common abbreviations may be predefined in mass_abbreviations.py (either locally or in the current working
            directory)

        - Use brackets to signify multiples of a given component (nested brackets are supported)

        - Isotopes may be specified using an isotope-element format within a bracket (e.g. carbon 13 would be specified
            as "(13C)" ). The mass of that isotope must be defined in the mass dictionary being used by the script
            (default NIST mass).

        - The charge may be specified in the formula, but care must be taken here. Charge must be specified in either
            sign-value (e.g. '+2') or within a bracket. Otherwise, the script may attempt to interpret the charge as a
            magnitude specifier of the previous block or as an isotope, and errors will be encountered.

        - A composition dictionary with the format `{'Element': number_of_that_element, ...}` may be provided instead
            of a string formula

        """
        if verbose is True:
            sys.stdout.write(f'Generating molecule object from input {string}\n')
        # split charge into value and sign
        self.charge, self.sign = interpret_charge(charge)
        self.mass_key = mass_key  # store mass dictionary that the script will use
        self.verbose = verbose
        if type(string) == dict:  # if a composition dictionary was provided
            self.composition = string
        elif type(string) == str:  # set string and interpret formula
            self.molecular_formula = string
        else:
            raise TypeError(f'The provided string type is not interpretable: {type(string)}')

        if self.verbose is True:
            self.print_details()

    def __repr__(self):
        return f'{self.__class__.__name__}({self.molecular_formula})'

    def __str__(self):
        return self.__repr__()

    def __contains__(self, item):
        if type(item) == str:
            return item in self._comp
        elif type(item) == list or type(item) == tuple:
            return all([element in self._comp for element in item])
        elif type(item) == dict:
            return all([
                element in self._comp and self._comp[element] >= num for element, num in item.items()
            ])
        elif isinstance(item, Molecule):
            return self.__contains__(item.composition)
        else:
            raise TypeError(f'The item {item} is not a recognized type for containment checks. Type: {type(item)}')

    def __iter__(self):
        for element in self._comp:
            yield element

    def __getitem__(self, item):
        return self._comp[item]

    def __eq__(self, other):
        if type(other) == dict:
            return other == self._comp
        elif isinstance(other, Molecule):
            return other.composition == self._comp
        return False

    def __ne__(self, other):
        return not self.__eq__(other)

    def __lt__(self, other):
        if type(other) == dict:
            return all([
                number < self._comp[element] for element, number in other.items()
            ])
        elif isinstance(other, Molecule):
            return all([
                number < self._comp[element] for element, number in other.composition.items()
            ])
        else:
            raise TypeError(f'Comparison of type {type(other)} to {self.__class__} is unsupported. ')

    def __le__(self, other):
        return self.__eq__(other) or self.__lt__(other)

    def __gt__(self, other):
        if type(other) == dict:
            return all([
                number > self._comp[element] for element, number in other.items()
            ])
        elif isinstance(other, Molecule):
            return all([
                number > self._comp[element] for element, number in other.composition.items()
            ])
        else:
            raise TypeError(f'Comparison to type {type(other)} to {self.__class__} is unsupported. ')

    def __ge__(self, other):
        return self.__eq__(other) or self.__gt__(other)

    def __getinitargs__(self):
        return (
            self.composition,
            f'{self.charge}{self.sign}',
            self.mass_key,
            self.verbose,
        )

    def __reduce__(self):
        """pickle support"""
        return (
            self.__class__,
            self.__getinitargs__(),
        )

    def __add__(self, other):
        """
        Several supported addition methods:
        If a valid molecular formula string is provided, that string will be added.
        If another Molecule class instance is provided, the provided instance will be
        added to the current instance.
        """
        if type(other) is str:
            other = composition_from_formula(other)
        elif isinstance(other, Molecule) is True:
            other = other.composition
        elif type(other) == dict:
            pass
        else:
            raise ValueError(f'Addition of {other} to {self} is invalid')
        new = copy.copy(self._comp)  # starter for new dictionary

        for key in other:
            try:
                new[key] += other[key]
            except KeyError:
                new[key] = other[key]
        return self.__class__(
            new,
            charge=f'{self.charge}{self.sign}'
        )

    def __radd__(self, other):
        return self.__add__(other)

    def __iadd__(self, other):
        if type(other) is str:
            other = composition_from_formula(other)
        elif isinstance(other, Molecule) is True:
            other = other.composition
        elif type(other) == dict:
            pass
        else:
            raise ValueError(f'Addition of {other} to {self} is invalid')
        new = copy.copy(self._comp)  # starter for new dictionary
        for key in other:
            try:
                new[key] += other[key]
            except KeyError:
                new[key] = other[key]
        self.composition = new
        return self

    def __sub__(self, other):
        """
        See __add__ for details.
        Subtract has a catch for a negative number of a given element
        (the minimum that can be reached is zero).
        """
        if type(other) is str:
            other = composition_from_formula(other)
        elif isinstance(other, Molecule) is True:
            other = other.composition
        elif type(other) == dict:
            pass
        else:
            raise ValueError(f'Addition of {other} to {self} is invalid')
        new = copy.copy(self._comp)  # starter for new dictionary

        for key in other:
            if new[key] - other[key] < 0 or key not in new:
                raise ValueError('Subtraction of {other[key]} {key} from {self} would yield a negative number of that '
                                 'element.')
            new[key] -= other[key]
        return self.__class__(
            new,
            charge=f'{self.charge}{self.sign}'
        )

    def __rsub__(self, other):
        return self.__sub__(other)

    def __isub__(self, other):
        if type(other) is str:
            other = composition_from_formula(other)
        elif isinstance(other, Molecule) is True:
            other = other.composition
        elif type(other) == dict:
            pass
        else:
            raise ValueError(f'Addition of {other} to {self} is invalid')
        new = copy.copy(self._comp)  # starter for new dictionary

        for key in other:
            if new[key] - other[key] < 0 or key not in new:
                raise ValueError('Subtraction of {other[key]} {key} from {self} would yield a negative number of that '
                                 'element.')
            new[key] -= other[key]
        self.composition = new
        return self

    def __mul__(self, other):
        """allows integer multiplication of the molecular formula"""
        if type(other) != int:
            raise ValueError(f'Non-integer multiplication of a {self.__class__.__name__} object is unsupported')
        new = copy.copy(self._comp)  # starter for new dictionary
        for key in new:
            new[key] = new[key] * other
        return self.__class__(
            new,
            charge=f'{self.charge}{self.sign}'
        )

    def __rmul__(self, other):
        return self.__mul__(other)

    def __imul__(self, other):
        if type(other) != int:
            raise ValueError(f'Non-integer multiplication of a {self.__class__.__name__} object is unsupported')
        new = copy.copy(self._comp)  # starter for new dictionary
        for key in new:
            new[key] = new[key] * other
        self.composition = new
        return self

    def __truediv__(self, other):
        """allows integer division of the molecular formula"""
        if type(other) != int:
            raise ValueError(f'Non-integer division of a {self.__class__.__name__} object is unsupported')
        new = copy.copy(self._comp)  # starter for new dictionary
        for key in new:
            newval = new[key] / other
            if newval.is_integer() is False:
                raise ValueError(f'Division of {new[key]} {key} by {other} yielded a non-integer number {newval}')
            new[key] = int(newval)
        return self.__class__(
            new,
            charge=f'{self.charge}{self.sign}'
        )

    def __itruediv__(self, other):
        if type(other) != int:
            raise ValueError(f'Non-integer division of a {self.__class__.__name__} object is unsupported')
        new = copy.copy(self._comp)  # starter for new dictionary
        for key in new:
            newval = new[key] / other
            if newval.is_integer() is False:
                raise ValueError(f'Division of {new[key]} {key} by {other} yielded a non-integer number {newval}')
            new[key] = int(newval)
        self.composition = new
        return self

    @property
    def composition(self):
        """Composition dictionary"""
        return self._comp

    @composition.setter
    def composition(self, dct):
        if type(dct) != dict:
            raise TypeError('The composition must be a dictionary')
        dct = copy.copy(dct)
        dct = abbreviations(dct)  # check for and convert abbreviations
        if 'charge' in dct:  # if charge was specified in the formula
            self.charge, self.sign = interpret_charge(dct['charge'])
            del dct['charge']
        check_in_mass_dict(dct)  # check in mass dictionary
        self._comp = dct  # set local dictionary

    @property
    def molecular_formula(self):
        """Molecular formula of the molecule"""
        out = ''
        # todo catch carbon and hydrogen isotopes first
        if 'C' in self.composition:  # carbon and hydrogen first according to hill formula
            out += f'C{self.composition["C"]}' if self.composition['C'] > 1 else 'C'
        if 'H' in self.composition:
            out += f'H{self.composition["H"]}' if self.composition['H'] > 1 else 'H'
        for key, val in sorted(self.composition.items()):  # alphabetically otherwise
            if key != 'C' and key != 'H':
                if key in mass_dict:
                    out += f'{key}{self.composition[key]}' if self.composition[key] > 1 else f'{key}'
                else:  # if an isotope
                    ele, iso = string_to_isotope(key)
                    out += f'({iso}{ele})'
                    out += f'{self.composition[key]}' if self.composition[key] > 1 else ''
        return out

    @molecular_formula.setter
    def molecular_formula(self, formula):
        self.composition = composition_from_formula(formula)
        self._mf = formula

    @property
    def molecular_formula_formatted(self):
        """returns the subscript-formatted molecular formula"""
        out = ''
        if 'C' in self.composition:
            out += f'C{to_subscript(self.composition["C"]) if self.composition["C"] > 1 else "C"}'
        if 'H' in self.composition:
            out += f'H{to_subscript(self.composition["H"]) if self.composition["H"] > 1 else "H"}'
        for key, val in sorted(self.composition.items()):
            if key not in ['C', 'H']:
                if key in mass_dict:
                    out += f'{key}{to_subscript(self.composition[key])}' if self.composition[key] > 1 else f'{key}'
                else:
                    ele, iso = string_to_isotope(key)
                    out += f'{to_superscript(iso)}{ele}'
                    out += f'{to_subscript(self.composition[key])}' if self.composition[key] > 1 else ''
        return out

    @property
    def sf(self):
        """legacy catch for shorthand 'string formula' attribute"""
        return self.molecular_formula

    @property
    def molecular_weight(self):
        """Molecular weight of the molecule"""
        mwout = 0
        for element, number in self.composition.items():
            try:
                mass = mass_dict[element]
                for isotope in mass:
                    if isotope == 0:
                        continue
                    # add every isotope times its natural abundance times the number of that element
                    mwout += mass[isotope][0] * mass[isotope][1] * number
            except KeyError:  # if isotope
                ele, iso = string_to_isotope(element)
                mwout += mass_dict[ele][iso][0] * number  # assumes 100% abundance if specified
        return mwout

    @property
    def mw(self):
        """legacy catch for shorthand molecular weight"""
        return self.molecular_weight

    @property
    def percent_composition(self):
        """Elemental percent composition of the molecule"""
        pcompout = {}  # percent composition dictionary
        for element, number in self.composition.items():
            try:
                mass = mass_dict[element]
                for isotope in mass:
                    if isotope == 0:
                        continue
                    if element not in pcompout:
                        pcompout[element] = 0.
                    # add mass contributed by that element
                    pcompout[element] += mass[isotope][0] * mass[isotope][1] * number
            except KeyError:  # if isotope
                ele, iso = string_to_isotope(element)
                pcompout[str(iso) + ele] = mass_dict[ele][iso][0] * number
        mw = self.molecular_weight
        for element in pcompout:  # determines the percent composition of each element
            try:
                pcompout[element] = pcompout[element] / mw
            except ZeroDivisionError:
                pcompout[element] = 0.
        return pcompout

    @property
    def pcomp(self):
        """legacy catch for shorthand percent composition"""
        return self.percent_composition

    @property
    def monoisotopic_mass(self):
        """An estimation of the exact mass given by the molecular formula. This is likely not accurate for high-mass
        species"""
        em = 0.
        for element, number in self.composition.items():
            try:
                em += mass_dict[element][0][0] * number
            except KeyError:
                ele, iso = string_to_isotope(element)
                em += mass_dict[ele][iso][0] * number
        # # accounts for the mass of an electron (uncomment if this affects your data)
        # if self.sign == '+':
        #    em -= (9.10938356*10**-28)*charge
        # if self.sign == '-':
        #    em += (9.10938356*10**-28)*charge
        return em / self.charge

    @property
    def standard_deviation_comp(self):
        """
        cacluates the standard deviation of the isotope pattern of the supplied composition
        this calculation is based on Rockwood and Van Orden 1996 doi: 10.1021/ac951158i
        """
        stdev = 0.
        for element, number in self.composition.items():
            meanmass = 0
            eledev = 0  # elemental deviation
            mass = mass_dict[element]
            for isotope in mass:  # calculate weighted average mass
                if isotope != 0:
                    meanmass += sum(mass[isotope])  # weighted average mass
            for isotope in mass:
                if mass != 0:
                    eledev += mass[isotope][1] * (mass[isotope][0] - meanmass) ** 2
            stdev += eledev * number
        return np.sqrt(stdev)

    def print_details(self):
        """prints the details of the generated molecule"""
        sys.stdout.write(f'{self}\n')
        sys.stdout.write(f'formula: {self.molecular_formula}\n')
        sys.stdout.write(f'molecular weight: {round(self.molecular_weight, 6)}\n')
        sys.stdout.write(f'monoisotopic mass: {round(self.monoisotopic_mass, 6)}\n')
        sys.stdout.write('\n')
        self.print_percent_composition()

    def print_percent_composition(self):
        """prints the percent composition in a reader-friendly format"""
        sys.stdout.write('elemental percent composition:\n')
        pcomp = self.percent_composition
        for element, percent in sorted(pcomp.items()):
            sys.stdout.write(f'{element}: {percent * 100.:6.4}%\n')


class IPMolecule(Molecule):
    _ipmethod = None
    _gausip = None  # gaussian isotope pattern storage
    _dropmethod = None

    def __init__(self,
                 string: (str, dict),
                 charge=1,
                 consolidate=3,
                 criticalerror=3 * 10 ** -6,
                 decpl=7,
                 dropmethod=None,
                 emptyspec=True,
                 groupmethod='weighted',
                 ipmethod='isospec',
                 keepall=False,
                 npeaks=5000,
                 resolution=5000,
                 threshold=0.01,
                 save=False,
                 verbose=VERBOSE,
                 precalculated=None,
                 ):
        """
        A class with many mass-spectrometric properties such as estimated exact masses, isotope patterns, error
        estimators, and basic plotting tools.

        :param str string: the molecule name to interpret. See Molecule documentation for more details
        :param int, str charge: the charge of the molecule (for mass spectrometric applications).
            This will affect any properties related to the mass to charge
            ratio. If the charge is specified in the input molecular formula, this will be
            overridden.

        :param int, float resolution: The resolution of the instrument to simulate when generating the gaussian isotope
            pattern. This also affects the bounds attribute.

        :param int consolidate: When using the consolidate drop method, consolidate peaks within 10^-*consolidate*
            of each other. See *dropmethod* for more details.

        :param float criticalerror:
            The critical error value used for warning the user of a potential calculation error.
            This only affects the ``print_details()`` function output. Default 3*10**-6 (3 parts per million)

        :param int decpl: The number of decimal places to track while calculating the isotope pattern.
            Decreasing this will improve efficiency but decrease accuracy. Options: integer.

        :param 'threshold', 'npeaks', 'consolidate' dropmethod: The peak drop method to use if desired.
            Using a peak dropping method will improve calculation times, but decrease the accuracy of the
            calculated isotope pattern. 'threshold' drops all peaks below a specified threshold value (specified using
            the *threshold* keyword argument). 'npeaks' keeps the top *n* peaks, specified by the *npeaks* keyword
            argument. 'consolidate' combines the intensity of peaks below the threshold value into the
            nearest peak (within the delta specified by the *consolidate* keyword argument, this method is the most
            accurate). The new peak *m/z* value is determined by the weighted average of the combined peaks. This will
            be repeated until the peak is above the threshold or there are no near peaks.

        :param bool emptyspec: Whether to use an empty spectrum obect. Disable this for very large molecules to
            improve calculation time.

        :param 'weighted', 'centroid' groupmethod: The grouping method to use when calculating the bar isotope pattern
            from the raw isotope pattern. Weighted calculates the peak locations using the weighted average of the *m/z*
            and intensity values. Centroid finds the center *m/z* value of a group of peaks.

        :param 'multiplicative', 'combinatorial', 'hybrid', 'cuda', ipmethod: The method to use for determining the isotope
            pattern. 'multiplicative' multiplies the existing list of intensities by each element. 'combinatorial' uses
            combinatorics and iterators to calculate each possible combination. 'hybrid' uses combinatorics to calcuate
            the pattern from each element, then multiplies those together

        :param bool keepall: Whether to keep all peaks calculated in the isotope pattern. When false, this will drop
            all intensities below 0.0001 after calculating the isotope pattern.

        :param int npeaks: The number of peaks to keep if *dropmethod* is 'npeaks'. See *dropmethod* for more details.

        :param float threshold: The threshold value determining whether or not to drop a peak. Only has an effect if
            *dropmethod* is not ``None``. See *dropmethod* for more details.

        :param bool verbose: Verbose output. Mostly useful when calculating for large molecules or while debugging.

        """
        # todo implement apply_threshold method for trimming resulting spectrum
        self.ipmethod = ipmethod
        self._spectrum_raw = None  # spectrum object holder
        self._raw = None  # raw isotope pattern
        self.bar_isotope_pattern = [[], []]
        self.criticalerror = criticalerror
        self.decpl = decpl
        self.dropmethod = dropmethod
        self.emptyspec = emptyspec
        self.consolidate = consolidate
        self.groupmethod = groupmethod
        self.keepall = keepall
        self.npeaks = npeaks
        self.resolution = resolution
        self.threshold = threshold
        self.save = save  # todo reimplement and detail in docstring

        if precalculated is not None:  # if precalculated values were provided, pull and set to prevent recalculation
            self._comp = precalculated['composition']
            self._spectrum_raw = precalculated['spectrum']
            self.bar_isotope_pattern = precalculated['barip']
            self._raw = precalculated['rawip']
            self._gausip = precalculated['gausip']

        Molecule.__init__(
            self,
            string,
            charge,
            verbose=verbose,
        )

        if save is True:
            self.save_to_jcamp()

    def __reduce__(self):
        return (
            self.__class__,
            self.__getinitargs__(),
        )

    def __getinitargs__(self):
        return (
            self.composition,
            self.charge,
            self.consolidate,
            self.criticalerror,
            self.decpl,
            self.dropmethod,
            self.emptyspec,
            self.groupmethod,
            self.ipmethod,
            self.keepall,
            self.npeaks,
            self.resolution,
            self.threshold,
            self.save,
            self.verbose,
            {  # precalculated values
                'composition': self.composition,
                'spectrum': self.spectrum_raw,
                'rawip': self.raw_isotope_pattern,
                'barip': self.bar_isotope_pattern,
                'gausip': self.gaussian_isotope_pattern if self._gausip is not None else None,
            },
        )

    @property
    def ipmethod(self):
        return self._ipmethod

    @ipmethod.setter
    def ipmethod(self, value):
        if value not in VALID_IPMETHODS:
            raise ValueError(f'The isotope pattern generation method {value} is not valid. ipmethod must be one '
                             f'of: {", ".join(VALID_IPMETHODS)}')
        self._ipmethod = value

    @property
    def dropmethod(self):
        return self._dropmethod

    @dropmethod.setter
    def dropmethod(self, value):
        if value not in VALID_DROPMETHODS:
            raise ValueError(f'The intensity dropping method {value} is not valid. dropmethod must be one '
                             f'of: {", ".join(VALID_DROPMETHODS)}')
        self._dropmethod = value

    @property
    def estimated_exact_mass(self):
        """determines the precise exact mass from the bar isotope pattern"""
        ind = self.bar_isotope_pattern[1].index(
            max(self.bar_isotope_pattern[1])
        )
        return self.bar_isotope_pattern[0][ind]

    @property
    def em(self):
        """Legacy attribute access: estimated exact mass"""
        return self.estimated_exact_mass

    @property
    def molecular_weight_estimated(self):
        """The molecular weight of the molecule estimated by the isotope pattern"""
        return pattern_molecular_weight(
            *self.raw_isotope_pattern,
            charge=self.charge,
        )

    @property
    def pmw(self):
        """Legacy retrieval of pattern molecular weight"""
        return self.molecular_weight_estimated

    @property
    def error(self):
        """Error of the generated isotope pattern"""
        return molecular_weight_error(
            calculated=self.molecular_weight_estimated,
            expected=self.molecular_weight,
        )

    @property
    def sigma(self):
        """Standard deviation of the isotope pattern"""
        return standard_deviation(self.fwhm)

    @property
    def nominal_mass(self):
        """the nominal mass of the molecule"""
        return int(round(self.em))

    @property
    def fwhm(self):
        try:  # try to return from estimated, unless uncalculated, use monoisotopic
            return self.estimated_exact_mass / self.resolution
        except (IndexError, ValueError):
            return self.monoisotopic_mass / self.resolution

    @property
    def barip(self):
        """Legacy attribute access"""
        return self.bar_isotope_pattern

    @property
    def raw_isotope_pattern(self):
        if self._raw is None:
            self._raw = self.spectrum_raw.trim()
        return self._raw

    @property
    def rawip(self):
        """Legacy attribute access"""
        return self.raw_isotope_pattern

    @property
    def spectrum_raw(self):
        return self._spectrum_raw

    @property
    def gaussian_isotope_pattern(self):
        if self._gausip is None:  # if it hasn't been calculated, generate
            self._gausip = gaussian_isotope_pattern(
                self.bar_isotope_pattern,
                self.fwhm
            )
        return self._gausip

    @property
    def gausip(self):
        """Legacy retrieval"""
        return self.gaussian_isotope_pattern

    @property
    def composition(self):
        return self._comp

    @composition.setter
    def composition(self, dct):
        if type(dct) != dict:
            raise TypeError('The composition must be a dictionary')
        if dct == self.composition:  # do nothing if the composition dictionary is the same as current
            return
        dct = copy.copy(dct)
        dct = abbreviations(dct)  # check for and convert abbreviations
        if 'charge' in dct:  # if charge was specified in the formula
            self.charge, self.sign = interpret_charge(dct['charge'])
            del dct['charge']
        check_in_mass_dict(dct)  # check in mass dictionary
        self._comp = dct  # set local dictionary
        self._calculate_ips()  # calculate isotope patterns
        # todo save to pickle

    @property
    def isotope_pattern_standard_deviation(self):
        """
        Cacluates the standard deviation of the isotope pattern of the supplied composition
        this calculation is based on Rockwood and Van Orden 1996 doi: 10.1021/ac951158i
        """
        return np.sqrt(
            sum([
                intensity * (mz - self.pmw) ** 2  # weighted distance from the estimated molecular weight
                for mz, intensity in zip(*self.raw_isotope_pattern)
            ])
        )

    @property
    def bounds(self):
        """Convenient attribute access to default bounds. Call calculate_bounds for additional options. """
        return self.calculate_bounds()

    @property
    def per_peak_bounds(self):
        """Convenient attribute access to per-peak bounds. Call calculate_bounds for additional options. """
        return self.calculate_bounds(perpeak=True)

    def calculate_bounds(
            self,
            conf: float = 0.95,
            perpeak: bool = False,
            threshold: float = 0.01
    ):
        """
        Calculates the *m/z* bounds of the isotope pattern of the molecule object based
        on a confidence interval and the *m/z* values of the bar isotope pattern.
        This can be used to automatically determine the integration bounds required to
        contain XX% of the counts associated with that molecule in a mass spectrum.

        :param conf: The confidence interval to use for calculating the bounds.
            e.g. *0.95* corresponds to a 95% confidence interval.
        :param perpeak: Whether or not to return the bounds required to integrate each
            peak of the isotope pattern individually.
            This can be useful in a very noisy mass spectrum to avoid
            baseline noise within the integration interval.
        :param threshold: The threshold used to determine whether a peak should be
            included in the bounds.
        :return: bounds.
            If *perpeak* is False, this will return a two item list
            corresponding to the start and end *m/z* bounds.
            If *perpeak* is True, returns a dictionary of bounds with
            the key format of
            ``dict[parent m/z value]['bounds'] = [start m/z, end m/z]``

        **Examples**

        To determine the integration bounds of C61H51IP3Pd: 

        ::

            >>> mol = IPMolecule('C61H51IP3Pd')
            >>> mol.calculate_bounds(0.95)
            [1104.9458115053008, 1116.3249999321531]

            >>> mol.calculate_bounds(0.99)
            [1104.8877964620444, 1116.3830149754094]

            >>> mol.calculate_bounds(0.95, True)
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
        if self.verbose is True:
            sys.stdout.write('Calculating bounds from simulated gaussian isotope pattern')
        threshold = threshold * max(self.bar_isotope_pattern[1])
        tempip = [[], []]
        for ind, inten in enumerate(self.bar_isotope_pattern[1]):  # checks for intensities above threshold
            if inten >= threshold:
                tempip[0].append(self.bar_isotope_pattern[0][ind])
                tempip[1].append(self.bar_isotope_pattern[1][ind])
        if perpeak is True:  # if per-peak bounds are called for
            out = {}
            for mz in tempip[0]:
                out[str(mz)] = {}
                out[str(mz)]['bounds'] = stats.norm.interval(conf, mz, scale=self.sigma)
        else:  # a general range that covers the entire isotope pattern
            out = [stats.norm.interval(conf, tempip[0][0], scale=self.sigma)[0],
                   stats.norm.interval(conf, tempip[0][-1], scale=self.sigma)[1]]
        if self.verbose is True:
            if perpeak is False:
                sys.stdout.write(': %.3f-%.3f' % (out[0], out[1]))
            sys.stdout.write(' DONE\n')
        return out

    def _calculate_ips(self):
        """Call to calculate isotope patterns based on the specified parameters"""
        # generates the raw isotope pattern (charge of 1)
        if self.ipmethod == 'combinatorics':
            calculator = isotope_pattern_combinatoric
        elif self.ipmethod == 'multiplicative':
            calculator = isotope_pattern_multiplicative
        elif self.ipmethod == 'hybrid':
            calculator = isotope_pattern_hybrid
        # elif self.ipmethod == 'cuda':
        #     calculator = isotope_pattern_cuda
        elif self.ipmethod == 'isospec':
            calculator = isotope_pattern_isospec
        else:
            raise ValueError(f'The isotope pattern method {self.ipmethod} is not valid')

        self._spectrum_raw = calculator(
            self.composition,
            decpl=self.decpl,
            verbose=self.verbose,
            dropmethod=self.dropmethod,
            threshold=self.threshold,
            npeaks=self.npeaks,
            consolidate=self.consolidate,
            fwhm=self.fwhm,
        )

        # apply charge
        self.spectrum_raw.charge = self.charge

        # generate bar isotope pattern based on the raw pattern
        self.bar_isotope_pattern = bar_isotope_pattern(
            self.raw_isotope_pattern,
            self.fwhm
        )

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
        # TODO fix this method (worthwhile?)
        #   - 2015-09-15 06 gives a bounds error
        yvals = []
        res = []
        maxy = float(max(exp[1]))
        if maxy == 0.:
            return 'could not calculate'
        for ind, val in enumerate(exp[1]):  # normalize y values
            yvals.append(float(val) / maxy * 100.)
        # avgy = sum(exp[1])/len(exp[1])
        for ind, mz in enumerate(exp[0]):
            if min(self.gausip[0]) < mz < max(self.gausip[0]):  # if within isotope pattern
                nspind = self.spectrum_raw.index(mz)  # calculate index
                if self.spectrum_raw.y[nspind] is not None:  # if the predicted intensity is not None
                    # difference between observed and predicted (residuals)
                    res.append(yvals[ind] - self.spectrum_raw.y[nspind])
                    # tot.append(self.spec.y[nspind]-avgy) # difference between predicted and mean
        # rsqrd = 1-(sumsquare(res)/sumsquare(tot)) # r-squared value (apparently not applicable to non-linear fits)
        return np.sqrt(sumsquare(res) / len(res))

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

    def load_from_pickle(self, customfile=None):
        """loads data from pickle"""
        raise NotImplementedError('This functionality has been temporarily disabled due to significant changes in the '
                                  'class. ')
        # TODO specify hierachy and pull if better method than specified
        if customfile is None:  # if no directory was specified, use current working directory
            customfile = os.path.join(
                os.getcwd(),
                'molecules',
                self.molecular_formula(self.comp) + '.mol',
            )
        if os.path.isfile(customfile) is True:
            if self.ipmethod.lower() == 'multiplicative':
                key = 'multiplicative'
            elif self.ipmethod.lower() == 'combinatorics':
                key = 'combinatorics'
            if self.dropmethod is not None:
                key += ' %s' % self.dropmethod
            subkey = self.decpl  # decimal places
            with open(customfile, 'rb') as targetfile:
                incoming = pickle.load(targetfile)
                if key in incoming and subkey in incoming[key]:
                    items = incoming[key][subkey]
                    strcharge = '%s%d' % (self.sign, self.charge)
                    if items['charge'] == strcharge:  # if the charge combination matches
                        print('Loading data from saved file %s' % customfile)
                        self.bar_isotope_pattern = items['bar isotope pattern']
                        self.raw_isotope_pattern = items['raw isotope pattern']
                        self.gausip = items['gaussian isotope pattern']
                        self.mw = items['mw']
                        self.mimass = items['monoisotopic mass']
                        self.em = items['estimated exact mass']
                        self.pcomp = items['percent composition']
                        self.error = items['error']
                        self.fwhm = items['full width at half max']
                        self.sigma = items['standard deviation']
                        self.sf = self.molecular_formula(self.comp)
                        return True
        return False  # if the exact match was not found, False

    def print_details(self):
        """prints the details of the generated molecule"""
        sys.stdout.write(f'{self}\n')
        sys.stdout.write(f'formula: {self.molecular_formula}\n')
        sys.stdout.write(f'molecular weight: {round(self.molecular_weight, self.decpl)}\n')
        sys.stdout.write(f'monoisotopic mass: {round(self.monoisotopic_mass, self.decpl)}\n')
        sys.stdout.write(f'estimated exact mass: {round(self.estimated_exact_mass, self.decpl)}\n')
        sys.stdout.write(f'error: {self.error:.3}\n')
        if abs(self.error) > self.criticalerror:
            sys.stdout.write(f'WARNING: Error is greater than {self.criticalerror}!\n')
        sys.stdout.write('\n')
        self.print_percent_composition()

    def plot_bar_pattern(self):
        """plots and shows the isotope bar pattern"""
        fwhm = self.em / self.resolution
        pl.bar(self.bar_isotope_pattern[0], self.bar_isotope_pattern[1], width=fwhm, align='center')
        pl.xlabel('m/z', style='italic')
        pl.ylabel('normalized intensity')
        pl.ticklabel_format(useOffset=False)
        pl.show()

    def plot_gaussian_pattern(self, exp=None):
        """plots and shows the simulated gaussian isotope pattern"""
        pl.plot(*self.gaussian_isotope_pattern, linewidth=1)
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
        pl.bar(self.raw_isotope_pattern[0], self.raw_isotope_pattern[1], width=self.fwhm)
        pl.xlabel('m/z', style='italic')
        pl.ylabel('normalized intensity')
        pl.ticklabel_format(useOffset=False)
        pl.show()

    def save_to_jcamp(self, name=None):
        """
        Saves the bar isotope pattern to JCAMP-DX file format
        Output type roughly based on the output from ChemCalc.org
        see http://www.jcamp-dx.org/protocols.html for details on the JCAMP-DX specifications.

        :param name: optional name for the output file (default is {molecular formula}.jdx)
        """
        if os.path.isdir(os.path.join(os.getcwd(), 'molecules')) is False:
            os.makedirs(os.path.join(os.getcwd(), 'molecules'))
        if name is None:  # if no name supplied, auto generate
            name = self.molecular_formula
            name += '.jdx'
        elif name.lower().endswith('.jdx') is False:
            name += '.jdx'

        if self.verbose is True:
            sys.stdout.write(f'Saving {name} to {os.path.join(os.getcwd(), "molecules")}')
            sys.stdout.flush()

        header = [  # comment lines to put before data
            # header items
            f'TITLE= {self.molecular_formula}',
            'JCAMP-DX= 5.01',
            'DATA TYPE= MASS SPECTRUM',
            'DATA CLASS= PEAK TABLE',
            f'ORIGIN= Calculated spectrum from PythoMS {self.__class__} class https://github.com/larsyunker/PythoMS',
            f'OWNER= {os.getlogin()}',
            f'SPECTROMETER/DATA SYSTEM= {self.__class__} class {self.ipmethod} method',
            f'CREATION DATE= {datetime.now().astimezone()}',
            'XUNITS= M/Z',
            'YUNITS= RELATIVE ABUNDANCE',
            f'NPOINTS= {len(self.bar_isotope_pattern[0])}',
            f'FIRSTX= {self.bar_isotope_pattern[0][0]}',
            f'LASTX= {self.bar_isotope_pattern[0][1]}',

            # user defined labels
            f'$Molecular weight= {self.molecular_weight}',
            f'$Resolution= {self.res}',
            f'$Threshold= {self.threshold if self.threshold is not None else ""}',
            f'$Error= {self.error:.2}',
            f'$Nominal mass = {self.nominal_mass}',
            f'$Monoisotopic mass= {self.monoisotopic_mass}',
            f'$Estimated exact mass= {self.estimated_exact_mass}',
        ]
        with open(os.path.join(os.getcwd(), "molecules", name), 'wt') as outfile:
            for line in header:  # write header lines
                if len(line) != 0:
                    outfile.write(f'##{line}\n')
            outfile.write('##PEAK TABLE= (XY..XY)\n')
            for mz, intensity in zip(*self.bar_isotope_pattern):  # write data lines
                outfile.write(f'{mz}, {intensity}\n')
            outfile.write('##END=\n')

    def save_to_pickle(self, name=None):
        """
        Saves the molecule's properties to pickle
        """
        if os.path.isdir(os.path.join(os.getcwd(), 'molecules')) is False:
            os.makedirs(os.path.join(os.getcwd(), 'molecules'))
        if name is None:  # if no name supplied, auto generate
            name = self.molecular_formula
            name += '.mol'
        elif name.lower().endswith('.mol') is False:
            name += '.mol'

        if self.verbose is True:
            sys.stdout.write(f'Saving {name} to {os.path.join(os.getcwd(), "molecules")}')
            sys.stdout.flush()

        with open(os.path.join(os.getcwd(), "molecules", name), 'wb') as outfile:
            pickle.dump(
                self,
                outfile
            )

        # todo differentiate between generation methods in the output files


if __name__ == '__main__':  # for testing and troubleshooting
    # st.printstart()
    mol = IPMolecule(
        'L2PdAr+I',
        # charge= 2, # specify charge (if not specified in formula)
        # res=1050000, # specify spectrometer resolution (default 5000)
        verbose=True,
        # decpl=10,
        # dropmethod='threshold',
        # threshold=0.00001,
        # ipmethod='hybrid',
        ipmethod='combinatorics',
        # keepall=True,
    )
    # mol.print_details()
    # st.printelapsed()
    # st.printprofiles()
