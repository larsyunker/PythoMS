# molecule module

The `molecule` module contains two primary classes: `Molecule` and
`IPMolecul`. These are classes which automatically calculate a variety
of molecular properties of provided molecular formulae. The `Molecule`
class provides has an interpreter for converting a string formula to
a molecular formula, and provides access to basic attributes such as
molecular weight. The `IPMolecule` class contains tools for automatically
generating predicted isotope patterns for the molecular formula, and
has attributes for estimated exact mass, nominal mass, error of the
formula, as well as some tools for comparing or plotting patterns.

# Molecule Class

The `Molecule` class was created to allow the user to access molecular
properties of a molecular formula.

## String Specification
**Notes regarding string specification**

- Common abbreviations may be predefined in mass_abbreviations.py
    (either locally or in the current working directory)

- Use brackets to signify multiples of a given component (nested
    brackets are supported)

- Isotopes may be specified using an isotope-element format within a
    bracket (e.g. carbon 13 would be specified as "(13C)" ). The mass
    of that isotope must be defined in the mass dictionary being used
    by the script (default NIST mass).

- The charge may be specified in the formula, but care must be taken
    here. Charge must be specified in either sign-value (e.g. '+2') or
    within a bracket. Otherwise, the script may attempt to interpret the
    charge as a magnitude specifier of the previous block or as an
    isotope, and errors will be encountered.

- A composition dictionary with the format
    `{'Element': number_of_that_element, ...}` may be provided instead
    of a string formula

## Instantiation

A molecule object can be defined by calling the class with a provided
molecular formula.

```
>>> from pythoms.molecule import Molecule
>>> mol = Molecule('C61H51IP3Pd')
>>> mol.monoisotopic_mass  # monoisotopic mass
1109.12832
>>> mol.mw  # molecular weight
1110.300954404405
>>> mol.print_percent_composition() # prints the percent composition
elemental percent composition:
C:  65.99%
H:   4.63%
I:  11.43%
P:  8.369%
Pd:  9.584%
```


## Abbreviations

Abbreviations may be used in the `Molecule` class, but those abbreviations
must be either defined in `pythoms.mass_abbreviations` or in
`mass_abbreviations.py` in the current working directory.
Additional entries can be added by using a unique trigger starting with
an upper case letter and followed by lower case letters. (A new
'element' is triggered by an upper case letter in the molecular formula.)
Each new entry should be a dictionary giving the molecular formula, e.g.
`{'Abbreviation': {'C': 2, 'H': 5}, ...}`

```
>>> from pythoms.mass_abbreviations import abbrvs
>>> abbrvs['Ar+']
{'C': 25, 'H': 21, 'P': 1}
>>> abbrvs['L']
{'C': 18, 'H': 15, 'P': 1}
```

The molecule class may be called using any defined abbreviation.

```
>>> mol = Molecule('L2PdAr+I')
>>> mol.molecular_formula
'C61H51IP3Pd'
```


## Brackets and isotopes

Both brackets and bracket nesting are supported. Supported bracket types
are: `(`, `)`, `[`, `]`, `{`, and `}`.

```
>>> mol = Molecule('Pd{Ar+}I(P(C6H5)3)')
>>> mol.molecular_formula
'C43H36IP2Pd'
```

Isotopes can be specified with the isotope preceeding the element symbol with
both enclosed in a bracket.  e.g.

```
>>> mol1 = Molecule('N(CH2CH3)4')
>>> mol1.molecular_weight
130.2514048325
>>> mol2 = Molecule('N([13C]H2CH3)4')
>>> mol2.molecular_weight
134.2218812385
```

The charge of the molecule can be specified either in brackets in the
molecular formula, or as a keyword argument (see below). The charge in
the molecular formula will override the keyword argument charge. The
charge only has an effect when using the `IPMolecule` class, but is
supported in `Molecule` to prevent issues when providing a charged
formula.

```
>>> mol = Molecule('L2PdAr+(2+)')
>>> mol.molecular_weight
983.3964814044051
```

# IPMolecule Class
A class with many mass-spectrometric properties such as estimated exact
masses, isotope patterns, error estimators, and basic plotting tools.

## Instatiation
```
>>> mol = IPMolecule('C61H51IP3Pd')
```

The isotope pattern molecule objects have a number of useful additional
attributes such as estimated exact mass, patterns, and error.
```
>>> mol.estimated_exact_mass  # estimated exact mass
1109.1303789568635
```

## Stored isotope patterns

There are three isotope patterns stored within each `IPMolecule`
instance: raw, bar, and gaussian. When accessed, these patterns
are returned as a paired list of x and y values (`[[x1, x2, ...], [y1, y2, ...]]`).
The user may then further process or manipulate the lists in any way
they see fit.

```
>>> mz_vals, intensity_vals = mol.bar_isotope_pattern
```

### `raw_isotope_pattern`

The `raw_isotope_pattern` is the calculated isotope pattern with no
combining algorithm applied. Depending on the specified decimal place
tracked (`decpl`) and the complexity of the pattern, there may be many peaks
separated by only very small differences in m/z. Additionally, unless a
dropping method is applied, very, very small intensities are retained.
Keeping these improves the overall accuracy, but can lead to long
computational times for large or isotopically complex molecules. This
is the isotope pattern that would be observed if one had access to a
spectrometer with infinite resolution and dynamic range. The pattern is
stored in case the user needs it, but it is recommended to use the
`bar_isotope_pattern` for visualization purposes.

```
>>> mol.raw_isotope_pattern  # calculated isotope pattern non-consolidated isotope pattern
[[1105.130443, ..., 1225.654769], [3.7321624588364437, 6.425e-320]]
>>> len(mol.raw_isotope_pattern[0])
19344
```

### `bar_isotope_pattern`

The `bar_isotope_pattern` is a pattern generated from the consolidation
of the raw isotope pattern intensities.  A grouping algorithm is
applied to the raw isotope pattern to combine intensities sufficiently
close to one another. The behaviour of this grouping can be tweaked
using the keyword arguments. See the `bar_isotope_pattern` method for
more details.

```
>>> mol.bar_isotope_pattern
[[1105.130443, ..., 1225.6547985235313], [2.285507124030963, ..., 3.9525e-320]]
>>> len(mol.bar_isotope_pattern[0])
121
```

### `gaussian_isotope_pattern`

The `gaussian_isotope_pattern` is a simulated "observed" isotope pattern
generated from the bar isotope pattern. An algorithm is used to
generate normal distributions with a full-width-at-half-maximum
determined from the user-specfied resolution. The normal distributions
are then combined to generate a pattern which will be the most similar
to that observed in a mass spectrometer. Since this process is
computationally arduous, this pattern is not generated on
initialization, but is generated when it is first accessed. This can
extend the accession computation times.

```
>>> mol.gaussian_isotope_pattern
[[1104.68, ..., 1226.11], [0.0, ..., 0.0]]
```

## Calculation methods and Error

The error of the generated isotope patterm may be retrieved via the
`error` attribute. This may be used as a measure of the accuracy of the
calculated isotope pattern, and is calculated by comparing the molecular
weight of the isotope pattern to the expected value. Generally, an
error at or below 10<sup>-6</sup> is deemed an acceptable pattern. The author
has found that while different calculation methods are more
computationally efficient, the error can increase greatly as a result.

```
>>> mol = IPMolecule('C61H51IP3Pd', ipmethod='multiplicative')  # call time 2981 ms
>>> mol.error
-1.2287137530123156e-15
>>> mol2 = IPMolecule('C61H51IP3Pd', ipmethod='combinatorics')  # call time 13120 ms
>>> mol2.error
-1.43349937851437e-15
>>> mol3 = IPMolecule('C61H51IP3Pd', ipmethod='hybrid')  # call time 2022
>>> mol3.error
-5.70102178665877e-6
```

The most visible effects on error come from application of a drop method.
This is where a method is applied to drop or consolidate intensity of
low-intensity peaks. There are several drop methods available to the user,
but user beware of their effect on error. By default, no drop methods
are applied. The dropmethod has no effect on the combinatoric method
(a drop method does not improve computational time in its case).

### `threshold` drop method

The `threshold` dropping method performs a check on every iteration
which drops any peaks below the specified threshold completely from the
spectrum. This can substantially improve computation time especially for
the multiplicative method, but also substantially increases error. This
becomes increasingly apparent for large or complicated patterns, where
initially low intensity values can have substantial effects on the
full pattern.

```
>>> kwargs = {'string': 'C61H51IP3Pd', 'dropmethod': 'threshold'}
>>> mol = IPMolecule(ipmethod='multiplicative', **kwargs)  # call time 22 ms
>>> mol.error
-3.605243759227462e-07
>>> mol3 = IPMolecule(ipmethod='hybrid', **kwargs)  # call time 527 ms
>>> mol3.error
-5.78989879993743e-6
```


### `npeaks` drop method

The `npeaks` dropping method retains the specified number of peaks only,
so the method will drop the `len(spectrum) - npeaks` lowest intensity
peaks. This can have a significant effect on the resulting spectrum.
This method is reportedly used by ChemCalc.

```
>>> kwargs = {'string': 'C61H51IP3Pd', 'dropmethod': 'npeaks'}
>>> mol = IPMolecule(ipmethod='multiplicative', **kwargs)  # call time 2908 ms
>>> mol.error
-1.2287137530123156e-15
>>> mol3 = IPMolecule(ipmethod='hybrid', **kwargs)  # call time 1799 ms
>>> mol3.error
-5.70102178665877e-6
```

### `consolidate` drop method

The `consolidate` drop method is similar to the `threshold` method in
that it checks for peak intensities below a certain threshold, but
consolidates the intensity of the dropped peaks into that of adjacent ones.
A weighted average is used to combine the intensities so as to best
preserve the spectrum accuracy. The `consolidate` keyword argument
is used to control how near adjacent peaks must be to be combined
(the script will check for adjacent peaks 10<sup>-`consolidate`</sup>).
If no adjacent peaks are found, the intensity is dropped from the spectrum.

```
>>> kwargs = {'string': 'C61H51IP3Pd', 'dropmethod': 'consolidate'}
>>> mol = IPMolecule(ipmethod='multiplicative', **kwargs)  # call time 48 ms
>>> mol.error
0.0221406794746035
>>> mol3 = IPMolecule(ipmethod='hybrid', **kwargs)  # call time 737 ms
>>> mol3.error
-5.94680076117617e-6
```

### Recommendations

Should the user encounter situations where the isotope pattern of large
or isotopically-complicated molecules need to be calculated, one will
need to weigh the tradeoff between computational speed or accuracy.
This will largely depend on the resolution of the instrument on which
the experimental spectrum was run, with high-resolution instruments
requiring rigorous pattern calculation, and lower-resolution instruments
requiring only a rough estimate of the pattern. The user will also want
to consider whether they need to evaluate the pattern qualitatively or
more precisely.

## Plotting

There are basic plotting methods provided within the `IPMolecule` class
for all three patterns so the user will be able to visually evaluate the
pattern that was generated. Should the user wish to create an isotope
pattern overlay figure, the `isotope pattern overlay.py` file contains a
script for this purpose.
