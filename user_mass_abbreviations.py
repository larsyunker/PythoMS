"""
User-specified mass abbreviations for the `Molecule` class may be specified here. These will be added to the default
mass abbreviations in `pythoms.mass_abbreviations` .

Dictionary of common abbreviations used in chemical formulas
A new 'element' is triggered by an upper case letter, so these abbreviations must begin with an upper case letter
followed by lower case letters, and the entire combination must not be the same as an elemental signifier.

Specified abbreviations should take the form of
'Abbrv': {'element': integer_number_of_atoms, ...}
"""

user_abbrvs = {
    'Abbreviation': {'C': 1, 'H': 1},  # example abbreviation
}
