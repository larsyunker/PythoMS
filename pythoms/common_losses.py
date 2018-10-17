"""
A common loss dictionary for use with msms interpreter assistant
values are stored with 5 decimal places, but are returned with n decimal places
"""

stored_dec = 5  # number of decimal places stored in the following dictionary
losses = {  # dictionary of common losses and their probable assignments
    14.003074: 'N',
    14.01565: 'CH2',
    15.02347: 'CH3',
    15.99491: 'O',
    17.00274: 'OH',
    18.01056: 'H2O',
    22.98977: 'Na',
    28.03130: 'C2H4',
    31.01839: 'MeO',
    31.97207: 'S',
    31.98983: 'O2',
    32.02621: 'MeOH',
    34.96885: 'Cl',
    38.96371: 'K',
    41.02655: 'MeCN',
    65.03912: 'Cp (cyclopentadienyl)',
    72.05751: 'THF',
    77.03912: 'Ph',
    78.04695: 'Benzene',
    78.91834: 'Br',
    83.95336: 'CH2Cl2',
    91.05478: 'C7H7 (CH2Ph, tropylium)',
    96.03753: 'PhF (fluorobenzene)',
    105.90349: 'Pd',
    122.08440: '4-methylpyridine',
    126.90447: 'I',
    135.11738: 'Cp* (pentamethyl cyclopentadienyl)',
    262.09114: 'PPh3',
    352.13809: 'Ar+ (C25H21P)',
}
