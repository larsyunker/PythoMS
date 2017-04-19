'''what is the molecular formula?'''
#formula = 'PdL2PhNC5H4CH3OTf'
#formula = 'PhBO3C5H9NEt4'
#formula = 'FC6H4BF3K'
#formula = 'Cs2CO3'
formula = 'Na2CO3'
#formula = 'Ar+IPF6'
#formula = 'TolB(OH)2'
#formula = 'H2O'

'''how much?'''
amount = 30

'''what is the unit?'''
unit = [
#'g',
#'mg',
#'ug',
#'mol',
#'mmol',
'umol',
]

from _classes._Molecule import Molecule

mol = Molecule(formula,
dropmethod='threshold',
)

if len(unit) > 1:
    raise ValueError('only one unit may be specified')
if unit[0].startswith('m'): # milli
    if unit[0][1] != 'o':
        order = 3
    else: # if 'mol'
        order = 0
elif unit[0].startswith('u'): # micro
    order = 6
else: # standard
    order = 0

amount = float(amount) # make sure it's a float
if unit[0].endswith('g'): # mass
    nmol = (amount/10**order)/mol.mw
    print mol, mol.mw
    print 'mass -> mol'
    print '%.2f mol' %(nmol)
    print '%.2f mmol' %(nmol*1000)
    print '%.2f umol' %(nmol*10**6)
elif unit[0].endswith('l'): #moles
    mass = (amount/10**order)*mol.mw
    print mol, mol.mw
    print 'mol -> mass'
    print '%.2f g' %mass
    print '%.2f mg' %(mass*1000)
    #print '%.2f ug' %(mass*10**6)