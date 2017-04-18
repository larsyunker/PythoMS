from _classes._Molecule import Molecule
from _classes._nist_mass import nist_mass as mass

#formula = 'C2H7N'
#mol = Molecule(formula)
#comp = dict(mol.comp)

comp = {'H': 7, 'C': 2, 'N': 1}

combinations = {}
for ele in comp:
    