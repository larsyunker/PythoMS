from _classes._XLSX import XLSX
from _classes._Molecule import Molecule
import time as t
import sys
wr = sys.stdout.write


ilist = range(0,101,2)
ilist[0] = 1
comb = []
mul = []
mw = []

base = 'C25H21PI' # Ar+I
base = 'Cys' # cystine

for i in ilist:
    target = '('+base+')'+str(i)
    wr('Multiplicative %d' %i)
    t0 = t.time()
    mol = Molecule(target,ipmethod='multiplicative',
        #verbose=True
        )
    mul.append(t.time() - t0)
    wr(' DONE\n')
    
    wr('Combinatorics  %d' %i)
    t0 = t.time()
    mol = Molecule(target,ipmethod='combinatorics',
        #verbose=True
        )
    comb.append(t.time() - t0)
    mw.append(mol.mw)
    
    wr(' DONE\n')
    
    if i > 10:
        break
    