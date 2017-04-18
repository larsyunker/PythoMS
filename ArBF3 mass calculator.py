from _classes._Molecule import Molecule

aryl = 'Pr'

print aryl+'BF3',Molecule(aryl+'BF3').em
for i in range(2,5):
    string = '('+aryl+'BF3)'+str(i)+'K'+str(i-1)
    mol = Molecule(string)
    print string, mol.em
print '('+aryl+'BF3)2'+'Cs',Molecule('('+aryl+'BF3)2'+'Cs').em
print '\n'
print aryl+'BF2OH',Molecule(aryl+'BF2OH').em
print aryl+'BF(OH)2',Molecule(aryl+'BF(OH)2').em
print aryl+'B(OH)3',Molecule(aryl+'B(OH)3').em
print '('+aryl+'B)2O3H',Molecule('('+aryl+'B)2O3H').em
print aryl+'BO2H',Molecule(aryl+'BO2H').em