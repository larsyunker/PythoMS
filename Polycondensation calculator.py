from pythoms.molecule import IPMolecule

kwargs = {
    'dropmethod': 'threshold'
}

Ar = 'Ar+'
X = 'I'
arunit = 'C6H4'
cap = 'Ph'

n = 6

print(Ar)
print('aromatic units')
print('n\tmass')
for i in range(1, n):
    print('%d\t%.2f' % (i, IPMolecule(Ar + (arunit * i) + X, **kwargs).estimated_exact_mass))

print('\nPd units')
print('n\tmass')
for i in range(1, n):
    print('%d\t%.2f' % (i, IPMolecule('L2Pd' + Ar + (arunit * i) + X, **kwargs).estimated_exact_mass))

print('\ncapped')
print('n\tmass')
for i in range(0, n + 3):
    print('%d\t%.2f' % (i, IPMolecule(cap + Ar + (arunit * i), **kwargs).estimated_exact_mass))
