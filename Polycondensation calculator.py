from pythoms.classes import Molecule

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
    print('%d\t%.2f' % (i, Molecule(Ar + (arunit * i) + X, **kwargs).em))

print('\nPd units')
print('n\tmass')
for i in range(1, n):
    print('%d\t%.2f' % (i, Molecule('L2Pd' + Ar + (arunit * i) + X, **kwargs).em))

print('\ncapped')
print('n\tmass')
for i in range(0, n + 3):
    print('%d\t%.2f' % (i, Molecule(cap + Ar + (arunit * i), **kwargs).em))
