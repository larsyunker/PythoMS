from PythoMS._classes._nist_mass import nist_mass

print('nist_mass = {')
for key in sorted(nist_mass.keys()):
    print("    '%s': {" % key)
    element = nist_mass[key]
    for subkey in sorted(element.keys()):
        if element[subkey][1] == 0.:
            print("        %d: (%f, %.1f), " % (subkey, element[subkey][0], element[subkey][1]))
            continue
        print("        %d: (%f, %f), " % (subkey, element[subkey][0], element[subkey][1]))
    print('    }, ')
print('} \n')

class Test(object):
    def __init__(self):
        self.lst = [1,2,3]
        self.genlst = self.gen

    def gen(self):
        for val in self.lst:
            yield val

    def __iter__(self):
        for i in self.lst:
            yield i


t = Test()
print(1 in t)