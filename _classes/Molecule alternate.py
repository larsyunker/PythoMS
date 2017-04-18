#def choose(n, k):
#    """
#    A fast way to calculate binomial coefficients by Andrew Dalke (contrib).
#    """
#    if 0 <= k <= n:
#        ntok = 1
#        ktok = 1
#        for t in xrange(1, min(k, n - k) + 1):
#            ntok *= n
#            ktok *= t
#            n -= 1
#        return ntok // ktok
#    else:
#        return 0
#
#print choose(4,0)


#import itertools

#test = []
#for item in itertools.combinations_with_replacement(lst,4):
#    test.append(item)

from _nist_mass import nist_mass as mass
from itertools import combinations_with_replacement as cwr
#from itertools import product as prod
import numpy as np
from scipy.special import factorial

comp = {'S':1,'O':4}

#def elecomb(ele,num):
def num_permu(lst,isos):
    """calculates the number of unique permutations of the given set of isotopes"""
    counts = [lst.count(x) for x in isos] # counts the number of each isotope in the set
    fact = np.prod([factorial(x) for x in counts]) # calculates the product of the factorials of those counts
    return int(factorial(len(lst))/fact) # return the factorial of the list length over the count factorial product
#
#    
#    isos = [] # isotopes list
#    #sums = {}
#    for isotope in mass[ele]: # for each isotope of that element
#        if isotope != 0 and mass[ele][isotope][1] != 0: # if intensity is nonzero, append isotope number
#            isos.append(isotope)
#    for comb in cwr(isos,num):
#        
#        if tot not in sums:
#            sums[tot] = [comb,1]
#        else:
#            sums[tot][1] += 1
#    out = {}
#    for item in sums:
#        out[sums[item][0]] = sums[item][1]
#    print out
#    for item in out:
#        #count = item.count(isos[2])
#        counts = [item.count(x) for x in isos]
#        print item,out[item],[choose(len(item),x) for x in counts]

#dct = {
#'S':[(32, 'S'), (33, 'S'), (34, 'S'), (36, 'S')],
#'O':[(16, 'O'), (17, 'O'), (18, 'O')]
#}

from _nist_mass import nist_mass as mass
from itertools import combinations_with_replacement as cwr

class mfiter(object):
    def __init__(self,comp):
        self.dct = {}
        for element in comp:
            self.dct[element] = []
            for isotope in mass[element]:
                if isotope != 0 and mass[element][isotope][1] != 0:
                    self.dct[element].append(isotope)
        
        
    def __iter__(self): return self
    
    def next(self):
        if self.start >= self.stop:
            raise StopIteration
        current = self.start * self.start
        self.start += 1
        return current

a = mfiter(comp)
#for i in a:
#    print i


isos = [] # isotopes list
for element in comp:
    for isotope in mass[element]:
        if isotope != 0 and mass[element][isotope][1] != 0:
            isos.append((isotope,element))
print isos

"""
generate an iterator for each element
generate the cwr combinations list
calculate the number of possible permutations from the cwr combinations using counts (possible?)

generate a meta-generator to iterate over every combination of every element
as each iteration proceeds, add the m/z and intensity vlaues to a spectrum object

# build an iterator that will combine each combination for an element with each other combination for another element

# faster to combine all elements and combinations and work outward, or to do each element individually and combine
"""

#for ele in comp:
#    elecomb(ele,comp[ele])