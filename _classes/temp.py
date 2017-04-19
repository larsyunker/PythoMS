import operator as op
from _ScriptTime import ScriptTime
import sympy as sym

st = ScriptTime(profile=True)

#@st.profilefn
#def ncr(n, r):
#    r = min(r, n-r)
#    if r == 0: return 1
#    numer = reduce(op.mul, xrange(n, n-r, -1))
#    denom = reduce(op.mul, xrange(1, r+1))
#    return numer//denom
#
#def product(iterable):
#    prod = 1
#    for n in iterable:
#        prod *= n
#    return prod
#
#def npr(n, r):
#    """
#    Calculate the number of ordered permutations of r items taken from a
#    population of size n.
#
#    >>> npr(3, 2)
#    6
#    >>> npr(100, 20)
#    1303995018204712451095685346159820800000
#    """
#    assert 0 <= r <= n
#    return product(range(n - r + 1, n + 1))
#
#def ncr(n, r):
#    """
#    Calculate the number of unordered combinations of r items taken from a
#    population of size n.
#
#    >>> ncr(3, 2)
#    3
#    >>> ncr(100, 20)
#    535983370403809682970
#    >>> ncr(100000, 1000) == ncr(100000, 99000)
#    True
#    """
#    assert 0 <= r <= n
#    if r > n // 2:
#        r = n - r
#    return npr(n, r) // factorial(r)

#@st.profilefn
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
#@st.profilefn
#def symchoose(n,k):
#    fn = sym.factorial(n)
#    fn /= sym.factorial(k)
#    fn /= sym.factorial(n-k)
#    return fn.evalf()

@st.profilefn
def cwr(n,k):
    """
    calculates the number of combinations with repitition
    n: number of things to choose from
    k: choose k of them
    """
    fn = sym.factorial(n+k-1)
    fn /= sym.factorial(k)
    fn /= sym.factorial(n-1)
    return fn.evalf()

from scipy.special import factorial
@st.profilefn
def cwrnp(n,k):
    return factorial(n+k-1)/(factorial(k)*factorial(n-1))

n = 2
k = 3900
print cwr(n,k)
print cwr(n,k)
for i in range(1000):
    cwr(n,k)
    cwrnp(n,k)

st.printprofiles()