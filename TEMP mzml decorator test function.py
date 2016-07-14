from _classes._mzML import mzML

mzml = mzML('HZ-140516_HOTKEYMSMS 1376 II')

@mzml._foreachscan
def dothis(spectrum):
    out = []
    p = mzml.cvparam(spectrum)
    out.append(float(p['MS:1000016']['value'])) # scan start time
    out.append(int(p['MS:1000505']['value'])) # base peak intensity
    return out

doit = dothis()

x = []
y = []
for val in doit:
    x.append(val[0])
    y.append(val[1])