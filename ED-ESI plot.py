"""
for use with EDESI pulling

new:
    I hope this works for Darien
    ---0.1 beta---

"""
import sys
from _classes._mzML import mzML
from tome_v02 import binnspectra,bincidspectra
import matplotlib
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt

filename = 'HZ-140516_HOTKEYMSMS 1376 II.raw' # raw or mzml file name
fillzeros = True # fills spectrum with zeros
decpl = 1 # number of decimal places to track
threshold = 1 # threshold intensity value
mzrange=None # mzrange to track
sr = 'all' # scan range to track
#
mzml = mzML(filename,verbose=True)

# pull time list, each spectrum at those time points, the scan range, and the mz range
#speclist,sr,mzrange = mzml.pullspectra(mzrange=mzrange,sr=sr)
msms,limits = mzml.pullmsmsspectra()
binned = {}
grouped = {}
voltage = {}
zvals = {} # intensity matrix

for ion in msms:
    binned[ion] = binnspectra(msms[ion], len(msms[ion]), dec=decpl, startmz=limits[ion][0], endmz=limits[ion][1])[0] # bins all spectra together
    grouped[ion],voltage[ion] = bincidspectra(msms[ion], dec=decpl, startmz=limits[ion][0], endmz=limits[ion][1], threshold=threshold, fillzeros=fillzeros)
    for ind,scan in enumerate(grouped[ion]):
        sys.stdout.write('\rExtracting z values #%d/%d (%.1f%%) ' %(ind+1,len(grouped),float(ind+1)/float(len(grouped))*100.))
        if zvals.has_key(ion) is False:
            zvals[ion] = []
        zvals[ion].append(scan[1]) # append only the intensities to list
    sys.stdout.write(' DONE\n')

"""
each key in msms, binned, grouped, and voltage is a single ion that is being tracked
subkeys in msms are timepoints, containing subkeys for collision energy, x, and y
subkeys in binned is a [[m/z],[intensities]] structure which is all spectra binned
subkeys in grouped are voltages, which each have a [[m/z],[intensities]] spectrum
voltage[ion] is a list of voltages

axes for ED ESI plots
for each ion:
    xvalues = grouped[ion][0][0]
    yvalues = voltage[ion]
    z values = zvals[ion]
"""
ContourPlotLV = ""
for ion in msms:
    print ion
    ContourPlotLV = plt.contour(grouped[ion][0][0], voltage[ion], zvals[ion])

#plt.plot(binned[u'1376.30005'][0], binned[u'1376.30005'][1])