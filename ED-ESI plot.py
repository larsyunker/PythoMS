import sys
import matplotlib
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
from scipy import *

from _classes._mzML import mzML
from _classes._Spectrum import Spectrum
from tome_v02 import binnspectra,bincidspectra

"""
for use with EDESI pulling

new:
    I hope this works for Darien
    ---0.1 beta---
    Implemented PyPlot EDESI Mapping using the _mzML framework
    Added Helper Functions to the EDESI Plot Production
    Cleaned up the implementation of EDESIPlot
    -> Requires the input of SummedSpectrum, ContourPlot, BreakdownPlot
       information
    -> Input currently is designed for by calls by _mzML framework
       as the SummedSpectrum input is dependent on the structure of 
       binned[ions]
    To Change -> Implement a generalized version of EDESIPlot for
                 use cases outside the _mzML framework
"""
#######DEBUG FUNCTION#########

def debugPrint(isDebug, printStr):
    if(isDebug):
       print printStr  

def anon(x):
   #Used for testing lambda functions
   print x

######HELPER FUNCTION#########
def GraphCleanUp(ax):
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()

def RoundToTen(val):
    return int(ceil(val/10.0))*10

def RoundToTenFloor(val):
    return int(floor(val/10.0))*10

def RoundToHundred(val):
    return int(ceil(val/100.0))*100

def round1sig(yTick):
    for i in range(0, len(yTick)):
        #stackoverflow: 3410976
        if (yTick[i] == 0):
            continue
        yTick[i] = round(yTick[i], -int(floor(log10(abs(yTick[i])))))
    return yTick

def Valueround1sig(yTick):
    return round(yTick, -int(floor(log10(abs(yTick)))))

def normalizeMSSpectrum(mzList, maxValue):
    relInt = []
    for i in mzList:
       relInt.append((float(i)/float(maxValue))*100)
       #print ((float(i)/float(maxValue))*100)
    return relInt

#####################################

####EDESI Plot Function##############
#########################################################################
# PlotEDESI 
# ------------------------------
# SummedSpectrum -> pass summed
#   -> sumSpecX, sumSpecY = binned[ion]
#
# ContourPlot -> grouped[ion][0][0], voltage[ion], zvals[ion]
#   -> List of lists, ContourPlot[0]: mz, [1]: voltage, [2]: Intensity
#
# BreakdownCurve -> mzml.pullspeciesdata(sp), voltage[ion]
#   -> traceBank
#   -> for i in BreakdownCurve[0]:
#         plot(i['raw'], voltage[ion])
#
# No Return Object -> Saves Figure to Current Working Directory (CWD)
#########################################################################
def PlotEDESI(SummedSpectrum, ContourPlot, BreakdownCurve=False, **kwargs):
    
    #EDESI Plot Production variable
    outputFile = './DYLYTest'
    EDESIkwargs = {
        'minFilter': 20, # minFilter intensity value
        'threshold': 1100, # threshold of peak height for Breakdown tracing
        'plotZoom': True # Construct Plot with Zoom in region of interest? (Autozoom)
        'debug' = True
    }
    if set(kwargs.keys()) - set(EDESIkwargs.keys()):
        for i in set(kwargs.keys()) - set(EDESIkwargs.keys()):
            raise KeyError('Unsupported keyword argument: {errKWarg}'.format(errKWarg = str(i)))
    EDESIkwargs.update(kwargs)


    line_colours = ('BlueViolet', 'Crimson', 'ForestGreen', 'Indigo', 'Tomato', 'Maroon')
    levels = arange(minFilter, max(array(ContourPlot[2]).max(axis=1)), 6)

    subplotVal = 1 #Contour Plot position
    if(BreakdownCurve):
        subplotVal = 2

    plt.figure(figsize=(8,9)) #Figure Size 8'', 9''
    matplotlib.rcParams.update({'font.size':16})
    matplotlib.rcParams['xtick.direction'] = 'out'
    matplotlib.rcParams['ytick.direction'] = 'out'

    gs1 = gridspec.GridSpec(subplotVal,2,width_ratios=[3,1], height_ratios=[1,3]) #Optimal Dimension
    gs1.update(wspace=0.025, hspace=0.025)

    ## Contour Plot -> ax2 ##
    debugPrint(debug, "Contour Start")
    
    ax2 = plt.subplot(gs1[subplotVal])
    
    ax2.contour(ContourPlot[0], ContourPlot[1], ContourPlot[2], levels, colors = 'k')
    ContouryTick = arange(0, Valueround1sig(max(ContourPlot[1]))+10, 10)
    #Labels the smallest value of the axis and then continue labelling with the arange increment ie 50, 750, 1500, 2250 rather than 50, 800, 1550, 2300 (tl;dr 50 does not affect the arange values)
    ax2.set_yticks(concatenate(([Valueround1sig(min(ContourPlot[1]))],ContouryTick[:-1]), axis=0))
    #ContourxTick = arange(0, Valueround1sig(max(binArr.keys()))+750, 750) #Unused
    ax2.set_ylabel('Collision Energy (V)')
    ax2.set_xlabel('m/z', style='italic')
    GraphCleanUp(ax2)
    
    debugPrint(debug, "Contour End")

    ## Summed Spectrum -> ax1 ##
    debugPrint(debug, "Sum Spec Start")
    
    ax1 = plt.subplot(gs1[0], sharex=ax2)
    sumSpecX, sumSpecY = SummedSpectrum
    ax1.plot(sumSpecX, sumSpecY)
    plt.title('Energy Dependent-ESI MS/MS Plot M={MINFIL}'.format(MINFIL=minFilter), y=1.08)
    plt.setp(ax1.get_xticklabels(), visible=False)
    ax1.set_ylabel('Intensity (%)')
    #stackoverflow: 11244514
    #ax1.set_yticks([0, 50, 100]) #Relative Intensity`
    GraphCleanUp(ax1)
    ax1.spines['bottom'].set_visible(False)
    ax1.xaxis.set_ticks_position('none')
    
    debugPrint(debug, "Sum Spec End")

    ## Zooming into particular m/z ##
    if (plotZoom):
        StartRange = min(map(float, sumSpecX))
        EndRange = max(map(float, sumSpecX))
        Range = (EndRange-StartRange)/10
        xMinVal = StartRange - Range
        xMaxVal = EndRange + Range
        xtickRange = arange(RoundToHundred(StartRange)-100, RoundToHundred(EndRange)+100, 100)
        ax2.set_xticks(xtickRange)
        ax2.set_xlim(xmin=xMinVal, xmax=xMaxVal)
        ax1.set_xlim(xmin=xMinVal, xmax=xMaxVal)

    ## Constructing Breakdown Curve -> ax3 ##
    traceBank = {}
    print BreakdownCurve
    if(BreakdownCurve):
        debugPrint(debug, "Breakdown Start")
        
        for mz in xrange(len(SummedSpectrum[0])):
            if SummedSpectrum[1][mz] > threshold:
               subDict = {'bounds': [int(SummedSpectrum[0][mz])], 'affin': '-', 'raw': [], 'spectrum': Spectrum(1)} 
               traceBank[int(SummedSpectrum[0][mz])] = subDict
        outTrace, TIC, rtime = mzml.pullspeciesdata(traceBank)
        ax3 = plt.subplot(gs1[3], sharey=ax2)
        for trace in traceBank:
            debugPrint(debug, outTrace)
            ax3.plot(outTrace[trace]['raw'], ContourPlot[1])
        plt.setp(ax3.get_yticklabels(), visible=False)
        ax3.set_xlabel('Intensity (%)')
        ax3.set_xticks([100])
        GraphCleanUp(ax3)
        
        debugPrint(debug, "Breakdown End")

    plt.savefig('../{OUTFILE}.png'.format(OUTFILE = outputFile+str(minFilter)), bbox_inches='tight')
##################################

###############################################################
#MAIN
###############################################################
#_mzML processing variables
filename = 'HZ-140516_HOTKEYMSMS 1376 II.raw' # raw or mzml file name
fillzeros = True # fills spectrum with zeros
decpl = 1 # number of decimal places to track
mzrange = None # mzrange to track
sr = 'all' # scan range to track
mzml = mzML(filename,verbose=True)

#EDESI Plot Production variable
minFilter = 20 # minFilter intensity value
threshold = 1156 # threshold of peak height for Breakdown tracing
plotBreakdown = True # Construct Plot with Breakdown?
plotZoom = True # Construct Plot with Zoom in region of interest? (Autozoom)

# pull time list, each spectrum at those time points, the scan range, and the mz range
#speclist,sr,mzrange = mzml.pullspectra(mzrange=mzrange,sr=sr)
msms,limits = mzml.pullmsmsspectra()
binned = {}
grouped = {}
voltage = {}
zvals = {} # intensity matrix

for ion in msms:
    binned[ion] = binnspectra(msms[ion], len(msms[ion]), dec=decpl, startmz=limits[ion][0], endmz=limits[ion][1])[0] # bins all spectra together
    grouped[ion],voltage[ion] = bincidspectra(msms[ion], dec=decpl, startmz=limits[ion][0], endmz=limits[ion][1], threshold=minFilter, fillzeros=fillzeros)
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
#ContourPlotLV = ""
for ion in msms:
    PlotEDESI(binned[ion], 
          [grouped[ion][0][0], voltage[ion], zvals[ion]], 
           )
    #print ion
    #levels = arange(minFilter, max(array(zvals[ion]).max(axis=1)), 6)
    #ContourPlotLV = plt.contour(grouped[ion][0][0], sorted(voltage[ion]), zvals[ion], levels, color='k')


