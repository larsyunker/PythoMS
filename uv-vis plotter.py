"""
UV-Vis plotting tool 

Extracts UV-Vis spectra from an mzML file and plots traces from start to end with spacing deltat
overrides can be specified to adjust the output

to do:
    
"""

# specify mzml file to extract from
filename = 'BTM-42 UVVis'

# time range to plot
# the plotted range will be [start,end) if the total time/deltat is a whole number
tr = [10,90.1]

# time spacing between the traces
deltat = 10

# override settings here
override = {
#'fs':16, # font size
#'lw':1.5, # line width of traces
#'size':[7.87,4.87], # image size [width,length] in inches
#'xrange':[500,700], # wavelength bounds (in nm)
#'yrange':[0,3], # absorbance bounds(in a.u.)
#'legloc':0, # legend location (see ttp://matplotlib.org/api/legend_api.html for more location codes)
}

if __name__ == '__main__':
    from _classes._mzML import mzML
    from scipy import arange
    from tome_v02 import takeclosest,plotuv
    mzml = mzML(filename) # initiate mzml object
    uvspecs,wavelengths = mzml.pulluvspectra(tr=tr) # pull uv spectra between the specified points
    timepoints = sorted(uvspecs.keys()) # extract timepoints
    times = arange(tr[0],tr[1],deltat) # evenly spaced times between start and end
    
    specin = []
    for time in times:
        ind = takeclosest(timepoints,time) # find the closest time to that
        specin.append(uvspecs[timepoints[ind]]['intensity']) # append that spectrum to the input list
    
    if override.has_key('outname') is False:
        override['outname'] = filename.split('.')[0] + ' UV-Vis'
    
    plotuv(wavelengths,specin,times=times,**override) # plot and save