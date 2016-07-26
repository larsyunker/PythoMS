"""
UV-Vis plotting tool 

Extracts UV-Vis spectra from an mzML file and plots traces from start to end with spacing deltat
overrides can be specified to adjust the output

new:
    updated to work with mzML 2.4
    ---1.1---

to do:
    
"""

# specify mzml file to extract from
filename = 'BTM-42 UVVis'

# time range to plot
# start point (either scan number or time point [time must be a float])
start = 10.
# end point
end = 90.1

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
    from tome_v02 import locateinlist,plotuv
    mzml = mzML(filename) # initiate mzml object
    fn = mzml.associate_to_function('UV') # determine which function contains UV-Vis data
    uvspecs = mzml.retrieve_scans(start,end,fn) # pull uv spectra
    wavelengths = list(uvspecs[0][0]) # wavelength list
    uvspecs = [y for x,y in uvspecs] # set uvspecs list to be only the y values
    timepoints = mzml.functions[fn]['timepoints'] # pull time points of the UV function
    l,r = locateinlist(timepoints,start,'greater'),locateinlist(timepoints,end,'lesser') # locate indicies of timepoints
    timepoints = timepoints[l:r+1] # trim time list accordingly
    times = arange(start,end,deltat) # evenly spaced times between start and end
    
    specin = []
    for time in times:
        ind = locateinlist(timepoints,time) # find the closest time to that
        specin.append(uvspecs[ind]) # append that spectrum to the input list
    
    if override.has_key('outname') is False:
        override['outname'] = filename.split('.')[0] + ' UV-Vis'
    
    plotuv(wavelengths,specin,times=times,**override) # plot and save