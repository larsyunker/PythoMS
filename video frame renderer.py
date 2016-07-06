"""
Grabs data from a mzML file and renders a series two-subplot images for
combination into a video sequence. 
The left subplot is a mass spectrum at the current time
The right subplot is a chart of intensity traces leading up to the current time

The script assumes the video will be rendered in 1080p @ 30fps

new:
    works with _mzML class now
    should be able to plot both - and + species
    normalizes to the summed intensity of the specified species
    no longer requires an excel file for parameters (set species information in script)
    incorporated species trace and spectrum binning
    switched spectrum binning to use NoneSpectrum class (avoids improper copying)
    ---9.0---
    renamed from "convert spec chrom png" to "video fram renderer"
    updated calls to functions
    made all text based on fontsize
    changed timepoint colour to blue (should be more visible)
    tweaked the spacing between the two plots
    ---9.1---
    ---9.2

to add:
    add functionality to track both + and - mode spectra (will have to modify pullspectra() in _mzML)
"""

# input *.raw filename
filename = '20160325EXPT2.raw'

# species to track and plot
# sub and superscripts can be denoted by TeX formatting
# e.g. an Ar group with a positive charge is denoted by 'Ar$^+$', and CH2 would be denoted CH$_2$ (see http://matplotlib.org/users/mathtext.html#subscripts-and-superscripts for more info)
# m/z bounds, the desired colour, and affinity must be specified for every species
sp = {
' ':{'bounds':[825.223,837.0],'colour':(0,128,0),'affin':'+'}, #[Ru](PPh$_3$)$_2$(PEt$_2$H)
'  ':{'bounds':[838.277,851.0],'colour':(255,0,0),'affin':'+'}, # [Ru](PPh$_3$)$_2$(NCPh)
'   ':{'bounds':[921.34,934.0],'colour':(0,0,255),'affin':'+'} # [Ru](PPh$_3$)$_2$(PPh$_2$H)
}

# set number of scans to sum
n = 1

# define scan range
# (scan range that you wish to render)
scr = [208,878]
#scr = [261,427]

# define mz range for spectrum image
# (these will be the bounds of the x-axis)
mz = [820.,939.]

# left-right scalar for scan info placement
# 0 is fuly right, 1 is fully left
infop = 0.1

# reaction start point (eg. catalyst injection)
# provide time in minutes
inj = 5.074

# provide timepoints for injections/additions
# give the form [['injection name',time in min],etc...]
timepoints = {
'phosphines added':5.074
}

# axis line width
axwidth = 1.0

# fontsize
fs = 13

# save an image every n scans 
save = 1

def spectrumtrace(filename,sp,scr='all',n=1,mz='all',inj=0.,save=1):
    """
    plots a mass spectrum and an intensity trace for every scan in a raw file based on the parameters supplied
    
    filename: the *.raw filename in the working directory to work from
    sp: dictionary of species to track and render
    scr: scan range to sum
        [start,end]
    n: number of scans to sum
        integer
    mz: mz range to track
        [m/z start, m/z end]
    inj: injection point (e.g. for catalyst injection)
        float
    save: save every # number of scans
        integer
    """
    import os,sys
    sys.path.append(os.path.dirname(os.path.realpath(__file__))+'/_classes')
    from tome_v02 import bindata,binnspectra
    from _mzML import mzML
    from _ScriptTime import ScriptTime
    from _Spectrum import Spectrum
    from _Colour import Colour
    from bisect import bisect_left ,bisect_right
    import pylab as pl
    import os,sys
    
    st = ScriptTime(profile=True)
    
    @st.profilefn
    def plotit(x,y,index):
        """
        generates plot
        input: x and y of mass spectrum, index of current time point
        """
        fig = pl.figure(figsize = (9.6,5.4),dpi= 100) # set figure size to 1920x1080
        font = {'fontname':'Arial'} #font parameters for axis/text labels
        tickfont = pl.matplotlib.font_manager.FontProperties(family='Arial',size=fs) # font parameters for axis ticks
    
        axl = fig.add_subplot(121) # left subplot (mass spectrum)
        axl.spines["right"].set_visible(False)
        axl.spines["top"].set_visible(False)
        axl.spines["bottom"].set_visible(False)
        axl.plot(x,y, 'k-', lw=0.75)
        
        pl.xlabel('m/z', style='italic',**font)
        pl.ylabel('Relative Intensity',**font)
        for axis in ["top","bottom","left","right"]:
            axl.spines[axis].set_linewidth(axwidth)
        axl.spines["bottom"].set_position(('axes',-0.01)) #offset x axis
        for label in axl.get_yticklabels():
            label.set_fontproperties(tickfont)
        for label in axl.get_xticklabels():
            label.set_fontproperties(tickfont)    
        pl.tick_params(axis='y', length=axwidth*3, width=axwidth, direction='out',right='off')
        pl.tick_params(axis='x', length=axwidth*3, width=axwidth, direction='out',top='off')
        for key in sp: #ind,val in enumerate(sp):
            l,r = bisect_left(x,sp[key]['bounds'][0]),bisect_right(x,sp[key]['bounds'][1]) # index location of selected peak in spectrum
            axl.plot(x[l:r],y[l:r],color = sp[key]['colour'], lw=1) # plot spectrum in colour for selected peaks
            axl.text(sp[key]['bounds'][0],1.01,key,fontsize=fs,color=sp[key]['colour'])
            
        pl.xlim(mz)
        pl.ylim([-0.001,1])
        
        axr = fig.add_subplot(122) # right subplot (chromatogram)
        pl.xlim([mintime,maxtime])
        pl.ylim([-0.001,1])
        axr.spines["right"].set_visible(False)
        axr.spines["top"].set_visible(False)
        pl.tick_params(axis='y', length=axwidth*3, width=axwidth, direction='out',right='off')
        pl.tick_params(axis='x', length=axwidth*3, width=axwidth, direction='out',top='off')
        
        for mode in mskeys: 
            sumkey = str(n)+'sum'+mode
            spkey = str(n)+'norm'
            for key in sp:
                if sp[key]['affin'] is mode: # pair species with appropriate rtime
                    axr.plot(rtime[sumkey][:index],sp[key][spkey][:index], linewidth = 1.0, label = key, color = sp[key]['colour'])
        pl.xlabel('time (min)',fontsize = fs, **font)
        pl.tick_params(axis='y',labelleft='off')
        for label in axr.get_yticklabels():
            label.set_fontproperties(tickfont)
        for label in axr.get_xticklabels():
            label.set_fontproperties(tickfont)
        
        for key in timepoints: # add vertical timepoint lines
            if maxtime+inj >= timepoints[key]: # and rtime[0] <= timepoints[ind][1]
                pl.axvline(x=(timepoints[key]-inj), ymin = 0, ymax = 1, linewidth=0.75, color = 'b', linestyle = ':')
                pl.text(timepoints[key]-inj,0.5,key, fontsize = fs, color = 'b', backgroundcolor = 'w', rotation = 'vertical', horizontalalignment='center',verticalalignment='center',alpha = 0.75,**font)
        textx = maxtime - (maxtime-mintime)*infop # calculate location for scan number and time text
        pl.text(textx,0.96, 'scan %i' %curspec,fontsize=fs,**font) # text for scan number
        pl.text(textx,0.92, '%.1f min' %maxtime,fontsize=fs,**font) # text for time
        pl.subplots_adjust(left = 0.07, right = 0.99, bottom = 0.095, top = 0.96, wspace = 0.06, hspace = 0.05) # hard coded subplot tightening
        #pl.tight_layout(pad=0.75) # automatically tighten subplots
        dpiset = 200 # 100 is 960x540, 150dpi is 1440x810, 200dpi is 1920x1080
        pl.savefig(os.getcwd()+r'\imgs\scan'+str(itr)[2:6]+'.png',figsize = (19.2,10.8),dpi=dpiset)
        pl.clf()
        pl.close()
    
    @st.profilefn
    def msfignorm(x,y):
        """
        Normalizes the height of a mass spectrum
        The height will be the sum of the heights of the base peaks in the window
        
        The function will normalize the y-values (assumes intensity) and return them
        """
        height = 0 # starting point
        for key in sp:
            l,r = bisect_right(x,sp[key]['bounds'][0]),bisect_left(x,sp[key]['bounds'][1]) # index location of selected peak in spectrum
            try:
                height += max(y[l:r]) # add maximum in selected region to height
            except ValueError: # if no intensity in region
                height += 0.01    
        
        for ind,val in enumerate(y): #normalizes all y values
            y[ind] = val/height
        
        return y
    
    def timelimits(index):
        """
        finds the appropriate time limits for the traces
        """
        mintime = 10000
        maxtime = -10000
        for mode in mskeys:
            sumkey = str(n)+'sum'+mode
            if sumkey in rtime.keys():
                if rtime[sumkey][0] < mintime:
                    mintime = rtime[sumkey][0]
                if index == 0:
                    index +=1
                if rtime[sumkey][index] > maxtime:
                    maxtime = rtime[sumkey][index]
        return mintime,maxtime
    
    st.printstart()
    
    # axis line width
    axwidth = 1.0
    # fontsize
    fs = 12
    # left-right scalar for scan info placement
    infop = 0.2
    
    if save < n: # if the script is told to save more often than it sums
        save = n
    
    mskeys = ['+','-']
    for key in sp: # append list places for chrom, summed chrom, and normalized chrom
        sp[key]['raw'] = []
        sp[key]['spectrum'] = Spectrum(3,startmz=sp[key]['bounds'][0],endmz=sp[key]['bounds'][1])
        sp[key]['%s' %(str(n)+'sum')] = []
        sp[key]['%s' %(str(n)+'norm')] = []
    
    mzml = mzML(filename) # load mzML class
    
    sp,TIC,rtime = mzml.pullspeciesdata(sp) # integrate species
    spec,sr,mz = mzml.pullspectra(mzrange=mz) # pull all spectra
    
    # run combine, regardless if called for (in order for keys to be correct
    #if n > 1: # run combine functions if n > 1
    sys.stdout.write('%s summing and normalizing species traces' %str(n))
    sumkey = str(n)+'sum'
    normkey = str(n)+'norm'
    sumsp = []
    for key in sp:
        sp[key][sumkey] = bindata(n,1,sp[key]['raw']) # bin each species
        sp[key]['colour'] = Colour(sp[key]['colour']).mpl # convert colour into matplotlib format
        for ind,val in enumerate(sp[key][sumkey]): # for normalization
            try:
                sumsp[ind] += val
            except IndexError:
                sumsp.append(val)
       
    for mode in mskeys: 
        sumkey = str(n)+'sum'+mode
        modekey = 'raw'+mode
        if modekey in rtime.keys(): # if there is data for that mode
            rtime[sumkey] = bindata(n,n,rtime[modekey])
            for ind,val in enumerate(rtime[sumkey]):
                rtime[sumkey][ind] = val - inj # shift time data to zero at injection point
            TIC[sumkey] = bindata(n,1,TIC[modekey])
            for key in sp: # for each species
                if sp[key]['affin'] in mskeys: # if species has affinity
                    spkey = str(n)+'sum'
                    sp[key][normkey] = []
                    for ind,val in enumerate(sp[key][spkey]):
                        sp[key][normkey].append(val/(sumsp[ind]+0.01)) #+0.01 to avoid div/0 errors
    sys.stdout.write(' DONE\n')
    sys.stdout.flush()

    spec = binnspectra(spec,n,startmz=mz[0],endmz=mz[1]) # bin mass spectra

    if os.path.isdir('imgs') == False: # check for /img directory and create if missing
        os.makedirs('imgs')
    for ind,val in enumerate(spec):
        curspec = ind*n+1
        if curspec >= scr[0] and curspec <= scr[1]: # if index is within scanrange to output
            sys.stdout.write('\rRendering scan #%i %.1f%% (scan range: %i to %i)' %(curspec,(float(curspec)-float(scr[0]))/(float(scr[1])-float(scr[0]))*100.,scr[0],scr[1]))
            val[1] = msfignorm(val[0],val[1]) # normalize spectrum
            
            itr = str(100000+curspec)
            mintime,maxtime = timelimits(ind)
            plotit(val[0],val[1],ind)
    sys.stdout.write(' DONE\n')
    st.printend()
    st.printprofiles()
        
    
if __name__ == '__main__':
    try:
        spectrumtrace(filename,sp,scr=scr,n=n,mz=mz,inj=inj,save=save)
    finally:
        import gc,sys
        sys.stdout.write('\nClearing memory')
        gc.collect()
        sys.stdout.write(' DONE\n')
    sys.stdout.write('fin.')