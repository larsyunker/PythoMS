"""
Tome v02 A compilation of all Lars' python scripts as callable functions

functions:
    autoresolution (estimates the resolution of a spectrum)
    bindata (bins a list of values)
    binnspectra (bins n mass spectra together into a single mass spectrum)
    bincidspectra (bins mass spectra together based on their collision voltage)
    filepresent (checks for a file or directory in the current working directory)
    find_all (finds all locations of files of a given name in a provided directory)
    linmag (generates a list of values which is linear in magnification)
    linramp (generates a list of values which is linear from start to finish)
    locateinlist (locates a value or the closest value to it in a sorted list)
    lyround (rounds a number given a particular base number)
    mag (calculates and returns the magnification of a given y value relative to the start)
    normalize (normalizes a list to a given value)
    plotms (plots a mass spectrum)
    sigmafwhm (cacluates sigma and fwhm from a resolution and a mass)
    strtolist (converts a string to a list)  
    version_input (uses the appropriate user input function depending on the python version)      

changelog:
    created mzML class and moved many functions to work within that class (removed several functions from Tome)
    added strtolist
    moved classes to separate files
    fullspeclist has been moved to _Spectrum class (there were issues with mutation of the original)
    calcindex has also been moved to _Spectrum class (it is used solely in that class)
    moved colours to _Colour class
    removed automz (now handled in the Molecule class)
    created bincidspectra to bin spectra with the same cid together
    removed loadwb, openpyxlcheck, pullparams (now included in XLSX class)
    generalized filepresent
    removed pwconvert (now included in mzML class)
    completely rewrote resolution
    rewrote resolution again to check multiple portions of the spectrum
    significant change to plotms
    moved alpha to XLSX class
    ---v02---
"""
# ----------------------------------------------------------
# -------------------FUNCTION DEFINITIONS-------------------
# ----------------------------------------------------------

def autoresolution(x,y,v=True):
    """
    determines the resolution of a provided spectrum
    
    x: list
        list of x values
    y: list
        list of y values (paired with x values
    v: bool
        verbose toggle
    
    resolution is based on the average resolution of 10 pseudo-random samples
    each sample spectrum is split into 10 sections, finding 10 peaks in order to calculate the resolution
    """
    def findsomepeaks(y,n=10):
        """roughly locates n peaks by maximum values in the spectrum and returns their index"""
        split = int(len(y)/n)
        start = 0
        end = start+split
        splity = []
        for i in range(n):
            splity.append(sci.asarray(y[start:end]))
            start += split
            end += split
        out = []
        for ind,section in enumerate(splity):
            maxy = max(section)
            if maxy == max(section[1:-1]): # if max is not at the edge of the spectrum
                out.append(sci.where(section==maxy)[0][0]+split*ind)
        return out
        
    def resolution(x,y,index=None):
        """
        Finds the resolution and full width at half max of a spectrum
        x: list of mz values
        y: corresponding list of intensity values
        index: index of maximum intensity (optional; used if the resolution of a specific peak is desired)
        
        returns resolution
        """
        y = sci.asarray(y) # convert to array for efficiency
        if index is None: # find index and value of maximum
            maxy = max(y)
            index = sci.where(y==maxy)[0][0]
        else:
            maxy = y[index]
        if maxy/(sum(y)/len(y)) < 10: # if intensity to average is below this threshold (rough estimate of signal/noise)
            return None
        halfmax = maxy/2
        indleft = int(index)-1 # generate index counters for left and right walking
        indright = int(index)+1
        while y[indleft] > halfmax: # while intensity is still above halfmax
            indleft -= 1
        while y[indright] > halfmax:
            indright += 1
        return x[index]/(x[indright]-x[indleft]) # return resolution (mz over full width at half max)
    
    import scipy as sci
    if v is True:
        import sys
        sys.stdout.write('\rEstimating resolution of the spectrum')
    
    inds = findsomepeaks(y) # find some peaks in the spectrum
    res = []
    for ind in inds: # for each of those peaks
        res.append(resolution(x,y,ind))
    res = [y for y in res if y is not None] # removes None values (below S/N)
    res = sum(res)/len(res) # calculate average
    if v is True:
        sys.stdout.write(': %.1f\n' %res)
    return res # return average

def bindata(n,v,lst):
    """
    Function for summing a supplied list of data (binning)
    
    input (n,v,lst,name)
    n is number of values to sum
    v is equal to n if an average value is required (e.g. for time values, usually it is equal to 1)
    lst is the list of values for combination
    """
    out = []
    delta = 0
    ttemp = 0
    for ind,val in enumerate(lst):
        delta += 1
        ttemp += val # add current value to growing sum
        if delta == n: # critical number is reached
            out.append(ttemp/v) # append sum to list
            delta = 0 # reset critical count and sum
            ttemp = 0
    return out

def binnspectra(dct,n,dec=3,startmz=50.,endmz=2000.):
    """
    sums n mass spectra together into a single spectrum
    
    dct:
        dictionary of timepoints (see mzML.pullspectra output)
        each timepoint is a dictionary with 'x' and 'y' lists
    n:
        number of scans to sum
    dec:
        how many decimals to keep
        default 3
        increasing this does not typically yield better spectra
    startmz:
        lowest mz to sum
        default 50
    endmz:
        highest mz to sum
        default 2000
    """
    import sys
    from _Spectrum import Spectrum
    out = []
    delta = 0
    spec = Spectrum(dec,startmz=startmz,endmz=endmz)
    for time in dct: # for each timepoint
        delta += 1
        sys.stdout.write('\rBinning spectrum #%i/%i  %.1f%%' %(delta,len(dct),float(delta)/float(len(dct))*100.))
        spec.addspectrum(dct[time]['x'],dct[time]['y']) # add spectrum
        if delta == n: # critical number is reached
            out.append(spec.trim(zeros=True)) # append list
            spec.resety() # reset y list in object
            delta = 0 # reset critical sum
    sys.stdout.write(' DONE\n')
    return out

def bincidspectra(dct,dec=3,startmz=50.,endmz=2000.,threshold=0,fillzeros=False):
    """
    bins mass spectra together based on their collision voltage
    
    lst:
        list of dictionarys corresponding to each scan (see mzML.pullspectra output)
    dec:
        how many decimals to keep
        default 3
        increasing this does not typically yield better spectra
    startmz:
        lowest mz to sum
        default 50
    endmz:
        highest mz to sum
        default 2000
    """
    from _Spectrum import Spectrum
    import sys
    binned = {}
    for time in dct:
        #sys.stdout.write('\rBinning spectrum by CID value #%i/%i  %.1f%%' %(ind+1,len(lst),float(ind+1)/float(len(lst))*100.))
        if binned.has_key(dct[time]['CE']) is False: # generate key and spectrum object if not present
            binned[dct[time]['CE']] = Spectrum(dec,startmz=startmz,endmz=endmz)
        else: # otherwise add spectrum
            binned[dct[time]['CE']].addspectrum(dct[time]['x'],dct[time]['y'])
    
    if threshold > 0 or fillzeros is True: # if manipulation is called for
        for vol in binned: # for each voltage
            sys.stdout.write('\rZero filling spectrum for %s eV' %`vol`)
            if threshold > 0:
                binned[vol].threshold(threshold) # apply threshold
            if fillzeros is True:
                binned[vol].fillzeros() # fill with zeros
        sys.stdout.write(' DONE\n')
    
    cv = [] # list for collision voltages
    specout = [] # list for spectra
    for vol,spec in sorted(binned.items()):
        sys.stdout.write('\rTrimming spectrum for %s eV' % `vol`)
        cv.append(vol) # append voltage to list
        specout.append(spec.trim()) # append trimmed spectrum to list
    sys.stdout.write(' DONE\n')
    sys.stdout.flush()
    return specout,cv
    
def filepresent(filename,ftype='file'):
    """
    function to check the presense of a file ('file') or directory ('dir')
    in the current working directory
    """
    import os
    if ftype == 'dir':
        if os.path.isdir(filename) == False:
            raise IOError('\nThe directory "%s" could not be located in the current working directory'%(filename))
    if ftype == 'file':
        if os.path.isfile(filename) == False:
            raise IOError('\nThe file "%s" could not be located in the current working directory'%(filename))

def find_all(fname,path):
    """
    Finds all files of a given name within a specified directory.
    Adapted from http://stackoverflow.com/questions/1724693/find-a-file-in-python
    
    Module dependancies: os
    """
    import os
    locations = []
    for root,dirs,files in os.walk(path):
        if fname in files:
            locations.append(os.path.join(root,fname))                   
    return locations

def linmag(vali,magstart,magend,dur):
    """
    funciton for generating a ramp of values that is linear in magnification
    vali: initial value (globally)
    magstart: magnification at the start of the ramp
    magend: magnification at the end of the ramp
    dur: number of steps (duration)
    """
    out = []
    for i in range(dur):
        out.append(float(vali)/((magend-magstart)/dur*i + magstart))
    return out

def linramp(valstart,valend,dur):
    """
    Function for generating a linear ramp of values
    valstart: value at the start of the ramp
    valend: value at the end of the ramp
    dur: number of steps (duration)
    """
    out = []
    for i in range(int(dur)):
        out.append( ((float(valend-valstart))/(float(dur)))*(i) + valstart )
    return out

def locateinlist(lst,value,bias='closest'):
    """
    Finds index in a sorted list of the value closest to a given value
    
    If two numbers are equally close, return the smallest number.
    based on http://stackoverflow.com/questions/12141150/from-list-of-integers-get-number-closest-to-a-given-value
    
    lst: list
        list of values to search
    value: float or int
        value number to find
    bias: 'lesser','greater', or 'closest'
        default 'left'
        'lesser' will return the position of the value just less than the provided value
        'greater' will return the position of the value just greater than the provided value
        'closest' will return the index of the nearest value to the one provided
    """                
    pos = self.bl(lst, value)
    if pos == 0: # if at start of list
        return pos
    elif pos == len(lst): # if insertion is beyond index range
        return pos -1 
    if lst[pos] == value: # if an exact match is found
        return pos
    if bias == 'greater': # return value greater than the value (bisect_left has an inherent bias to the right)
        return pos
    if bias == 'lesser': # return value lesser than the provided
        return pos -1
    if bias == 'closest': # check differences between index and index-1 and actual value, return closest
        adjval = abs(lst[pos-1] - value)
        curval = abs(lst[pos] - value)
        if adjval < curval: # if the lesser value is closer
            return pos-1
        if adjval == curval: # if values are equidistant
            return pos-1
        else:
            return pos

def lyround(x,basen):
    """
    Function for rounding a given number using a specific base
    based on http://stackoverflow.com/questions/2272149/round-to-5-or-other-number-in-python
    """
    base = basen**(int(len(str(int(x))))-1)
    return int(base * round(float(x)/base))

def mag(initial,final):
    """
    calculate magnification for a value
    """
    return float(initial)/float(final)

def normalize(lst,maxval):
    """
    function for normalizing a list of values to a given value
    input:
        lst: list of values
        maxval: value to normalize to
    """
    listmax = max(lst)
    for ind,val in enumerate(lst):
        lst[ind] = float(val)/float(listmax)*maxval
    return lst

def plotms(realspec,simdict={},**kwargs):
    """
    plots a formatted mass spectrum and optionally overlays predicted isotope patterns
    
    realspec: the spectrum to plot
        list of [[mzvalues],[intvalues]]
    
    simdict: molecular formulae to predict and overlay
        can be handed several options
        - a dictionary of formulas (with each formula being its own dictionary with colour and alpha keys)
        - a list of formulas (all overlays will default to black with 0.5 alpha)
        - a single formula (black with 0.5 alpha)
    
    kwargs:
        many parameters can be changed to tweak the appearance of the plot
        see the settings dictionary in this function for details     
    """
    def localmax(x,y,xval,lookwithin=1):
        """finds the local maximum within +/- lookwithin of the xval"""
        l,r = bl(x,xval-lookwithin),br(x,xval+lookwithin)
        return max(y[l:r])
    
    def trimspectrum(x,y,left,right):
        """trims a spectrum to the left and right bounds"""
        l,r = bl(x,left),br(x,right) # find indicies
        return x[l:r],y[l:r] # trim spectrum
    
    def estimatedem(x,y,em,simmin,simmax,lookwithin=1):
        """estimates the exact mass of a peak"""
        l,r = bl(x,simmin-lookwithin),br(x,simmax+lookwithin) # narrow range to that of the isotope pattern
        print em
        print x[l:r]
        locmax = max(y[l:r]) # find local max in that range
        for ind,val in enumerate(y):
            if val == locmax: # if the y-value equals the local max
                print x[ind]
                if ind >= l and ind <= r: # and if the index is in the range (avoids false locations)
                    return x[ind]
        difleft = abs(em-simmin)
        difright = abs(em-simmax)
        print difleft,difright
        return '>%.1f' %max(difleft,difright) # if no match is found, return maximum difference
    
    def checksimdict(dct):
        """
        checks the type of simdict, converting to dictionary if necessary
        also checks for alpha and colour keys and adds them if necessary (defaulting to key @ 0.5)
        """
        if type(dct) is not dict:
            if type(dct) is str:
                dct = {dct:{}}
            elif type(dct) is list or type(dct) is tuple:
                tdct = {}
                for i in dct:
                    tdct[i] = {}
                dct = tdct
        for species in dct:
            if dct[species].has_key('colour') is False:
                dct[species]['colour'] = 'k'
            if dct[species].has_key('alpha') is False:
                dct[species]['alpha'] = 0.5
        return dct
            
    import sys
    from _classes._Colour import Colour
    from _classes._Molecule import Molecule
    from tome_v02 import autoresolution,normalize
    import pylab as pl
    from bisect import bisect_left as bl
    from bisect import bisect_right as br
    
    settings = { # default settings
    'mz':'auto', # m/z bounds for the output spectrum
    'outname':'spectrum', # name for the output file
    'output':'save', # 'save' or 'show' the figure
    'simtype':'bar', # simulation overlay type ('bar' or 'gaussian')
    'spectype':'continuum', # spectrum type ('continuum' or 'centroid')
    'maxy':'max', # max or value
    'norm':True, # True or False
    'simnorm':'spec', # top, spec, or value
    'xlabel':True, # show x label
    'ylabel':True, # show y label
    'xvalues':True, #show x values
    'yvalues':True, # show y values
    'showx':True, # show x axis
    'showy':True, # how y axis
    'offsetx':True, # offset x axis (shows low intensity species better)
    'fs':16, # font size
    'lw':1.5, # line width for the plotted spectrum
    'axwidth':1.5, # axis width 
    'simlabels':False, # show labels isotope for patterns
    'bw':'auto', # bar width for isotope patterns (auto does 2*fwhm)
    'specfont':'Arial', # the font for text in the plot
    'size':[7.87,4.87], # size in inches for the figure
    'dpiout':300, # dpi for the output figure
    'exten':'png', # extension for the output figure
    'res':False, # output the resolution of the spectrum
    'delta':False, # output the mass delta between the spectrum and the isotope patterns
    'stats':False, # output the goodness of match between the spectrum and the predicted isotope patterns,
    'speccolour':'k', # colour for the spectrum to be plotted
    'padding':'auto', # padding for the output plot
    'verbose':True # verbose setting
    }
    
    if set(kwargs.keys()) - set(settings.keys()): # check for invalid keyword arguments
        string = ''
        for i in set(kwargs.keys()) - set(settings.keys()):
            string += ` i`
        raise KeyError('Unsupported keyword argument(s): %s' %string)
    
    settings.update(kwargs) # update settings from keyword arguments
    
    res = autoresolution(realspec[0],realspec[1]) # calculate resolution
    
    simdict = checksimdict(simdict) # checks the simulation dictionary
    for species in simdict: # generate Molecule object and set x and y lists
        simdict[species]['colour'] = Colour(simdict[species]['colour'])
        simdict[species]['mol'] = Molecule(species, res=res) 
        if settings['simtype'] == 'bar':
            simdict[species]['x'],simdict[species]['y'] = simdict[species]['mol'].barip
        if settings['simtype'] == 'gaussian':
            simdict[species]['mol'].gaussianisotopepattern()
            simdict[species]['x'],simdict[species]['y'] = simdict[species]['mol'].gausip
        
    if settings['mz'] == 'auto': # automatically determine m/z range
        if settings['verbose'] is True:
            sys.stdout.write('Automatically determining m/z window')
        mz = [10000000,0]
        for species in simdict:
            simdict[species]['bounds'] = simdict[species]['mol'].bounds() # calculate bounds
            if simdict[species]['bounds'][0] < mz[0]:
                mz[0] = simdict[species]['bounds'][0]-1
            if simdict[species]['bounds'][1] > mz[1]:
                mz[1] = simdict[species]['bounds'][1]+1
        if mz == [10000000,0]:
            mz = [min(realspec[0]),max(realspec[0])]
        settings['mz'] = mz
        if settings['verbose'] is True:
            sys.stdout.write(': %i - %i\n' %(int(mz[0]),int(mz[1])))
            sys.stdout.flush()
    
    realspec[0],realspec[1] = trimspectrum(realspec[0],realspec[1],settings['mz'][0]-1,settings['mz'][1]+1) # trim real spectrum for efficiency
    
    if settings['norm'] is True: # normalize spectrum
        realspec[1] = normalize(realspec[1],100.)
    
    for species in simdict: # normalize simulations
        if settings['simnorm'] == 'spec': # normalize to maximum around exact mass
            simdict[species]['y'] = normalize(simdict[species]['y'],localmax(realspec[0],realspec[1],simdict[species]['mol'].em,simdict[species]['mol'].fwhm))
        elif settings['simnorm'] == 'top': # normalize to top of the y value
            if settings['maxy'] == 'max':
                raise ValueError('Simulations con only be normalized to the top of the spectrum when the maxy setting is a specific value')
            simdict[species]['y'] = normalize(simdict[species]['y'],settings['maxy'])
        elif type(settings['simnorm']) is int or type(settings['simnorm']) is float: # normalize to specified value
            simdict[species]['y'] = normalize(simdict[species]['y'],settings['simnorm'])
        if settings['delta'] is True:
            est = estimatedem(realspec[0],realspec[1],simdict[species]['mol'].em,min(simdict[species]['x']),max(simdict[species]['x'])) # try to calculate exact mass
            if type(est) is float:
                simdict[species]['delta'] = simdict[species]['mol'].em - est
            else:
                simdict[species]['delta'] = est
    
    pl.clf() # clear and close figure if open
    pl.close()
    fig = pl.figure(figsize = tuple(settings['size']))
    ax = fig.add_subplot(111)
    
    ax.spines["right"].set_visible(False) #hide right and top spines
    ax.spines["top"].set_visible(False)
    
    if settings['showx'] is False: 
        ax.spines["bottom"].set_visible(False) # hide bottom axis
    if settings['showy'] is False:
        ax.spines["left"].set_visible(False) # hide left axis
    
    for axis in ["top","bottom","left","right"]:
        ax.spines[axis].set_linewidth(settings['axwidth'])
    
    if settings['offsetx'] is True: # offset x axis
        ax.spines["bottom"].set_position(('axes',-0.01))  
    
    font = {'fontname':settings['specfont'],'fontsize':settings['fs']} #font parameters for axis/text labels
    tickfont = pl.matplotlib.font_manager.FontProperties(family=settings['specfont'],size=settings['fs']) # font parameters for axis ticks
    
    ax.set_xlim(settings['mz']) # set x bounds
    
    if settings['maxy'] == 'max': # set y bounds
        ax.set_ylim((0,max(realspec[1])))
        top = max(realspec[1])
    elif type(settings['maxy']) is int or type(settings['maxy']) is float:
        ax.set_ylim((0,settings['maxy']))
        top = settings['maxy']
    
    if settings['simtype'] == 'bar': # generates zeros for bottom of bars (assumes m/z spacing is equal between patterns)
        for species in simdict: 
            simdict[species]['zero'] = []
            for i in simdict[species]['x']:
                simdict[species]['zero'].append(0.)
        for species in simdict: # for each species
            for subsp in simdict: # look at all the species
                if subsp != species: # if it is not itself
                    ins = bl(simdict[subsp]['x'],simdict[species]['x'][-1]) # look for insertion point
                    if ins > 0 and ins < len(simdict[subsp]['x']): # if species highest m/z is inside subsp list
                        for i in range(ins): # add intensity of species to subsp zeros
                            simdict[subsp]['zero'][i] += simdict[species]['y'][-ins+i]
    if settings['res'] is True and settings['spectype'] != 'centroid': #include resolution if specified (and spectrum is not centroid)
        ax.text(mz[1],top*0.95,'resolution: '+str(round(res))[:-2],horizontalalignment='right',**font)
    
    for species in simdict: # plot and label bars
        if settings['simtype'] == 'bar':
            if settings['bw'] == 'auto':
                bw = simdict[species]['mol'].fwhm*2
            else:
                bw = settings['bw']
            ax.bar(simdict[species]['x'], simdict[species]['y'], bw, alpha = simdict[species]['alpha'], color = simdict[species]['colour'].mpl, linewidth=0, align='center',bottom=simdict[species]['zero'])
        elif settings['simtype'] == 'gaussian':
            ax.plot(simdict[species]['x'], simdict[species]['y'], alpha = simdict[species]['alpha'], color = simdict[species]['colour'].mpl, linewidth=settings['lw'])
            ax.fill(simdict[species]['x'], simdict[species]['y'], alpha = simdict[species]['alpha'], color = simdict[species]['colour'].mpl, linewidth=0)
        if settings['simlabels'] is True or settings['stats'] is True or settings['delta'] is True: # if any labels are to be shown
            string = ''
            bpi = simdict[species]['y'].index(max(simdict[species]['y'])) # index of base peak
            if settings['simlabels'] is True: # species name
                string += species
                if settings['stats'] is True or settings['delta'] is True: # add return if SER or delta is called for
                    string += '\n'
            if settings['stats'] is True: # standard error of regression
                string += 'SER: %.2f ' %simdict[species]['mol'].compare(realspec)
            if settings['delta'] is True: # mass delta
                string += 'mass delta: '
                if type(simdict[species]['delta']) is float:
                    string += '%.3f' %simdict[species]['delta']
                else:
                    string += '%s' %simdict[species]['delta']
            ax.text(simdict[species]['x'][bpi],top*(1.01),string, color = simdict[species]['colour'].mpl, horizontalalignment='center', **font)
    
    if settings['spectype'] == 'continuum':
        ax.plot(realspec[0], realspec[1], linewidth=settings['lw'], color=Colour(settings['speccolour']).mpl)
    elif settings['spectype'] == 'centroid':
        dist = []
        for ind,val in enumerate(realspec[0]): # find distance between all adjacent m/z values
            if ind == 0:
                continue
            dist.append(realspec[0][ind]-realspec[0][ind-1])
        dist = sum(dist)/len(dist) # average distance
        ax.bar(realspec[0], realspec[1], dist*0.75, linewidth=0, color=Colour(settings['speccolour']).mpl, align='center', alpha=0.8)

    # show or hide axis values/labels as specified
    if settings['yvalues'] is False: # y tick marks and values
        ax.tick_params(axis='y', labelleft='off',length=0)
    if settings['yvalues'] is True: # y value labels
        ax.tick_params(axis='y', length=settings['axwidth']*3, width=settings['axwidth'], direction='out',right='off')
        for label in ax.get_yticklabels():
            label.set_fontproperties(tickfont)
    if settings['ylabel'] is True: # y unit
        if top == 100: # normalized
            ax.set_ylabel('relative intensity', **font)
        else: # set to counts
            ax.set_ylabel('intensity (counts)', **font)
            
    if settings['xvalues'] is False:  # x tick marks and values
        ax.tick_params(axis='x', labelbottom='off',length=0)
    if settings['xvalues'] is True: # x value labels
        ax.tick_params(axis='x', length=settings['axwidth']*3, width=settings['axwidth'] ,direction='out',top = 'off')
        for label in ax.get_xticklabels():
            label.set_fontproperties(tickfont) 
    if settings['xlabel'] is True: # x unit
        ax.set_xlabel('m/z', style='italic', **font)
    
    if settings['padding'] == 'auto':
        pl.tight_layout(pad=0.5) # adjust subplots
        if settings['simlabels'] is True or settings['stats'] is True or settings['delta'] is True: 
            pl.subplots_adjust(top = 0.90) # lower top if details are called for
    elif type(settings['padding']) is list and len(settings['padding']) == 4:
        pl.subplots_adjust(left=settings['padding'][0], right=settings['padding'][1], bottom=settings['padding'][2], top=settings['padding'][3])
    
    if settings['output'] == 'save': # save figure
        outname = '' # generate tag for filenaming
        for species in simdict:
            outname+=' '+species
        outname = settings['outname'] + outname + '.' + settings['exten']
        pl.savefig(outname, dpi=settings['dpiout'], format=settings['exten'], transparent=True)
        if settings['verbose'] is True:
            sys.stdout.write('Saved figure as:\n"%s"\nin the working directory' %outname)
    
    elif settings['output'] == 'show': # show figure
        pl.show()

def plotuv(wavelengths,intensities,**kwargs):
    """
    plots a UV-Vis spectrum
    input:
    wavelengths: list
        list of wavelengths
    intensities: list or list of lists
        list of intensities (matching wavelengths list)
        can also supply several lists of intensities (to plot a progressive UV plot)
    kwargs:
        see settings for supported kwargs and what they do
    """
    settings = { # default settings for the function
    'outname':'UV-Vis spectrum', # name for the output file
    'fs':16, # font size
    'lw':1.5, # line width for the plotted spectrum
    'axwidth':1.5, # axis width 
    'size':[7.87,4.87], # size in inches for the figure
    'dpiout':300, # dpi for the output figure
    'exten':'png', # extension for the output figure
    'specfont':'Arial', # the font for text in the plot
    # colours to use for multiple traces in the same spectrum (feel free to specify your own)
    'colours':['#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99','#8dd3c7','#ffffb3','#bebada','#fb8072','#80b1d3','#fdb462','#b3de69','#fccde5','#d9d9d9','#bc80bd','#ccebc5',], 
    'xrange':None, # the limits for the x axis
    'yrange':None, # the limits for the y axis
    'times':None, # time points for each provided trace (for legend labels)
    'output':'save', # 'save' or 'show' the figure
    'padding':None, # padding for the output plot
    'verbose':True, # chatty
    'legloc':0, # legend location (see http://matplotlib.org/api/legend_api.html location codes)
    }
    if set(kwargs.keys()) - set(settings.keys()): # check for invalid keyword arguments
        string = ''
        for i in set(kwargs.keys()) - set(settings.keys()):
            string += ` i`
        raise KeyError('Unsupported keyword argument(s): %s' %string)
    
    settings.update(kwargs) # update settings from keyword arguments
    
    import sys
    import pylab as pl
    from _classes._Colour import Colour
    pl.clf() # clear and close figure if open
    pl.close()
    fig = pl.figure(figsize = tuple(settings['size']))
    ax = fig.add_subplot(111)
    
    ax.spines["right"].set_visible(False) #hide right and top spines
    ax.spines["top"].set_visible(False)
    
    font = {'fontname':settings['specfont'],'fontsize':settings['fs']} #font parameters for axis/text labels
    tickfont = pl.matplotlib.font_manager.FontProperties(family=settings['specfont'],size=settings['fs']) # font parameters for axis ticks
    
    if type(intensities[0]) is float: # if the function has only been handed a single spectrum
        intensities = [intensities]
    
    # determine and set limits for axes
    if settings['xrange'] is None: # auto determine x limits
        settings['xrange'] = [min(wavelengths),max(wavelengths)]
    if settings['yrange'] is None: # auto determine y limits
        settings['yrange'] = [0,0]
        for spec in intensities:
            if max(spec) > settings['yrange'][1]:
                settings['yrange'][1] = max(spec)
    ax.set_xlim(settings['xrange']) # set x bounds
    ax.set_ylim(settings['yrange']) # set y bounds
    
    # apply font and tick parameters to axes
    ax.tick_params(axis='x', length=settings['axwidth']*3, width=settings['axwidth'] ,direction='out',top = 'off')
    for label in ax.get_xticklabels():
        label.set_fontproperties(tickfont) 
    ax.tick_params(axis='y', length=settings['axwidth']*3, width=settings['axwidth'], direction='out',right='off')
    for label in ax.get_yticklabels():
        label.set_fontproperties(tickfont)
    for axis in ["top","bottom","left","right"]:
        ax.spines[axis].set_linewidth(settings['axwidth'])
    
    if settings['times'] is not None:
        if len(settings['times']) != len(intensities):
            raise IndexError('The numer of times provided do not match the number of traces provided.')
    
    for ind,spec in enumerate(intensities): # plot traces
        if settings['times'] is not None:
            string = 't = '+str(round(settings['times'][ind],1))+'m'
            ax.plot(wavelengths,spec,label=string,color=Colour(settings['colours'][ind]).mpl,linewidth=settings['lw'])
        else:
            ax.plot(wavelengths,spec,color=Colour(settings['colours'][ind]).mpl,linewidth=settings['lw'])
    
    if settings['times'] is not None:
        ax.legend(loc=0,frameon=False)
    
    ax.set_xlabel('wavelength (nm)', **font)
    ax.set_ylabel('absorbance (a.u.)', **font)
    
    if settings['padding'] is None:
        pl.tight_layout(pad=0.5) # adjust subplots
    elif type(settings['padding']) is list and len(settings['padding']) == 4:
        pl.subplots_adjust(left=settings['padding'][0], right=settings['padding'][1], bottom=settings['padding'][2], top=settings['padding'][3])
    
    if settings['output'] == 'save': # save figure
        outname = settings['outname'] + '.' + settings['exten']
        pl.savefig(outname, dpi=settings['dpiout'], format=settings['exten'], transparent=True)
        if settings['verbose'] is True:
            sys.stdout.write('Saved figure as:\n"%s"\nin the working directory' %outname)
    
    elif settings['output'] == 'show': # show figure
        pl.show()

def sigmafwhm(res,em):
    """
    determines the full width at half max and sigma for a normal distribution
    res is the resolution of the instrument
    em is the mass being calculated
    """
    import math
    fwhm = em/res
    sigma = fwhm/(2*math.sqrt(2*math.log(2))) # based on the equation FWHM = 2*sqrt(2ln2)*sigma
    return fwhm,sigma

def strtolist(string):
    """
    converts a string to a list with more flexibility than string.split()
    """
    out = []
    temp = ''
    brackets = ['(',')','[',']','{','}']
    for char in list(string):
        if char not in brackets and char != ',':
            temp += char
        if char == ',':
            try:
                out.append(int(temp))
            except ValueError:
                out.append(float(temp))
            temp = ''
    if len(temp) != 0: # if there is a weird ending character
        try:
            out.append(int(temp))
        except ValueError:
            out.append(float(temp))
    return out

def version_input(string):
    """checks the python version and uses the appropriate version of user input"""
    import sys
    if sys.version.startswith('2.7'):
        return raw_input('%s' %string)
    if sys.version.startswith('3.'):
        return input('%s' %string)
    else:
        raise EnvironmentError('The version_input method encountered an unsupported version of python.')