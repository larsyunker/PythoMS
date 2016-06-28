"""
Tome v02 A compilation of all Lars' python scripts as callable functions

functions:
    alpha (returns column index from numerical index for excel files)
    bindata (bins a list of values)
    binnspectra (bins n mass spectra together into a single mass spectrum)
    bincidspectra (bins mass spectra together based on their collision voltage)
    filepresent (checks for a file or directory in the current working directory)
    find_all (finds all locations of files of a given name in a provided directory)
    linmag (generates a list of values which is linear in magnification)
    linramp (generates a list of values which is linear from start to finish)
    lyround (rounds a number given a particular base number)
    mag (calculates and returns the magnification of a given y value relative to the start)
    msipoverlay (generates a figure overlaying predicted isotope patterns on top of experimental data)
    normalize (normalizes a list to a given value)
    resolution (calculates the resolution of a spectrum peak)
    sigmafwhm (cacluates sigma and fwhm from a resolution and a mass)
    strtolist (converts a string to a list)        

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
    ---v02---


"""
# ----------------------------------------------------------
# -------------------FUNCTION DEFINITIONS-------------------
# ----------------------------------------------------------

def alpha(index):
    """
    function for taking an index and converting it to the corresponding column used by excel
    works up to column ZY (column # 701)
    """
    import string
    alphabet = list(string.ascii_lowercase)
    out = ''
    n = (index+1)/26
    if n != 0:
        out += alphabet[n-1]
        out += alphabet[(index+1)%26-1]
    else:
        out = alphabet[index]
    return out.upper()

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
    #for ind,scan in enumerate(lst): # for each pair
    #    delta += 1
    #    sys.stdout.write('\rBinning spectrum #%i/%i  %.1f%%' %(ind+1,len(lst),float(ind+1)/float(len(lst))*100.))
    #    spec.addspectrum(scan['x'],scan['y']) # add spectrum
    #    if delta == n: # critical number is reached
    #        out.append(spec.trim(zeros=True)) # append list
    #        spec.resety() # reset y list in object
    #        delta = 0 # reset critical sum
    #sys.stdout.write(' DONE\n')
    #return out

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
        if binned.has_key(dct[time]['CE']) is False: # generate key if not present
            binned[dct[time]['CE']] = Spectrum(dec,startmz=startmz,endmz=endmz)
        binned[dct[time]['CE']].addspectrum(dct[time]['x'],dct[time]['y'])
    #for ind,scan in enumerate(lst):
    #    sys.stdout.write('\rBinning spectrum by CID value #%i/%i  %.1f%%' %(ind+1,len(lst),float(ind+1)/float(len(lst))*100.))
    #    if binned.has_key(scan['col']) is False: # generate key if not present
    #        binned[scan['col']] = Spectrum(dec,startmz=startmz,endmz=endmz)
    #    binned[scan['col']].addspectrum(scan['x'],scan['y']) # add spectrum existing
    
    #sys.stdout.write(' DONE\n')
    
    if threshold > 0 or fillzeros is True: # if manipulation is called for
        for vol in binned:
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
    
    def estimatedem(x,y,em,lookwithin=1):
        """estimates the exact mass of a peak"""
        l,r = bl(x,em-lookwithin),br(x,em+lookwithin)
        # !!! figure out how to find the exact max value of y that was found
    
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
            
    import os,sys
    sys.path.append(os.path.dirname(os.path.realpath(__file__))+'/_classes')
    from _Colour import Colour
    from _Molecule import Molecule
    from tome_v02 import resolution,normalize
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
    'norm':False, # True or False
    'simnorm':'spec', # top, spec, or value
    'axlabels':[True,True], # show/hide axes labels
    'axvalues':[True,True], # show/hide axes values and tick marks
    'axshow':[True,True], # show/hide axes
    'offsetx':True, # offset x axis (shows low intensity species better)
    'fs':16, # font size
    'lw':1.5, # line width for the plotted spectrum
    'axwidth':1.5, # axis width 
    'simlabels':False, # show labels isotope for patterns
    'bw':2, # bar width for isotope patterns
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
    
    res = resolution(realspec[0],realspec[1]) # calculate resolution
    
    simdict = checksimdict(simdict) # checks the simulation dictionary
    for species in simdict: # generate Molecule object and set x and y lists
        simdict[species]['colour'] = Colour(simdict[species]['colour'])
        simdict[species]['mol'] = Molecule(species, res=res) 
        if settings['simtype'] == 'bar':
            simdict[species]['x'],simdict[species]['y'] = simdict[species]['mol'].barip
        if settings['simtype'] == 'gaussian':
            simdict[species]['x'],simdict[species]['y'] = simdict[species]['mol'].gausip
        
    if settings['mz'] == 'auto': # automatically determine m/z range
        if settings['verbose'] is True:
            sys.stdout.write('Automatically determining m/z window: ')
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
            sys.stdout.write('%i - %i\n' %(int(mz[0]),int(mz[1])))
            sys.stdout.flush()
    
    realspec[0],realspec[1] = trimspectrum(realspec[0],realspec[1],settings['mz'][0],settings['mz'][1]) # trim real spectrum for efficiency
    
    if settings['norm'] is True: # normalize spectrum
        realspec[1] = normalize(realspec[1],100.)
    
    for species in simdict: # normalize simulations
        if settings['simnorm'] == 'spec': # normalize to maximum around exact mass
            simdict[species]['y'] = normalize(simdict[species]['y'],localmax(realspec[0],realspec[1],simdict[species]['mol'].em,simdict[species]['mol'].fwhm))
        elif settings['simnorm'] == 'top': # normalize to top of the y value
            simdict[species]['y'] = normalize(simdict[species]['y'],settings['maxy'])
        elif type(settings['simnorm']) is int or type(settings['simnorm']) is float: # normalize to specified value
            simdict[species]['y'] = normalize(simdict[species]['y'],settings['simnorm'])
    
    pl.clf() # clear and close figure if open
    pl.close()
    fig = pl.figure(figsize = tuple(settings['size']))
    ax = fig.add_subplot(111)
    
    ax.spines["right"].set_visible(False) #hide right and top spines
    ax.spines["top"].set_visible(False)
    
    if settings['axshow'][0] is False: 
        ax.spines["bottom"].set_visible(False) # hide bottom axis
    if settings['axshow'][1] is False:
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

    if settings['res'] is True: #include resolution if specified
        ax.text(mz[1],top,'resolution: '+str(res),horizontalalignment='right',**font)
    if settings['delta'] is not False: #display mass delta if specified
        ax.text(mz[1],top*0.9,'mass delta: '+str(round(settings['delta'],3)),horizontalalignment='right',**font)
    
    for species in simdict: # plot and label bars
        if settings['simtype'] == 'bar':
            ax.bar(simdict[species]['x'], simdict[species]['y'], simdict[species]['mol'].fwhm*settings['bw'], alpha = simdict[species]['alpha'], color = simdict[species]['colour'].mpl, linewidth=0, align='center',bottom=simdict[species]['zero'])
        elif settings['simtype'] == 'gaussian':
            ax.plot(simdict[species]['x'], simdict[species]['y'], alpha = simdict[species]['alpha'], color = simdict[species]['colour'].mpl, linewidth=settings['lw'])
            ax.fill(simdict[species]['x'], simdict[species]['y'], alpha = simdict[species]['alpha'], color = simdict[species]['colour'].mpl, linewidth=0)
        if settings['simlabels'] is True: # if labels are to be shown
            bpm = max(simdict[species]['y']) # simulation base peak intensity
            bpi = simdict[species]['y'].index(bpm) # index of base peak
            if settings['stats'] is True:
                simdict['SER'] = simdict[species]['mol'].compare(realspec)
                ax.text(simdict[species]['x'][bpi],top*(1.01),'%s SER: %.2f'%(species,simdict[species]['SER']), color = simdict[species]['colour'].mpl, horizontalalignment='center', **font)
            else:
                ax.text(simdict[species]['x'][bpi],top*(1.01),species, color = simdict[species]['colour'].mpl, horizontalalignment='center', **font)
    
    if settings['spectype'] == 'continuum':
        ax.plot(realspec[0], realspec[1], linewidth=settings['lw'], color=Colour(settings['speccolour']).mpl)
    elif settings['spectype'] == 'centroid':
        ax.bar(realspec[0], realspec[1], 1.0, linewidth=0, color=Colour(settings['speccolour']).mpl)

    # show or hide axis values/labels as specified
    if settings['axvalues'][1] is False: # y tick marks and values
        ax.tick_params(axis='y', labelleft='off',length=0)
    if settings['axvalues'][1] is True: # y value labels
        ax.tick_params(axis='y', length=settings['axwidth']*3, width=settings['axwidth'], direction='out',right='off')
        for label in ax.get_yticklabels():
            label.set_fontproperties(tickfont)
    if settings['axlabels'][1] is True: # y unit
        if top == 100: # normalized
            ax.set_ylabel('relative intensity', **font)
        else: # set to counts
            ax.set_ylabel('intensity (counts)', **font)
            
    
    if settings['axvalues'][0] is False:  # x tick marks and values
        ax.tick_params(axis='x', labelbottom='off',length=0)
    if settings['axvalues'][0] is True: # x value labels
        ax.tick_params(axis='x', length=settings['axwidth']*3, width=settings['axwidth'] ,direction='out',top = 'off')
        for label in ax.get_xticklabels():
            label.set_fontproperties(tickfont) 
    if settings['axlabels'][0] is True: # x unit
        ax.set_xlabel('m/z', style='italic', **font)
    
    if settings['padding'] == 'auto':
        pl.tight_layout(pad=0.5) # adjust subplots
        pl.subplots_adjust(top = 0.95)
    elif type(settings['padding']) is list and len(settings['padding']) == 4:
        pl.subplots_adjust(left=settings['padding'][0], right=settings['padding'][1], bottom=settings['padding'][2], top=settings['padding'][3])
    
    if settings['output'] == 'save': # save figure
        outname = '' # generate tag for filenaming
        for species in simdict:
            outname+=' '+species
        outname = settings['outname'] + outname + '.' + settings['exten']
        pl.savefig(outname, dpi=settings['dpiout'], format=settings['exten'], transparent=True)
        if settings['verbose'] is True:
            sys.stdout.write('Saved figure as \n"%s"\nin the working directory' %outname)
    
    if settings['output'] == 'show': # show figure
        pl.show()

def resolution(x,y,index=None):
    """
    Finds the resolution and full width at half max of a spectrum
    x: list of mz values
    y: corresponding list of intensity values
    index: index of maximum intensity (optional; used if the resolution of a specific peak is desired)
    
    returns resolution
    """
    import scipy as sp
    y = sp.asarray(y) # convert to array for efficiency
    if index is None: # find index and value of maximum
        maxy = max(y)
        index = sp.where(y==maxy)[0][0]
    else:
        maxy = y[index]
    halfmax = maxy/2
    indleft = int(index)-1 # generate index counters for left and right walking
    indright = int(index)+1
    while y[indleft] > halfmax: # while intensity is still above halfmax
        #if y[indleft] > y[indleft]+1: # if intensity start increasing (half max was not reached)
        #    raise ValueError('Half max for the left side of the peak could not be determined')
        indleft -= 1
    while y[indright] > halfmax:
        #if y[indright] > y[indright]-1:
        #    raise ValueError('Half max for the right side of the peak could not be determined')
        indright += 1
    return x[index]/(x[indright]-x[indleft]) # return resolution (mz over full width at half max)

#def resolution(x,y,maxindex,realmax):
#    """
#    Function for finding the resolution of a peak given its index in a list and the maximum value
#    x: list of mz values
#    y: corresponding list of intensity values
#    maxindex: index of maximum intensity
#    realmax: value of maximum intensity
#    
#    returns: [halfmax,halfleft,leftmz,halfright,rightmz,delta mz, resolution]
#    """
#    out = []
#    out.append(realmax/2) #0 halfmax
#    for i in range(maxindex):
#        if y[maxindex-i] <= out[0]:
#            out.append(y[maxindex-i]) #1 halfleft
#            out.append(x[maxindex-i]) #2 leftmz
#            break
#    for i in range(len(x)-maxindex):
#        if y[maxindex+i] <= out[0]:
#            out.append(y[maxindex+i]) #3 halfright
#            out.append(x[maxindex+i]) #4 rightmz
#            break
#    out.append(out[4]-out[2]) #5 dmz
#    out.append(x[maxindex]/out[-1]) #6 resolution
#    return out

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
