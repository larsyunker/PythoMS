"""
  isotope pattern check program v011 beta (dependant on tome)
  overlays a supplied predicted isotope pattern onto an acquired spectrum
  
new/changed:
    colour function can now interpret hex codes
    added overrides for figure size and output type
    new small function to import wb (the script may no longer require the loading of openpyxl outside of this function)
    removed import of os (this appears to not be used anymore)
    ---009---
    switched to tome_v02
    converted simdict entries to be dictionaries
    ---010.0---
    reworked settings dictionaries to be easier to read (and changed script/tome code to work with the new keys)
    ---010.1---
    rewrittten to work with MolecularFormula class
    removed functionality that loads isotope patterns from excel file
    converted into a standalone function
    cleaned up code to use exp (experimental) instead of realmz,realint (also got rid of a few extraneous variables)
    ---011.0---
    added mass delta calculator
    fixed formula interpreter in _MolecularFormula
    fixed normalization to max with close isotope patterns
    added bar stacking of intensities if there is pattern overlap (assumes m/z spacing is equal between patterns)
    added 'detailed' setting for more verbose isotope pattern checks
    ---011.1---
    updated to work with NoneSpectrum 2.1 functionality
    fixed snipping (should now work at either end of the spectrum)
    updated to work with Molecule 2.3
    ---011.2---
    added standard error output (using the compare function in Molecule)
    consolidated multiple loops into one
    Molecule objects are stored in simdict
    ---011.3---
    updated to work with the latest functionality of Molecule
    ---011.4 building

to add:
    convert to use XLSX class
    get rid of realmax necessity for figure generation (redefine in sets['ymax'])
    create a more generic version for naming the output file in the tome function
    databridge for masslynx conversion (is this even useful?)
"""
# provide the experimental spectrum xlsx
# the script will automatically use the first sheet
spectrum = '2016-06-13 09 Pd cation with boronic acid'

# supply simulated spectrum in dictionary format
# 'molecular formula':['colour','alpha','charge']
simdict = {
#'L2PdAr+OMe':{'colour':'#fc8d59','alpha':0.75,'charge':1},
#'L2PdAr+':{'colour':'#e41a1c','alpha':0.75,'charge':2}
#'(MeOC6H4BO)2OH':{'colour':'#e41a1c','alpha':0.75,'charge':1}
'L2PdPh':{'colour':'#4f81bd','alpha':0.75}
}

# choose the desired m/z range for the figure
# [m/z left,m/z right] or set to 'Auto'
mz = 'Auto'
#mz = [1008,1028]

# choose a figure type for auto settings
# options: 'pub', 'pubsvg', 'inset', 'insetsvg', 'thesis', 'user' , 'detailed'
# if you want to modify parameters, change the "user" preset in the figuresetting function below
setting = 'detailed'

# size override (change for quick adjustment of the output figure size)
# input format: [xwidth,ywidth] or None
# values in inches
sovr = None

# file format override (change for quick adjustment of the output figure type)
# valid inputs: 'png' or 'svg' or None
tovr = None

# specify bar plot ('bar') or predicted gaussian ('gaussian') distribution
otype = 'gaussian'

def ipoverlay(specfile,simdict,mz='Auto',setting='user',sovr=None,tovr=None,otype='bar',outname='ipoverlay'):
    """
    Function for overlaying one or more predicted isotope patterns over an experimental mass spectrum
    
    exp: experimental spectrum
        [[m/z values],[intensity values]]
    
    simdict: dictionary of all the desired compounds to show an isotope pattern for
        simdict = {'formula of compound':{'color':###,'alpha':0.##,'charge':1}, ...}
    
    mz: override mz window from 'Auto'
        [m/z start,m/z end)
    
    setting:
        dictionary key defined in figuresetting()
    
    sovr: override output figure size
        [xwidth,ywidth] (values in inches)
    
    tovr: override output format
        'png' or 'svg'
    
    otype: overlay type
        'bar' or 'gaussian'
    """
    def figuresetting(typ):
        """
        Returns figure settings based on provided key for use with msfigure
        
        ymax: the maximum y axis value of the output chart
            '100' (for use with normalization to 100)
            'counts' (actual counts)
            ### (specific count value)    
        
        norm: how the script should normalize the simulations
            'max' (normalizes to the maximum value in the simulation mz range)
            'top' (normalizes to the maximum y value of the chart)
            ### (specified value)
        
        axlabel: show axis labels?
            True/False for [x axis,y axis]
        
        axvalues: show axis values?
            True/False for [x axis,y axis]
        
        axshow: show axis?
            True/False for [x axis,y axis]
        
        fs: desired fontsize
            integer value in points
        
        lw: linewidth for the plotted spectrum
            integer or float value (recommend 1.5 as a default)
        
        axwidth: linewidth for the axes
            integer or float value (recommend setting the same as linewidth)
        
        label: show simulation labels
            True/False
        
        bw: prdicted pattern bar width
            integer multiple of the FWHM of the base peak in the isotope pattern
        
        res: show calculated resolution
            [True/False,None] (None is replaced by the resolution and is required)
        
        specfont: specify font to be used
            Scott's preferred font for publications is Arial
        
        size: figure size (in inches)
            [width,height] (floating point values)
        
        dpiout: dots per inch for the output png
            integer value (recommend using 300 for high resolution)
        
        exten: extension for the output file
            'png' for flat images
            'svg' for use with vector programs
        
        delta: calculate mass delta between predicted and experimental mass
            [True/False,None] (None is replaced by the delta value and is required)
        
        stats: calculate and output the standard error of regression between the experimental and predicted isotope patterns
            True/False
        
        Additional preset can be added in dictionary form provided that they have all of the above keys
        """
        ipsetting = {
        'pub':{'ymax':'100','norm':'max','axlabels':[True,True],'axvalues':[True,True],'axshow':[True,True],'fs':16,'lw':1.5,'axwidth':1.5,'simlabels':False,'bw':2,'res':[False,None],'specfont':'Arial','size':[7.87,4.87],'dpiout':300,'exten':'png','delta':[False,None],'stats':False},
        'pubsvg':{'ymax':'100','norm':'max','axlabels':[True,True],'axvalues':[True,True],'axshow':[True,True],'fs':16,'lw':1.5,'axwidth':1.5,'simlabels':False,'bw':2,'res':[False,None],'specfont':'Arial','size':[7.87,4.87],'dpiout':300,'exten':'svg','delta':[False,None],'stats':False},
        'inset':{'ymax':'100','norm':'max','axlabels':[True,True],'axvalues':[True,True],'axshow':[True,True],'fs':12,'lw':1.5,'axwidth':1.5,'simlabels':False,'bw':2,'res':[False,None],'specfont':'Arial','size':[2.8,2.8],'dpiout':300,'exten':'png','delta':[False,None],'stats':False},
        'insetsvg':{'ymax':'100','norm':'max','axlabels':[True,False],'axvalues':[True,False],'axshow':[True,False],'fs':12,'lw':1.5,'axwidth':1.5,'simlabels':False,'bw':2,'res':[False,None],'specfont':'Arial','size':[2.8,2.8],'dpiout':300,'exten':'svg','delta':[False,None],'stats':False},
        'thesis':{'ymax':'counts','norm':'max','axlabels':[True,True],'axvalues':[True,True],'axshow':[True,True],'fs':16,'lw':1.5,'axwidth':1.5,'simlabels':False,'bw':2,'res':[False,None],'specfont':'Arial','size':[6.5,4.0],'dpiout':300,'exten':'png','delta':[False,None],'stats':False},
        'detailed':{'ymax':'counts','norm':'max','axlabels':[True,True],'axvalues':[True,True],'axshow':[True,True],'fs':10,'lw':1.5,'axwidth':1.0,'simlabels':True,'bw':2,'res':[True,None],'specfont':'Arial','size':[5,3],'dpiout':300,'exten':'png','delta':[True,None],'stats':True},
        'user':{'ymax':'counts','norm':'max','axlabels':[True,True],'axvalues':[True,True],'axshow':[True,True],'fs':10,'lw':1.5,'axwidth':1.0,'simlabels':False,'bw':2,'res':[False,None],'specfont':'Arial','size':[6.5,4.0],'dpiout':300,'exten':'png','delta':[True,None],'stats':False}
        }
        try:
            sys.stdout.write('Using figure preset "%s"\n' %(setting))
            sys.stdout.flush()
            return ipsetting[typ]
        except KeyError:
            sys.exit('\nThe specified figure setting "%s" is not defined.\nPlease check your spelling' %setting)
    
    
    
    # ----------------------------------------------------------
    # -------------------PROGRAM BEGINS-------------------------
    # ----------------------------------------------------------
    
    import os,sys
    sys.path.append(os.path.dirname(os.path.realpath(__file__))+'/_classes')
    from tome_v02 import resolution,normalize,msipoverlay
    from _Molecule import Molecule
    from bisect import bisect_left,bisect_right
    from _XLSX import XLSX
    
    xlfile = XLSX(specfile)
    
    exp = xlfile.pullspectrum(xlfile.wb.get_sheet_names()[0])[0] # pull first sheet in workbook
    
    sets = figuresetting(setting) # pull settings for specified preset
    
    if sovr != None: # if overrides are provided, change values in sets dictionary
        sys.stdout.write('Overriding preset figure size. New: width %f" height %f".\n' %(sovr[0],sovr[1]))
        sets['size'] = sovr
    if tovr != None:
        sys.stdout.write('Overriding preset output format. Using %s.\n' %tovr)
        sets['exten'] = tovr
    
    for species in simdict: # generate Molecule objects for each species
        simdict[species]['mf'] = Molecule(species)
        
    if mz == 'Auto': # automatically determine m/z range
        mz = [50000,0]
        for species in simdict:
            bounds = simdict[species]['mf'].bounds()
            if bounds[0] < mz[0]:
                mz[0] = bounds[0]-1
            if bounds[1] > mz[1]:
                mz[1] = bounds[1]+1
        sys.stdout.write('Automatically determining m/z window: %i - %i\n' %(int(mz[0]),int(mz[1])))
        sys.stdout.flush()
    
    snipleft,snipright = bisect_left(exp[0],mz[0]),bisect_right(exp[0],mz[1]) # find indicies of mz window in real spectrum
    exp = [exp[0][snipleft:snipright],exp[1][snipleft:snipright]] #trim to mz range
    
    try:
        expmax = max(exp[1]) #find max value in mz window
    except ValueError:
        sys.exit('\nThere does not appear to be any spectrum data in the area of the isotope pattern.\nThis usually occurs when the wrong isotope pattern was chosen.\nPlease check your input.')
    maxindex = exp[1].index(expmax) # index of max value
    
    sets['res'][1] = resolution(exp[0],exp[1],maxindex,expmax) #calculate resolution
    
    if sets['ymax'] == '100': # normalize spectrum to 100
        normalize(exp[1],100.)
    
    deltas = [[],[]] # list for  mass deltas
    for species in simdict: # calculate attributes of each species
        simdict[species]['mf'].res = sets['res'][1]
        if otype == 'bar':
            simdict[species]['x'] = simdict[species]['mf'].barip[0]
            simdict[species]['y'] = simdict[species]['mf'].barip[1]
        elif otype == 'gaussian':
            simdict[species]['mf'].gaussianisotopepattern()
            simdict[species]['x'] = simdict[species]['mf'].gausip[0]
            simdict[species]['y'] = simdict[species]['mf'].gausip[1]
        simdict[species]['em'] = simdict[species]['mf'].em # append exact mass of predicted species
        simdict[species]['max'] = max(simdict[species]['y']) # append max int to dictionary
        if sets['stats'] is True:
            simdict[species]['SER'] = simdict[species]['mf'].compare(exp)
        deltas[0].append(exp[0][maxindex]-simdict[species]['x'][simdict[species]['y'].index(max(simdict[species]['y']))])
        deltas[1].append(abs(deltas[0][-1]))
    
        if sets['norm'] == 'max' :#normalize simulations to local max
            lookwithin = 1 # mz range about isotope max to look within
            mzmax = simdict[species]['x'][simdict[species]['y'].index(simdict[species]['max'])] # index mz value of max intensity
            lind,rind = bisect_right(exp[0],mzmax-lookwithin),bisect_left(exp[0],mzmax+lookwithin) # find index of mz experimental values within given range
            locmax = max(exp[1][lind:rind]) # find the local maximum within that range
            simdict[species]['y'] = normalize(simdict[species]['y'],locmax) # normalize isotope pattern to the max in this region
        if sets['norm'] == 'top': # normalize simulations to top of chart
            if sets['ymax'] == '100':
                simdict[species]['y'] = normalize(simdict[species]['y'],100.)
            if sets['ymax'] == 'counts': 
                simdict[species]['y'] = normalize(simdict[species]['y'],expmax)
            if type(sets['ymax']) is int:
                simdict[species]['y'] = normalize(simdict[species]['y'],sets['ymax'])
        if type(sets['norm']) is int: # normalize simulations to specified value
            simdict[species]['y'] = normalize(simdict[species]['y'],sets['norm'])
    sets['delta'][1] = deltas[0][deltas[1].index(min(deltas[1]))] # set minimum mass delta as mass delta
    msipoverlay(exp,simdict,sets,mz,spectrum,expmax,otype=otype)


if __name__ == '__main__':
    ipoverlay(spectrum,simdict,mz,setting,sovr,tovr,otype,spectrum)
    import gc
    gc.collect()