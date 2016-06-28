"""
  isotope pattern check program v011 beta (dependant on tome)
  overlays a supplied predicted isotope pattern onto an acquired spectrum
  
new/changed:
    completely rewritten so that the actual plotting could be a standalone function
    this script now just handles the appropriate functions to accomplish a plotted mass spectrum with isotope pattern overlays
    tweaking of the figure is now handled by keyword arguments
    ---12.0 alpha

to add:
    reincorporate mass delta
    get rid of realmax necessity for figure generation (redefine in sets['ymax'])
    create a more generic version for naming the output file in the tome function
    databridge for masslynx conversion (is this even useful?)
"""
# provide the experimental spectrum xlsx
# the script will automatically use the first sheet
spectrum = '2016-06-13 09 Pd cation with boronic acid'

# supply simulated spectrum in dictionary format
# 'molecular formula':{'colour':(R,G,B) or hex,'alpha':float}
simdict = {
#'L2PdAr+OMe':{'colour':'#fc8d59','alpha':0.75},
#'L2PdAr+(2+)':{'colour':'#e41a1c','alpha':0.75},
#'(MeOC6H4BO)2OH':{'colour':'#e41a1c','alpha':0.75},
'L2PdPh':{'colour':'#4f81bd','alpha':0.75}
}

# choose a figure type for auto settings
# options: 'pub', 'pubsvg', 'inset', 'insetsvg', 'thesis', 'detailed'
# additional presets can be added in the presets() function below
setting = 'pub'

# preset settings can be overridden here (see presets() function for details)
override = {
#'mz': [1008,1028], # modify the m/z bounds of the figure
#'simtype': 'gaussian', # generate a gaussian spectrum
#'exten':'svg' # change to scalable vector graphic
#'size':[4,3] # change the size of the image
}

def presets(typ):
    """
    Returns a preset modification to the default settings in plotms()
    Supported keyword arguments:
    
    mz: the m/z bounds
        default: 'auto'
        specific bounds can be set in list or tuple form [m/z start,m/z end]
    
    outname: the name for the output file
        default: 'spectrum'
    
    simtype: the type of plot for isotope pattern overlays
        default: 'bar'
        'bar': a bar plot
        'gaussian': a predicted gaussian distribution centered around the bar intensities
    
    spectype: the type of spectrum being handed to the script
        default: 'continuum'
        'continuum' will connect each y value to the next
        'centroid' will plot a bar for each y value
    
    maxy: the maximum y axis value
        default: 'max'
        'max': based on the maximum intensity of the supplied spectrum
        value: a specified value
    
    norm: whether the script should normalize the supplied spectrum
        default: False
        True/False
    
    simnorm: how to normalize the isotope patterns
        default: 'spec'
        'spec': normalizes to the local maximum in the spectrum
        'top': normalizes to the maximum y value
        value: normalizes to that value
    
    axlabel: show axis labels
        default: [True,True]
        True/False for [x axis,y axis]
    
    axvalues: show axis values
        default: [True,True]
        True/False for [x axis,y axis]
    
    axshow: show axis
        default: [True,True]
        True/False for [x axis,y axis]
    
    offsetx: offsets the x axis to better show low-intensity species
        default: True
        True/False
    
    fs: desired fontsize
        default: 16
        integer value in points
    
    lw: linewidth for the plotted spectrum
        default: 1.5
        integer or float value
   
    bw: predicted pattern bar width
        default: 2
        integer multiple of the FWHM of the base peak in the isotope pattern
    
    axwidth: linewidth for the axes
        default: 1.5
        integer or float value (recommend setting the same as linewidth)
    
    simlabels: show simulation labels
        default: False
        True/False
    
    res: show calculated resolution
        default: False
        True/False
    
    specfont: specify font to be used
        default: 'Arial'
    
    size: figure size (in inches)
        default: [7.87,4.87]
        [width,height] (floating point values)
    
    dpiout: dots per inch for the output png
        default: 300
        integer
    
    exten: extension for the output file
        default: 'png'
        'png' for flat images
        'svg' vector image
    
    delta: calculate mass delta between predicted and experimental mass
        default: False
        True/False
    
    stats: calculate and output the standard error of regression between the experimental and predicted isotope patterns
        default: False
        True/False
    
    speccolour: the colour that the plotted spectrum should be
        default: 'k' (black)
        can be handed an (R,G,B) tuple or hex string
    
    padding: the subplot padding for the spectrum
        default: 'auto'
        can be handed [L,R,B,T] padding values
    
    Additional presets can be added in dictionary form with whatever setting changes are desired
    """
    ipsetting = {
    'pub':{'norm':True},
    'pubsvg':{'norm':True,'exten':'svg'},
    'inset':{'norm':True,'axlabels':[True,False],'axvalues':[True,False],'axshow':[True,False],'fs':12,'size':[2.8,2.8]},
    'insetsvg':{'norm':True,'axlabels':[True,False],'axvalues':[True,False],'axshow':[True,False],'fs':12,'size':[2.8,2.8],'exten':'svg'},
    'thesis':{'size':[6.5,4.0]},
    'detailed':{'fs':10,'simlabels':True,'res':True,'size':[5.,3.],'delta':True,'stats':True},
    }
    try:
        sys.stdout.write('Using figure preset "%s"\n' %(setting))
        sys.stdout.flush()
        return ipsetting[typ]
    except KeyError:
        raise KeyError('\nThe specified figure setting "%s" is not defined.\nPlease check your spelling' %setting)

if __name__ == '__main__':
    import os,sys
    if os.path.dirname(os.path.realpath(__file__))+'/_classes' not in sys.path:
        sys.path.append(os.path.dirname(os.path.realpath(__file__))+'/_classes')
    from _XLSX import XLSX
    from tome_v02 import plotms
    xlfile = XLSX(spectrum) # load excel file
    
    exp = xlfile.pullspectrum(xlfile.wb.get_sheet_names()[0])[0] # load spectrum from first sheet in workbook
    
    keywords = presets(setting) # pull preset
    keywords.update({'outname':xlfile.bookname[:-5]}) # set default output filename
    keywords.update(override) # apply any user overrides
    
    plotms(exp,simdict,**keywords)
    import gc
    gc.collect()