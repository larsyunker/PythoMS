"""
  isotope pattern check program v011 beta (dependant on tome)
  overlays a supplied predicted isotope pattern onto an acquired spectrum
  
CHANGELOG:
- added databridge to extract spectra directly from a raw or mzml.gz file (validated)
---12.3
"""
# provide the experimental spectrum xlsx
# the script will automatically use the first sheet
spectrum = 'L2PdArCl'
# number of lines to skip in the excel file
# (e.g. if there are spectrum details above the actual spectrum values)
skiplines = 0

# sheet name in the excel file (if this is not specified, the script will use the first sheet in the file)
# sheetname = 'JW-RH-08-203 +'

# provide species to be simulated in dictionary format
# 'molecular formula':{'colour': ... ,'alpha':0-1}
# colour can be (R,G,B), (C,M,Y,K), or 'hex'
simdict = {
    # 'Ar+I':{'colour':'1f78b4','alpha':0.5},
    # 'L2PdAr+I':{'colour':'#6a3d9a','alpha':0.5},
    'L2PdAr+Cl': {'colour': '#006838', 'alpha': 0.5},
    # 'L2PdAr+(2+)':{'colour':'#e41a1c','alpha':0.5},
    # 'Pd(PPh2)Ar+':{'colour':'696969','alpha':0.5},
    # 'L2PdAr+C6H4CH3':{'colour':'#ff7f00','alpha':0.5},
    # 'L2PdAr+IMeOH':{'colour':(146,102,194),'alpha':0.5},
    # '(Ar+I)2PF6':{'colour':'#1f78b4','alpha':0.5},
    # 'L2PdAr+CH3C6H4MeOH':{'colour':'#ff7f00','alpha':0.5}
    # 'L2PdAr+OH':{'colour':'66c2a5','alpha':0.5},
    # '[NBu4]3PF6I':{'colour':'k','alpha':0.5},
    # 'L2Pd2(Ar+)2(OH)2(2+)':{'colour':'b','alpha':0.5},
    # 'L2Pd2(Ar+)2OHOMe(2+)':{'colour':'r','alpha':0.5},
    # 'L2Pd2(Ar+)2(OMe)2(2+)':{'colour':'g','alpha':0.5},
    # 'L2PdAr+OMe':{'colour':'b','alpha':0.5},
    # 'C54H57O5P2Ru':{'colour':'b','alpha':0.5},
    # 'C27H33N2O4':{'colour':'696969','alpha':0.5},
    # 'TiCp2MeCNOMe':{'colour':'b','alpha':0.5}
}

# choose a figure type for auto settings
# options: 'pub', 'pubsvg', 'inset', 'insetsvg', 'thesis', 'detailed'
# additional presets can be added in the presets() function below
setting = 'pub'

# preset settings can be overridden here (see presets() function for details)
override = {
    # 'norm':False,
    # 'bw':0.1, # bar width (in units of m/z)
    # 'exten':'svg', # change to scalable vector graphic
    # 'mz': [1015,1025], # modify the m/z bounds of the figure
    # 'offsetx':False, # apply a slight offset to the x axis
    # 'simtype': 'gaussian', # generate a gaussian spectrum
    # 'size':[10,4], # change the size of the image [width,height]
    # 'specfont':'Calibri', # change the font of the labels
    'stats': False,
    # 'xlabel':False,
    'ylabel': False,
    # 'yvalues':False,
    # 'showy':False,
    # 'norm':False,
    # 'maxy':120,
    # 'simnorm':120,
    # 'spectype':'centroid',
    # 'normwindow':1.,
    # 'fs':8,
    # 'output':'show',
}


def presets(typ):
    """
    Returns a preset modification to the default settings in plotms()
    Supported keyword arguments:
    
    axwidth: linewidth for the axes
        default: 1.5
        integer or float value (recommend setting the same as linewidth)
    
    bw: predicted pattern bar width
        default: 'auto'
        'auto' makes the bars 2* the full width at half max
        if a value is supplied, the bars will be that wide (in units of m/z)
    
    delta: show calculated mass delta between predicted and experimental mass
        default: False
        True/False
    
    dpiout: dots per inch for the output png
        default: 300
        integer
    
    exten: extension for the output file
        default: 'png'
        'png' for flat images
        'svg' vector image
    
    fs: desired fontsize
        default: 16
        integer value in points
    
    maxy: the maximum y axis value
        default: 'max'
        'max': based on the maximum intensity of the supplied spectrum
        value: a specified value
    
    lw: linewidth for the plotted spectrum
        default: 1.5
        integer or float value
    
    mz: the m/z bounds
        default: 'auto'
        specific bounds can be set in list or tuple form [m/z start,m/z end]
    
    norm: whether the script should normalize the supplied spectrum
        default: True
        True/False
    
    normwindow: what m/z window width should the autonormalization look
        default 'fwhm' (look within the full width at half max)
        otherwise, hand this an m/z value
    
    offsetx: offsets the x axis to better show low-intensity species
        default: True
        True/False
    
    outname: the name for the output file
        default: 'spectrum'
    
    padding: the subplot padding for the spectrum
        default: 'auto'
        can be handed [L,R,B,T] padding values
    
    res: show calculated resolution
        default: False
        True/False
    
    showy: show y axis line
        default: True
    
    showx: show x axis line
        default: True
    
    simlabels: show simulation labels
        default: False
        True/False
    
    simnorm: how to normalize the isotope patterns
        default: 'spec'
        'spec': normalizes to the local maximum in the spectrum
        'top': normalizes to the maximum y value
        value: normalizes to that value
    
    simtype: the type of plot for isotope pattern overlays
        default: 'bar'
        'bar': a bar plot
        'gaussian': a predicted gaussian distribution centered around the bar intensities
    
    size: figure size (in inches)
        default: [7.87,4.87]
        [width,height] (floating point values)
    
    speccolour: the colour that the plotted spectrum should be
        default: 'k' (black)
        can be handed an (R,G,B) or (C,M,Y,K) tuple or hex string
    
    specfont: specify font to be used
        default: 'Arial'
    
    spectype: the type of spectrum being handed to the script
        default: 'continuum'
        'continuum' will connect each y value to the next
        'centroid' will plot a bar for each y value
    
    stats: show the standard error of regression between the experimental and predicted isotope patterns
        default: False
        True/False
    
    verbose: chatty
        bool
    
    xlabel: show x label
        default: True
    
    xvalues: show x values
        default: True
    
    ylabel: show y label
        default: True
    
    yvalues: show y values
        default: True
    
    
    Additional presets can be added in dictionary form with whatever setting changes are desired
    """
    ipsetting = {
        'pub': {},
        'pubsvg': {'exten': 'svg'},
        'inset': {'ylabel': False, 'yvalues': False, 'showy': False, 'fs': 12, 'size': [2.8, 2.8]},
        'insetsvg': {'ylabel': False, 'yvalues': False, 'showy': False, 'fs': 12, 'size': [2.8, 2.8], 'exten': 'svg'},
        'thesis': {'size': [6.5, 4.0]},
        'detailed': {'norm': False, 'fs': 10, 'simlabels': True, 'res': True, 'size': [6.5, 4.0], 'delta': True,
                     'stats': True},
    }
    try:
        sys.stdout.write('Using figure preset "%s"\n' % (setting))
        sys.stdout.flush()
        return ipsetting[typ]
    except KeyError:
        raise KeyError('\nThe specified figure setting "%s" is not defined.\nPlease check your spelling' % setting)


if __name__ == '__main__':
    import sys
    from PythoMS._classes._XLSX import XLSX
    from tome_v02 import plotms

    keywords = presets(setting)  # pull preset kwargs

    if spectrum.lower().endswith('.mzml.gz') or spectrum.lower().endswith('.raw'):  # if supplied with a mass spec file
        from PythoMS._classes._mzML import mzML

        mzml = mzML(spectrum)
        exp = mzml.sum_scans()
        keywords.update({'outname': mzml.filename.split('.')[0]})  # set default output filename

    else:  # otherwise assume that it is an excel file
        from PythoMS._classes._XLSX import XLSX

        xlfile = XLSX(spectrum, verbose=True)  # load excel file
        try:  # if sheet name was specified
            sname = sheetname
        except NameError:  # otherwise use the first sheet
            sname = xlfile.wb.get_sheet_names()[0]
        exp = xlfile.pullspectrum(sname, skiplines=skiplines)[0]  # load spectrum from first sheet in workbook
        keywords.update({'outname': xlfile.bookname[:-5] + ' (' + sname + ')'})  # set default output filename

    keywords.update(override)  # apply any user overrides
    plotms(exp, simdict, **keywords)
    import gc

    gc.collect()
