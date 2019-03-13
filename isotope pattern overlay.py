import sys
import os
from pythoms.xlsx import XLSX
from pythoms.tome import plot_mass_spectrum
from pythoms.mzml import mzML

"""
Overlays a predicted isotope pattern over an acquired spectrum. 

To use this script: 
1. specify a file path to the experimental spectrum (XLSX, RAW, or mzML)
    a. specify a sheet name if applicable (otherwise the first sheet will be used in the excel file)
    b. specify the number of header lines (if not applicable, set to 0)
2. Define the species to be overlain on the spectrum. These are specified in the format 
    'MOLECULARFORMULA: {'colour': COLOUR, 'alpha': ALPHA},
3. Choose a preset (see the preset method below for options and details)
4. Specify overrides which deviate from the preset as desired. 
5. Run the script (or execute in console). 

"""
# current directory
curdir = 'C:\\Temp'

# provide the experimental spectrum xlsx or RAW/mzml file
# the script will automatically use the first sheet
spectrum = 'LY-2014-06-12 14.mzML.gz'
# number of lines to skip in the excel file
# (e.g. if there are spectrum details above the actual spectrum values)
skiplines = 0
column_offset = 1

# sheet name in the excel file (if this is not specified, the script will use the first sheet in the file)
sheetname = None

# provide species to be simulated in dictionary format
# 'molecular formula':{'colour': ... ,'alpha':0-1}
# colour can be (R,G,B), (C,M,Y,K), or 'hex'
simdict = {
    'Ar+I': {'colour': '1f78b4', 'alpha': 0.5},
    'L2PdAr+I': {'colour': '#6a3d9a', 'alpha': 0.5},
}

# choose a figure type for auto settings
# options: 'pub', 'pubsvg', 'inset', 'insetsvg', 'thesis', 'detailed'
# additional presets can be added in the presets() function below
setting = 'pub'

# preset settings can be overridden here (see presets() function for details)
override = {
    # 'norm': False,
    # 'bw': 0.1,  # bar width (in units of m/z)
    # 'exten': 'svg',  # change to scalable vector graphic
    # 'mz': [1015, 1025],  # modify the m/z bounds of the figure
    # 'offsetx':False, # apply a slight offset to the x axis
    # 'simtype': 'gaussian', # generate a gaussian spectrum
    # 'size':[10,4],  # change the size of the image [width,height]
    # 'specfont': 'Calibri', # change the font of the labels
    'stats': False,
    'resolution': 5000,
    # 'xlabel': False,
    # 'ylabel': False,
    # 'yvalues': False,
    # 'showy': False,
    # 'norm': False,
    # 'maxy': 120,
    # 'simnorm': 100,
    # 'spectype': 'centroid',
    # 'normwindow': 1.,
    # 'fs': 8,
    # 'output':'show',
    'ipmol_kwargs': {'dropmethod': 'threshold'},
}


def presets(typ):
    """
    Returns a preset modification to the default settings in plotms()
    Supported keyword arguments:

    :param list realspec: A paired list of x and y values of the form ``[[x values],[y values]]``
    :param dict simdict: This can either be a molecular formula to predict the isotope pattern of (string),
        a list of formulae, or a dictionary of the form
        ``simdict = {'formula1':{'colour':<hex or name or RGB tuple>, 'alpha':float}, ...}``.
        If this is dictionary is left empty, no isotope patterns will be overlaid on the output
        spectrum.
    :param list mz: The *m/z* bounds for the output spectrum. Default: 'auto', but can be supplied
        with a tuple or list of length 2 of the form ``[x start, x end]``.
    :param str outname: Name of the file to be saved.
    :param str output: Save ('save') or show ('show') the figure.
    :param str simtype: The type for the isotope pattern simulation overlay. Options: 'bar' or 'gaussian'.
    :param str spectype: The type of spectrum being handed to the function. Options: 'continuum' or 'centroid'.
    :param float maxy: The maximum y value for the spectrum. Options: 'max' or specify a value
    :param bool norm: Normalize the spectrum. Options: bool
    :param str, float simnorm: Normalize the isotope pattern simulations to what value. Options: 'top', 'spec', or
        specify a value. Top will normalize the patterns to ``maxy``, and will only function if maxy is not 'max'.
        Spec will normalize the patterns to the maximum spectrum y value within the x bounds of the
        simulated pattern.
        Specifying a value will normalize all isotope patterns to that value.
    :param bool xlabel: Whether to show the label for the *m/z* axis.
    :param bool ylabel: Whether to show the y-axis label.
    :param bool xvalues: Whether to show the values of the x-axis.
    :param bool yvalues: Whether to show the values of the y-axis.
    :param bool showx: Whether to show the x-axis line.
    :param bool showy: Whether to show the y-axis line.
    :param bool offsetx: Whether to offset the x-axis slightly.
        Enabling this shows makes it easier to see low intensity peaks.
    :param int fs: Font size to use for labels.
    :param float lw: Line width for the plotted spectrum.
    :param float axwidth: Line width for the axes and tick marks. Default 1.5
    :param bool simlabels: Whether to show the names of the simulated isotope patterns.
        The names will be exactly as supplied in ``simdict``.
    :param float bw: The width of the bar in *m/z* for bar isotope patterns. Options: 'auto' or float.
        This only has an affect if *simtype* is 'bar'.
        Auto make the bars equal to 2 times the full width at half max of the peak they are simulating.
    :param str specfont: The font to use for text in the plot. The specified font must be accepted by matplotlib.
    :param list size: The size in inches for the output figure. This must be a list of length 2 of the form
        ``[width,height]``.
    :param int dpiout: The dots per inch for the output figure.
    :param str exten: The file extension for the output figure. Options: 'png', 'svg', or other supported by matplotlib.
    :param float resolution: Override the auto-resolution calculation with a specified instrument resolution
    :param bool res_label: Whether to output the resolution of the spectrum onto the figure.
    :param bool delta: Whether to calculate and output the mass delta between the exact mass predicted by the isotope
        pattern simulation and the location of the maximum intensity within the bounds specified by *normwindow*.
    :param bool stats: Whether to calculate and output the goodness of fit between the predicted isotope pattern and
        the supplied spectrum. This functionality is still a work in progress.
    :param speccolour: The colour for the real spectrum , # colour for the spectrum to be plotted
    :param list padding: This allows the user to specify the subplot padding of the output figure.
        Options: 'auto' or list of the form ``[left,right,bottom,top]`` scalars.
    :param bool verbose: Verbose option for the script. Options: bool.
    :param float normwindow: The *m/z* window width within with too look for a maximum intensity value.
        This will only have an effect if *delta* is ``True``.
        Options: 'fwhm' for full width at half max or float.
    :param dict annotations: Annotations for the spectrum in dictionary form: ``{'thing to print':[x,y],}``.
    :param normrel: The maximum value for normalization. This can be used to globally set the top value for normalizing
        simulated isotope patterns. This is used most often to show the lack of an isotope pattern in the shown area.
    :param ipmol_kwargs: Keyword arguments to use for IPMolecule calls. See IPMolecule for more details.
    :param kwargs: catch for unused kwargs
    
    
    Additional presets can be added in dictionary form with whatever setting changes are desired
    """
    ipsetting = {
        'pub': {},
        'pubsvg': {  # formatted for publishing
            'exten': 'svg',
        },
        'inset': {  # sized for an inset in another figure
            'ylabel': False,
            'yvalues': False,
            'showy': False,
            'fs': 12,
            'size': [2.8, 2.8],
        },
        'insetsvg': {  # inset in svg format
            'ylabel': False,
            'yvalues': False,
            'showy': False,
            'fs': 12,
            'size': [2.8, 2.8],
            'exten': 'svg',
        },
        'thesis': {  # sized to fit the width of a standard word document
            'size': [6.5, 4.0],
        },
        'detailed': {  # for detailed analysis with additional stats outputs
            'norm': False,
            'fs': 10,
            'simlabels': True,
            'res_label': True,
            'size': [6.5, 4.0],
            'delta': True,
            'stats': True,
        },
    }
    try:
        sys.stdout.write('Using figure preset "%s"\n' % (setting))
        sys.stdout.flush()
        return ipsetting[typ]
    except KeyError:
        raise KeyError('\nThe specified figure setting "%s" is not defined.\nPlease check your spelling' % setting)


if __name__ == '__main__':
    os.chdir(curdir)  # change to current working directory

    keywords = presets(setting)  # pull preset kwargs

    if spectrum.lower().endswith('.mzml.gz') or spectrum.lower().endswith('.raw'):  # if supplied with a mass spec file
        mzml = mzML(spectrum)
        exp = mzml.sum_scans()
        keywords.update({'outname': mzml.filename.split('.')[0]})  # set default output filename

    else:  # otherwise assume that it is an excel file
        xlfile = XLSX(spectrum, verbose=True)  # load excel file
        if sheetname is None:  # otherwise use the first sheet
            sheetname = xlfile.wb.sheetnames[0]
        exp = xlfile.pullspectrum(sheetname, skiplines=skiplines)[0]  # load spectrum from first sheet in workbook
        keywords.update({  # set default output filename
            'outname': f'{xlfile.bookname[:-5]} ({sheetname})',
        })

    keywords.update(override)  # apply any user overrides
    plot_mass_spectrum(exp, simdict, **keywords)
    import gc

    gc.collect()
