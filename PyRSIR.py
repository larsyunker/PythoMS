"""
 Processes a .raw masslynx MS file or previously processed excel file with the parameters provided
 
 PyRSIR (Python Reconstructed Single Ion Recording; previously SOAPy/PyRSIM)
 pronounced "piercer"

If you use this python script to convert to mzML, you should cite this paper
(of the folks who wrote the msconvert program)

Chambers, M.C. Nature Biotechnology 2012, 30, 918-920
doi: 10.1038/nbt.2377

        TO USE THIS PROGRAM
1) change the working directory to the directory containing the desired *.raw or *.mzML file
2) supply the file name of the *.raw or *.mzML file in quotations in the filename parameter below
3) create an excel file containing a sheet named "parameters" with the desired peaks outlied
   (see below for more details)
   WARNING! This program will keep all data values, but will remove any charts present in the file.
4) supply the file name of the *.xlsx file in quotation in the xlsx parameter below
5) set the number of scans to sum (any positive integer or list of positive integers) in the n parameter below
   (if no summing is desired, set n = 1)

The species names, start m/z, and end m/z values should be contained within a sheet named "parameters"
The first row is ignored for labelling convenience of the user.

Column #1: name of species (this will be the column heading in the output sheets)
    (not needed if the molecular formula is specified)
Column #2: molecular formula of this species 
    (not required, but can bypass the need for start and end values)
Column #3: what spectrum to find this species in (+,-,UV)
Column #4: start value (m/z or wavelength)
Column #5: end value (m/z or wavelength)
"""

# # input *.raw filename
# filename = 'ND-04082014-HIYAMA.raw'
#
# # Excel file to read from and output to (in *.xlsx format)
# xlsx = 'ND-04082014-HIYAMA.xlsx'
#
# # set number of scans to sum (integer or list of integers)
# n = [3]

import sys
import pylab as pl

from pythoms.tome import bindata
from pythoms.scripttime import ScriptTime
from pythoms.mzml import mzML
from pythoms.spectrum import Spectrum
from pythoms.molecule import IPMolecule
from pythoms.xlsx import XLSX


def pyrsir(
        filename,
        xlsx,
        n,
        plot=True,  # plot the data for a quick look
        verbose=True,  # chatty
        bounds_confidence=0.99,  #
        combine_spectra=True,  # whether or not to output a summed spectrum
        return_data=False,  #
):
    """
    A method for generating reconstructed single ion monitoring traces.

    :param filename: path to mzML or raw file to process
    :param xlsx: path to excel file with correctly formatted columns
    :param n: number of scans to sum together (for binning algorithm)
    :param plot: whether to plot and show the data for a quick look
    :param verbose: chatty mode
    :param bounds_confidence: confidence interval for automatically generated bounds (only applicable if molecular
        formulas are provided).
    :param combine_spectra: whether to output a summed spectrum
    :param return_data: whether to return data (if the data from the function is required by another function)
    :return:
    """
    def check_integer(val, name):
        """
        This function checks that the supplied values are integers greater than 1
        
        A integer value that is non-negative is required for the summing function.
        Please check your input value. 
        """
        if type(val) != list and type(val) != tuple:  # if only one value given for n
            val = [val]
        for num in val:
            if type(num) != int:
                sys.exit('\nThe %s value (%s) is not an integer.\n%s' % (name, str(num), check_integer.__doc__))
            if num < 1:
                sys.exit('\nThe %s value (%s) is less than 1.\n%s' % (name, str(num), check_integer.__doc__))
        return val

    def plots():
        """
        Function for generating a set of plots for rapid visual assessment of the supplied n-level
        Outputs all MS species with the same sum level onto the same plot
        requirements: pylab as pl
        """
        pl.clf()  # clears and closes old figure (if still open)
        pl.close()
        nplots = len(n) + 1

        # raw data
        pl.subplot(nplots, 1, 1)  # top plot

        for mode in mskeys:
            modekey = 'raw' + mode
            if modekey in rtime.keys():
                pl.plot(rtime[modekey], tic[modekey], linewidth=0.75, label='TIC')  # plot tic
                for key in sp:  # plot each species
                    if sp[key]['affin'] is mode:
                        pl.plot(rtime[modekey], sp[key]['raw'], linewidth=0.75, label=key)
        pl.title('Raw Data')
        pl.ylabel('Intensity')
        pl.tick_params(axis='x', labelbottom='off')

        # summed data
        loc = 2
        for num in n:
            pl.subplot(nplots, 1, loc)
            sumkey = str(num) + 'sum'
            for mode in mskeys:
                modekey = str(num) + 'sum' + mode
                if modekey in rtime.keys():
                    pl.plot(rtime[modekey], tic[modekey], linewidth=0.75, label='TIC')  # plot tic
                    for key in sp:
                        if sp[key]['affin'] is mode:  # if a MS species
                            pl.plot(rtime[modekey], sp[key][sumkey], linewidth=0.75, label=key)
            pl.title('Summed Data (n=%i)' % (num))
            pl.ylabel('Intensity')
            pl.tick_params(axis='x', labelbottom='off')
            loc += 1
        pl.tick_params(axis='x', labelbottom='on')
        pl.show()

    def output():
        """
        Writes the retrieved and calculated values to the excel workbook using the XLSX object
        """
        if newpeaks is True:  # looks for and deletes any sheets where the data will be changed
            if verbose is True:
                sys.stdout.write('Clearing duplicate XLSX sheets.')
            delete = []
            for key in newsp:  # generate strings to look for in excel file
                delete.append('Raw Data (' + sp[key]['affin'] + ')')
                for num in n:
                    delete.append(str(num) + ' Sum (' + sp[key]['affin'] + ')')
                    delete.append(str(num) + ' Normalized (' + sp[key]['affin'] + ')')
            delete.append('Isotope Patterns')
            xlfile.removesheets(delete)  # remove those sheets
            if verbose is True:
                sys.stdout.write(' DONE.\n')

        if verbose is True:
            sys.stdout.write('Writing to "%s"' % xlfile.bookname)
            sys.stdout.flush()

        for mode in mskeys:  # write raw data to sheets
            modekey = 'raw' + mode
            if modekey in rtime.keys():
                sheetname = 'Raw Data (' + mode + ')'
                xlfile.writersim(sp, rtime[modekey], 'raw', sheetname, mode, tic[modekey])

        for num in n:  # write summed and normalized data to sheets
            sumkey = str(num) + 'sum'
            normkey = str(num) + 'norm'
            for mode in mskeys:
                modekey = 'raw' + mode
                if modekey in rtime.keys():
                    if max(n) > 1:  # if data were summed
                        sheetname = str(num) + ' Sum (' + mode + ')'
                        xlfile.writersim(sp, rtime[sumkey + mode], sumkey, sheetname, mode,
                                         tic[sumkey + mode])  # write summed data
                    sheetname = str(num) + ' Normalized (' + mode + ')'
                    xlfile.writersim(sp, rtime[sumkey + mode], normkey, sheetname, mode)  # write normalized data

        for key, val in sorted(sp.items()):  # write isotope patterns
            if sp[key]['affin'] in mskeys:
                xlfile.writemultispectrum(
                    sp[key]['spectrum'][0],  # x values
                    sp[key]['spectrum'][1],  # y values
                    key,  # name of the spectrum
                    xunit='m/z',  # x unit
                    yunit='Intensity (counts)',  # y unit
                    sheetname='Isotope Patterns',  # sheet name
                    chart=True,  # output excel chart
                )

        if rd is None:
            for key, val in sorted(chroms.items()):  # write chromatograms
                xlfile.writemultispectrum(chroms[key]['x'], chroms[key]['y'], chroms[key]['xunit'],
                                          chroms[key]['yunit'], 'Function Chromatograms', key)

        uvstuff = False
        for key in sp:  # check for UV-Vis spectra
            if sp[key]['affin'] is 'UV':
                uvstuff = True
                break
        if uvstuff is True:
            for ind, val in enumerate(tic['rawUV']):  # normalize the UV intensities
                tic['rawUV'][ind] = val / 1000000.
            xlfile.writersim(sp, rtime['rawUV'], 'raw', 'UV-Vis', 'UV', tic['rawUV'])  # write UV-Vis data to sheet

        if sum_spectra is not None:  # write all summed spectra
            for fn in sum_spectra:
                specname = '%s %s' % (mzml.functions[fn]['mode'], mzml.functions[fn]['level'])
                if 'target' in mzml.functions[fn]:
                    specname += ' %.3f' % mzml.functions[fn]['target']
                specname += ' (%.3f-%.3f)' % (mzml.functions[fn]['window'][0], mzml.functions[fn]['window'][1])
                xlfile.writemultispectrum(
                    sum_spectra[fn][0],  # x values
                    sum_spectra[fn][1],  # y values
                    specname,  # name of the spectrum
                    xunit='m/z',  # x unit
                    yunit='Intensity (counts)',  # y unit
                    sheetname='Summed Spectra',  # sheet name
                    chart=True,  # output excel chart
                )

        if verbose is True:
            sys.stdout.write(' DONE\n')

    def prepformula(dct):
        """looks for formulas in a dictionary and prepares them for pullspeciesdata"""
        for species in dct:
            if 'affin' not in dct[species]:  # set affinity if not specified
                fn = dct[species]['function']
                if mzml.functions[fn]['type'] == 'MS':
                    dct[species]['affin'] = mzml.functions[fn]['mode']
                if mzml.functions[fn]['type'] == 'UV':
                    dct[species]['affin'] = 'UV'
            if 'formula' in dct[species]and dct[species]['formula'] is not None:
                try:
                    dct[species]['mol'].res = res  # sets resolution in Molecule object
                except NameError:
                    res = int(mzml.auto_resolution())
                    dct[species]['mol'].res = res
                # dct[species]['mol'].sigma = dct[species]['mol'].sigmafwhm()[1]  # recalculates sigma with new resolution
                dct[species]['bounds'] = dct[species]['mol'].bounds  # caclulates bounds
        return dct

    # ----------------------------------------------------------
    # -------------------PROGRAM BEGINS-------------------------
    # ----------------------------------------------------------

    if verbose is True:
        stime = ScriptTime()
        stime.printstart()

    n = check_integer(n, 'number of scans to sum')  # checks integer input and converts to list

    if type(xlsx) != dict:
        if verbose is True:
            sys.stdout.write('Loading processing parameters from excel file')
            sys.stdout.flush()
        xlfile = XLSX(xlsx, verbose=verbose)
        sp = xlfile.pullrsimparams()
    else:  # if parameters were provided in place of an excel file
        sp = xlsx

    mskeys = ['+', '-']
    for key in sp:
        if 'formula' in sp[key] and sp[key]['formula'] is not None:  # if formula is specified
            sp[key]['mol'] = IPMolecule(sp[key]['formula'])  # create Molecule object
            sp[key]['bounds'] = sp[key]['mol'].calculate_bounds(
                bounds_confidence
            )  # generate bounds from molecule object with this confidence interval
    if verbose is True:
        sys.stdout.write(' DONE\n')

    rtime = {}  # empty dictionaries for time and tic
    tic = {}
    rd = False
    for mode in mskeys:  # look for existing positive and negative mode raw data
        try:
            modedata, modetime, modetic = xlfile.pullrsim('Raw Data (' + mode + ')')
        except KeyError:
            continue
        except UnboundLocalError:  # catch for if pyrsir was not handed an excel file
            continue
        if verbose is True:
            sys.stdout.write('Existing (%s) mode raw data were found, grabbing those values.' % mode)
            sys.stdout.flush()
        rd = True  # bool that rd is present
        modekey = 'raw' + mode
        sp.update(modedata)  # update sp dictionary with raw data
        for key in modedata:  # check for affinities
            if 'affin' not in sp[key]:
                sp[key]['affin'] = mode
        rtime[modekey] = list(modetime)  # update time list
        tic[modekey] = list(modetic)  # update tic list
        if verbose is True:
            sys.stdout.write(' DONE\n')

    # sp = prepformula(sp)
    newpeaks = False
    if rd is True:
        newsp = {}
        sum_spectra = None
        for key in sp:  # checks whether there is a MS species that does not have raw data
            if 'raw' not in sp[key]:
                newsp[key] = sp[key]  # create references in the namespace
        if len(newsp) is not 0:
            newpeaks = True
            if verbose is True:
                sys.stdout.write('Some peaks are not in the raw data, extracting these from raw file.\n')
            ips = xlfile.pullmultispectrum(
                'Isotope Patterns')  # pull predefined isotope patterns and add them to species
            for species in ips:  # set spectrum list
                sp[species]['spectrum'] = [ips[species]['x'], ips[species]['y']]
            mzml = mzML(filename)  # load mzML class
            sp = prepformula(sp)  # prep formula etc for summing
            newsp = prepformula(newsp)  # prep formula species for summing
            for species in newsp:
                if 'spectrum' not in newsp[species]:
                    newsp[species]['spectrum'] = Spectrum(3, newsp[species]['bounds'][0], newsp[species]['bounds'][1])
            newsp = mzml.pull_species_data(newsp)  # pull data
        else:
            if verbose is True:
                sys.stdout.write('No new peaks were specified. Proceeding directly to summing and normalization.\n')

    if rd is False:  # if no raw data is present, process mzML file
        mzml = mzML(filename, verbose=verbose)  # load mzML class
        sp = prepformula(sp)
        sp, sum_spectra = mzml.pull_species_data(sp, combine_spectra)  # pull relevant data from mzML
        chroms = mzml.pull_chromatograms()  # pull chromatograms from mzML
        rtime = {}
        tic = {}
        for key in sp:  # compare predicted isotope patterns to the real spectrum and save standard error of the regression
            func = sp[key]['function']
            if mzml.functions[func]['type'] == 'MS':  # determine mode key
                if combine_spectra is True:
                    sp[key]['spectrum'] = sum_spectra[sp[key]['function']].trim(
                        xbounds=sp[key]['bounds'])  # extract the spectrum object
                mode = 'raw' + mzml.functions[func]['mode']
            if mzml.functions[func]['type'] == 'UV':
                mode = 'rawUV'
            if mode not in rtime:  # if rtime and tic have not been pulled from that function
                rtime[mode] = mzml.functions[func]['timepoints']
                tic[mode] = mzml.functions[func]['tic']
            if 'formula' in sp[key] and sp[key]['formula'] is not None:
                sp[key]['match'] = sp[key]['mol'].compare(sp[key]['spectrum'])
        if combine_spectra is True:
            for fn in sum_spectra:
                sum_spectra[fn] = sum_spectra[fn].trim()  # convert Spectrum objects into x,y lists

    # if max(n) > 1: # run combine functions if n > 1
    for num in n:  # for each n to sum
        if verbose is True:
            sys.stdout.write('\r%d Summing species traces.' % num)
        sumkey = str(num) + 'sum'
        for key in sp:  # bin each species
            if sp[key]['affin'] in mskeys or mzml.functions[sp[key]['function']][
                'type'] == 'MS':  # if species is MS related
                sp[key][sumkey] = bindata(num, sp[key]['raw'])
        for mode in mskeys:
            sumkey = str(num) + 'sum' + mode
            modekey = 'raw' + mode
            if modekey in rtime.keys():  # if there is data for that mode
                rtime[sumkey] = bindata(num, rtime[modekey], num)
                tic[sumkey] = bindata(num, tic[modekey])
    if verbose is True:
        sys.stdout.write(' DONE\n')
        sys.stdout.flush()
    # else:
    #    for key in sp: # create key for normalization
    #        sp[key]['1sum'] = sp[key]['raw']

    for num in n:  # normalize each peak's chromatogram
        if verbose is True:
            sys.stdout.write('\r%d Normalizing species traces.' % num)
            sys.stdout.flush()
        sumkey = str(num) + 'sum'
        normkey = str(num) + 'norm'
        for mode in mskeys:
            modekey = 'raw' + mode
            if modekey in rtime.keys():  # if there is data for that mode
                for key in sp:  # for each species
                    if sp[key]['affin'] in mskeys or mzml.functions[sp[key]['function']][
                        'type'] == 'MS':  # if species has affinity
                        sp[key][normkey] = []
                        for ind, val in enumerate(sp[key][sumkey]):
                            # sp[key][normkey].append(val/(mzml.function[func]['tic'][ind]+0.01)) #+0.01 to avoid div/0 errors
                            sp[key][normkey].append(
                                val / (tic[sumkey + sp[key]['affin']][ind] + 0.01))  # +0.01 to avoid div/0 errors
    if verbose is True:
        sys.stdout.write(' DONE\n')

    if return_data is True:  # if data is to be used by another function, return the calculated data
        return mzml, sp, rtime, tic, chroms

    # import pickle #pickle objects (for troubleshooting)
    # pickle.dump(rtime,open("rtime.p","wb"))
    # pickle.dump(tic,open("tic.p","wb"))
    # pickle.dump(chroms,open("chroms.p","wb"))
    # pickle.dump(sp,open("sp.p","wb"))

    output()  # write data to excel file

    if verbose is True:
        sys.stdout.write('\rUpdating paramters')
        sys.stdout.flush()
    xlfile.updatersimparams(sp)  # update summing parameters
    if verbose is True:
        sys.stdout.write(' DONE\n')

    if verbose is True:
        sys.stdout.write('\rSaving "%s" (this may take some time)' % xlfile.bookname)
        sys.stdout.flush()
    xlfile.save()
    if verbose is True:
        sys.stdout.write(' DONE\n')

    if verbose is True:
        if verbose is True:
            sys.stdout.write('Plotting traces')
        plots()  # plots for quick review
        if verbose is True:
            sys.stdout.write(' DONE\n')
    if verbose is True:
        stime.printelapsed()


# import sys
#
# if len(sys.argv) > 1:  # if script was initiated from the command line, pull parameters from there
#     try:
#         pyrsir(sys.argv[1], sys.argv[2], sys.arg[3])
#     except IndexError:
#         raise IOError(
#             'The pyrsim function requires three inputs:\n- The raw filename\n- The excel parameters file\n- The number of scans to sum')

if __name__ == '__main__':
    curdir = input('Directory: ')  # cur
    filename = input('Raw or mzML filename: ')  # input filename
    xlsx = input('Excel file name: ')  # excel file name with summing parameters
    try:
        n = [int(i) for i in input('Scans to bin (separate multiple values with commas): ').split(',')]  # scans to bin
    except ValueError:
        print('Input value error, applying no binning')
        n = [1]
    pyrsir(filename, xlsx, n)
    sys.stdout.write('fin.')
    sys.stdout.flush()
    import gc

    gc.collect()
