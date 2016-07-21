"""
 Processes a .raw masslynx MS file or previously processed excel file with the parameters provided
 
 PyRSIR (Python Reconstructed Single Ion Recording; previously SOAPy/PyRSIM)
 pronounced "piercer"
 
 version 027
 new:
    updated to work with mzML v 2.4
    ---27.6 incompatible with mzML v2.4 or greater

to add/fix:
    update mzml to work with calibration
    create functionality for per-peak summing (create daughter dictionary?)
        if peaks overlap, combine
    
    fix integ() in mzML to work if there are no data points at the edge of the spectrum
    add output for images of isotope patterns and embed in excel (?possible?)
    ----

If you use this python script to process data, you should cite this paper
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

# input *.raw filename
filename = 'LY-2106-07-19 04'

# Excel file to read from and output to (in *.xlsx format)
xlsx = '2016-07-19 04 AgNO3 in sequential (unusable)'

# set number of scans to sum (integer or list of integers)
n = [3,5]

# ----------------------------------------------------------
# -------------------FUNCTION DEFINITIONS-------------------
# ----------------------------------------------------------

def pyrsir(filename,xlsx,n):    
    def checkinteger(val,name):
        """
        This function checks that the supplied values are integers greater than 1
        
        A integer value that is non-negative is required for the summing function.
        Please check your input value. 
        """
        import sys
        if type(val) != list and type(val) != tuple: # if only one value given for n
            val = [val]
        for num in val:
            if type(num) != int:
                sys.exit('\nThe %s value (%s) is not an integer.\n%s' %(name,str(num),checkinteger.__doc__))
            if num < 1:
                sys.exit('\nThe %s value (%s) is less than 1.\n%s' %(name,str(num),checkinteger.__doc__))
        return val
    
    def plots():
        """
        Function for generating a set of plots for rapid visual assessment of the supplied n-level
        Outputs all MS species with the same sum level onto the same plot
        requirements: pylab as pl
        """
        import pylab as pl
        pl.clf() # clears and closes old figure (if still open)
        pl.close()
        nplots = len(n)+1
        
        # raw data
        pl.subplot(nplots,1,1) # top plot
        
        for mode in mskeys:
            modekey = 'raw'+mode
            if modekey in rtime.keys():
                pl.plot(rtime[modekey],tic[modekey], linewidth = 0.75, label = 'TIC') #plot tic
                for key in sp: # plot each species
                    if sp[key]['affin'] is mode:
                        pl.plot(rtime[modekey],sp[key]['raw'], linewidth=0.75, label=key)
        pl.title('Raw Data')
        pl.ylabel('Intensity')
        pl.tick_params(axis='x',labelbottom='off')
        
        # summed data
        loc = 2
        for num in n:
            pl.subplot(nplots,1,loc)
            sumkey = str(num)+'sum'
            for mode in mskeys:
                modekey = str(num)+'sum'+mode
                if modekey in rtime.keys():
                    pl.plot(rtime[modekey],tic[modekey], linewidth = 0.75, label = 'TIC') #plot tic
                    for key in sp:
                        if sp[key]['affin'] is mode: #if a MS species
                            pl.plot(rtime[modekey],sp[key][sumkey], linewidth=0.75, label=key)
            pl.title('Summed Data (n=%i)' %(num))
            pl.ylabel('Intensity')
            pl.tick_params(axis='x',labelbottom='off')
            loc+=1
        pl.tick_params(axis='x',labelbottom='on')
        pl.show()
  
    def output():
        """
        Writes the retrieved and calculated values to the excel workbook using the XLSX object
        """
        if newpeaks is True: # looks for and deletes any sheets where the data will be changed
            sys.stdout.write('Clearing duplicate XLSX sheets.')
            delete = []
            for key in newsp: # generate strings to look for in excel file
                delete.append('Raw Data ('+sp[key]['affin']+')')
                for num in n:
                    delete.append(str(num)+' Sum ('+sp[key]['affin']+')')
                    delete.append(str(num)+' Normalized ('+sp[key]['affin']+')')
            delete.append('Isotope Patterns')
            xlfile.removesheets(delete) # remove those sheets
            sys.stdout.write(' DONE.\n')
        
        sys.stdout.write('Writing to "%s"' %xlfile.bookname)
        sys.stdout.flush()
                
        for mode in mskeys: # write raw data to sheets
            modekey = 'raw'+mode
            if modekey in rtime.keys():
                sheetname = 'Raw Data ('+mode+')'
                xlfile.writersim(sp,rtime[modekey],'raw',sheetname,mode,tic[modekey])

        for num in n: # write summed and normalized data to sheets
            sumkey = str(num)+'sum'
            normkey = str(num)+'norm'
            for mode in mskeys:
                modekey = 'raw'+mode
                if modekey in rtime.keys():
                    if max(n) > 1: # if data were summed
                        sheetname = str(num)+' Sum ('+mode+')'
                        xlfile.writersim(sp,rtime[sumkey+mode],sumkey,sheetname,mode,tic[sumkey+mode]) # write summed data
                    sheetname = str(num)+' Normalized ('+mode+')'
                    xlfile.writersim(sp,rtime[sumkey+mode],normkey,sheetname,mode) # write normalized data
        
        for key,val in sorted(sp.items()): # write isotope patterns
            if sp[key]['affin'] in mskeys:
                xlfile.writemultispectrum(sp[key]['spectrum'][0],sp[key]['spectrum'][1],'m/z','intensity','Isotope Patterns',key)
        
        if rd is None:
            for key,val in sorted(chroms.items()): # write chromatograms
                xlfile.writemultispectrum(chroms[key]['x'],chroms[key]['y'],chroms[key]['xunit'],chroms[key]['yunit'],'Function Chromatograms',key)
        
        uvstuff = False
        for key in sp: # check for UV-Vis spectra
            if sp[key]['affin'] is 'UV':
                uvstuff = True
                break
        if uvstuff is True:
            for ind,val in enumerate(tic['rawUV']): # normalize the UV intensities
                tic['rawUV'][ind] = val/1000000.
            xlfile.writersim(sp,rtime['rawUV'],'raw','UV-Vis','UV',tic['rawUV']) # write UV-Vis data to sheet
        
        if sumspec is not None:
            xlfile.writespectrum(sumspec[0],sumspec[1],'Summed Spectrum','m/z','counts')
            
        sys.stdout.write(' DONE\n')        
           
    def prepformula(dct):
        """looks for formulas in a dictionary and prepares them for pullspeciesdata"""
        for species in dct:
            if dct[species]['formula'] is not None:
                try:
                    dct[species]['mol'].res = res # sets resolution in Molecule object
                except NameError:
                    res = int(mzml.autoresolution())
                    dct[species]['mol'].res = res
                dct[species]['mol'].sigma = dct[species]['mol'].sigmafwhm()[1] # recalculates sigma with new resolution
                dct[species]['bounds'] = dct[species]['mol'].bounds(0.95) # caclulates bounds
            dct[species]['spectrum'] = Spectrum(3,dct[species]['bounds'][0],dct[species]['bounds'][1]) # generates a Spectrum object with those bounds
        return dct
    
    # ----------------------------------------------------------
    # -------------------PROGRAM BEGINS-------------------------
    # ----------------------------------------------------------
    global tome_v02,_ScriptTime,_mzML,_Spectrum,_Molecule,_XLSX
    from tome_v02 import bindata
    from _classes._ScriptTime import ScriptTime
    from _classes._mzML import mzML
    from _classes._Spectrum import Spectrum
    from _classes._Molecule import Molecule
    from _classes._XLSX import XLSX
    
    stime = ScriptTime()
    stime.printstart()
    
    n = checkinteger(n,'number of scans to sum') # checks integer input and converts to list
    
    sys.stdout.write('Loading processing parameters from excel file')
    sys.stdout.flush()
    xlfile = XLSX(xlsx)
    sp = xlfile.pullrsimparams()
    
    mskeys = ['+','-']
    for key in sp: # append list places for chrom, summed chrom, and normalized chrom
        sp[key]['raw'] = []
        if sp[key]['formula'] is not None: # if formula is specified
            sp[key]['mol'] = Molecule(sp[key]['formula']) # create Molecule object
            #sp[key]['bounds'] = sp[key]['mol'].bounds(0.99) # generate bounds from molecule object with this confidence interval
        if sp[key]['affin'] in mskeys:
            #sp[key]['spectrum'] = Spectrum(3,startmz=sp[key]['bounds'][0],endmz=sp[key]['bounds'][1])
            for num in n:
                sp[key]['%s' %(str(num)+'sum')] = []
                sp[key]['%s' %(str(num)+'norm')] = []
    sys.stdout.write(' DONE\n')
    
    
    rtime = {} # empty dictionaries for time and tic
    tic = {}
    rd = None
    for mode in mskeys: # look for existing positive and negative mode raw data
        try:
            rd = xlfile.wb.get_sheet_by_name('Raw Data ('+mode+')')
        except KeyError:
            continue
        modekey = 'raw'+mode
        if rd is not None: # if raw data is present, grab for processing
            sys.stdout.write('Existing (%s) mode raw data were found, grabbing those values.'%mode)
            sys.stdout.flush()
            rtime[modekey] = [] # generate empty lists required for data processing
            tic[modekey] = []
            for col,colval in enumerate(rd.columns):
                for row,rowval in enumerate(colval):
                    if row == 0: #skip first row
                        continue
                    elif colval[0].value == 'Time': # if column is Time, append to that list
                        rtime[modekey].append(rd.cell(row = (row+1), column = (col+1)).value)
                    elif colval[0].value == 'tic': # if column is tic, append to that list
                        tic[modekey].append(rd.cell(row = (row+1), column = (col+1)).value)
                    else: # all other columns
                        sp['%s' %str(colval[0].value)]['raw'].append(rd.cell(row = (row+1), column = (col+1)).value)
                if colval[0].value not in ['Time','tic']: # define affinity of species
                    sp['%s' %str(colval[0].value)]['affin'] = mode
            sys.stdout.write(' DONE\n')
            
    
    newpeaks = False
    if rd is not None:
        newsp = {}
        sumspec = None
        for key in sp: # checks whether there is a MS species that does not have raw data
            if len(sp[key]['raw']) is 0 and sp[key]['affin'] is not 'UV': #!!!!!! check whether the not UV if is needed
                newsp[key] = sp[key] # create references in the namespace
        if len(newsp) is not 0:
            newpeaks = True
            sys.stdout.write('Some peaks are not in the raw data, extracting these from raw file.\n')
            ips = xlfile.pullmultispectrum('Isotope Patterns') # pull predefined isotope patterns and add them to species
            for species in ips: # set spectrum list
                sp[species]['spectrum'] = [ips[species]['x'],ips[species]['y']]
            mzml = mzML(filename) # load mzML class
            newsp = prepformula(newsp) # prep formula species for summing
            for species in newsp:
                if newsp[species].has_key('spectrum') is False:
                    newsp[species]['spectrum'] = Spectrum(3,newsp[species]['bounds'][0],newsp[species]['bounds'][1])
            newsp = mzml.pullspeciesdata(newsp) # pull data
        else:
            sys.stdout.write('No new peaks were specified. Proceeding directly to summing and normalization.\n')
        
    if rd is None: # if no raw data is present, process mzML file
        mzml = mzML(filename) # load mzML class
        sp = prepformula(sp)
        sp,sumspec = mzml.pullspeciesdata(sp,True) # pull relevant data from mzML
        #sp,tic,rtime,sumspec = mzml.pullspeciesdata(sp,True) # pull relevant data from mzML
        chroms = mzml.pullchromdata() # pull chromatograms from mzML
        rtime = {}
        tic = {}
        for key in sp: # compare predicted isotope patterns to the real spectrum and save standard error of the regression
            func = sp[key]['function']
            if mzml.functions[func]['type'] == 'MS':
                mode = 'raw'+mzml.functions[func]['mode']
            if mzml.functions[func]['type'] == 'UV':
                mode = 'rawUV'
            if mode not in rtime:
                rtime[mode] = mzml.functions[func]['timepoints']
                tic[mode] = mzml.functions[func]['tic']
            if sp[key]['formula'] is not None:
                sp[key]['match'] = sp[key]['mol'].compare(sp[key]['spectrum'])
    
    if max(n) > 1: # run combine functions if n > 1
        for num in n: # for each n to sum
            sys.stdout.write('\r%s Summing species traces.' %str(n)[1:-1])
            sumkey = str(num)+'sum'
            for ind,key in enumerate(sp): # bin each species
                if sp[key]['affin'] in mskeys: # if species is MS related
                    sp[key][sumkey] = bindata(num,1,sp[key]['raw'])
            for mode in mskeys: 
                sumkey = str(num)+'sum'+mode
                modekey = 'raw'+mode
                if modekey in rtime.keys(): # if there is data for that mode
                    rtime[sumkey] = bindata(num,num,rtime[modekey])
                    tic[sumkey] = bindata(num,1,tic[modekey])
        sys.stdout.write(' DONE\n')
        sys.stdout.flush()
    
    for num in n: # normalize each peak's chromatogram
        sys.stdout.write('\r%s Normalizing species traces.' %str(n)[1:-1])
        sys.stdout.flush()
        sumkey = str(num)+'sum'
        normkey = str(num)+'norm'
        for mode in mskeys:
            modekey = 'raw'+mode
            if modekey in rtime.keys(): # if there is data for that mode
                for key in sp: # for each species
                    func = sp[key]['function']
                    if sp[key]['affin'] in mskeys: # if species has affinity
                        sp[key][normkey] = []
                        for ind,val in enumerate(sp[key][sumkey]):
                            #sp[key][normkey].append(val/(mzml.function[func]['tic'][ind]+0.01)) #+0.01 to avoid div/0 errors
                            sp[key][normkey].append(val/(tic[sumkey+sp[key]['affin']][ind]+0.01)) #+0.01 to avoid div/0 errors
    sys.stdout.write(' DONE\n')
    
    
    
    #import pickle #pickle objects (for troubleshooting)
    #pickle.dump(rtime,open("rtime.p","wb"))
    #pickle.dump(tic,open("tic.p","wb"))
    #pickle.dump(chroms,open("chroms.p","wb"))
    #pickle.dump(sp,open("sp.p","wb"))
    
    output() # write data to excel file
    xlfile.updatersimparams(sp) # update summing parameters
    
    sys.stdout.write('\rSaving "%s" (this may take some time)' %xlfile.bookname)
    sys.stdout.flush()
    xlfile.save()
    sys.stdout.write(' DONE\n') 
    
    sys.stdout.write('Plotting traces')
    plots() # plots for quick review
    sys.stdout.write(' DONE\n')
    
    stime.printelapsed()           

import sys
if len(sys.argv) > 1: # if script was initiated from the command line, pull parameters from there
    try:
        pyrsir(sys.argv[1],sys.argv[2],sys.arg[3])
    except IndexError:
        raise IOError('The pyrsim function requires three inputs:\n- The raw filename\n- The excel parameters file\n- The number of scans to sum')

if __name__ == '__main__':    
    pyrsir(filename,xlsx,n)
    sys.stdout.write('fin.')
    sys.stdout.flush()
    import gc
    gc.collect()