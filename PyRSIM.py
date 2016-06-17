"""
 Processes a .raw masslynx MS file or previously processed excel file with the parameters provided
 
 PyRSIM (Python Reconstructed Single Ion Monitoring) (previously SOAPy)
 version 027
 new:
     switched wb loading to function definition
     species output into excel is now sorted alphabetically
     restructured sp in dictionary to be a dictionary with keys for bounds, nsum, nnorm, etc
     now can do any number of n sums/norms in a single execution
     ---026.2---
     spectrum summing is incorporated into pullMSdata
     isotope patterns of each species are summed and saved in a separate sheet
     plot output changed to plot all summed traces (ignores norm traces)
     ---026.3---
     changed pullparams
     affinity of species (positive mode, negative mode, UV) is now set in the excel file in a new column
     all species to interpret are now given in the parameters sheet (no extra sheet required for UV-Vis data)
     now uses pullmzMLdata (rewritten pullMSdata to account for affinity)
     removed ext definition (assumption is that mzML files are supplied)
     new name for Raw Data excel sheets to account for both positive and negative mode
     works with TQD output
     can now process positive, negative mode as well as UV-Vis simultaneously
     excel output now groups species according to their affinity
     separated UV-Vis spec summing script from bulk script (will now have to be run separately)
     ---026.4---
     updated pullparams to be more efficient
     changed xml.dom.minidom import to be as xdm
     fixed plots function to only plot MS species (plots all + and - species on the same plot)
     detects when new peaks have been added and reprocesses mzML
     new scantype function determines what type of scan each spectrum is (MS+,MS-,UV,MSMS)
     fixed output functions to work when one MS mode is not present
     validated with TQD functions (extracts and outputs in the function chromatograms sheet)
     ---026.5---
     pullchromdata now outputs a dictionary of dictionaries (keys for x, y, xunit, and yunit)
     rewrote chromatogram output to work with dictionary
     chromatogram output is now sorted
     set up preliminary coding for dealing with MSMS spectra
     ---026.6---
     switched to use of mzML class and tome_v02
     updated and enabled command line initiation (it should work)
     added strtolist to handle list input from command line
     ---027.0---
     pwconvert now checks operating system
     switched to use of ScriptTime class
     ---027.1---
     switched import away from * (split up classes into separate files, etc.)
     switched to use of NoneSpectrum class for spectrum building
     ---027.2---
     renamed to PyRSIM
     updated input parameters
     pulls parameters from XLSX object
     creates Molecule objects if formula is specified
     if formula is specified, generates summing bounds from molecule
     updates parameters sheet with details if some cells were left blank
     removed Dandy call
     added error checking against isotope patterns (this might be broken)
     outputs standard error of the regression to the excel file if the pattern was compared
     updated to work with the latest versions of Spectrum, ScriptTime, and mzML
     moved rsim output to XLSX
     ---027.3---

to add/fix:
    update mzml to work with calibration
    finish off conversion to using XLSX
    create functionality for per-peak summing (create daughter dictionary?)
        if peaks overlap, combine
    
    fix integ() in mzML to work if there are no data points at the edge of the spectrum
    figure out how to find name of MSMS spectra
    functionality to extract and sum MSMS scans (rsim MSMS spectra)
    add output for images of isotope patterns and embed in excel (?possible?)
    figure out how to reinstate profiles (perhaps not useful anymore?)
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
#filename = 'LY-2016-05-18 09.raw'
filename = 'LY-2016-02-17 14.raw'
# Excel file to read from and output to (in *.xlsx format)
#xlsx = '2016-05-18 09'
xlsx = 'Book2 - Copy'
# set number of scans to sum (integer or list of integers)
n = [3,5]

# ----------------------------------------------------------
# -------------------FUNCTION DEFINITIONS-------------------
# ----------------------------------------------------------

def pyrsim(filename,xlsx,n):
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
        #pl.clf() # clears and closes old figure (if still open)
        #pl.close()
        nplots = len(n)+1
        
        # raw data
        pl.subplot(nplots,1,1) # top plot
        
        for mode in mskeys:
            modekey = 'raw'+mode
            if modekey in rtime.keys():
                pl.plot(rtime[modekey],TIC[modekey], linewidth = 0.75, label = 'TIC') #plot TIC
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
                    pl.plot(rtime[modekey],TIC[modekey], linewidth = 0.75, label = 'TIC') #plot TIC
                    for key in sp:
                        if sp[key]['affin'] is mode: #if a MS species
                            pl.plot(rtime[modekey],sp[key][sumkey], linewidth=0.75, label=key)
            pl.title('Summed Data (n=%i)' %(num))
            pl.ylabel('Intensity')
            pl.tick_params(axis='x',labelbottom='off')
            loc+=1
        pl.show()
  
    def output():
        """
        Function for writing the generated values to the excel file provided
        Each processed set of values is given its own sheet named accordingly
        Any function chromatograms are grouped into a single sheet called "Function Chromatograms"
        The output should not disturb any existing data, but will remove any charts present in the workbook
        
        requirements: openpyxl, lxml
        """
        if newpeaks is True:
            sys.stdout.write('Clearing duplicate xlsx sheets.')
            delete = []
            for key in sptemp: # generate strings to look for in excel file
                delete.append('Raw Data ('+sp[key]['affin']+')')
                for num in n:
                    delete.append(str(num)+' Sum ('+sp[key]['affin']+')')
                    delete.append(str(num)+' Normalized ('+sp[key]['affin']+')')
            delete.append('Isotope Patterns')
            #delete = set(delete)
            xlfile.removesheets(delete) # remove those sheets
            sys.stdout.write(' DONE.\n')
        
        sys.stdout.write('Writing to %s' %xlsx)
        sys.stdout.flush()
                
        for mode in mskeys: # write raw data to sheets
            modekey = 'raw'+mode
            if modekey in rtime.keys():
                sheetname = 'Raw Data ('+mode+')'
                xlfile.writersim(sp,rtime[modekey],'raw',sheetname,mode,tic=TIC[modekey])

        for num in n: # write summed and normalized data to sheets
            sumkey = str(num)+'sum'
            normkey = str(num)+'norm'
            for mode in mskeys:
                modekey = 'raw'+mode
                if modekey in rtime.keys():
                    if max(n) > 1: # if data were summed
                        sheetname = str(num)+' Sum ('+mode+')'
                        xlfile.writersim(sp,rtime[sumkey+mode],sumkey,sheetname,mode,tic=TIC[sumkey+mode]) # write summed data
                    sheetname = str(num)+' Normalized ('+mode+')'
                    xlfile.writersim(sp,rtime[sumkey+mode],normkey,sheetname,mode) # write normalized data
        
        if 'Isotope Patterns' not in xlfile.wb.get_sheet_names(): #check for preexisting isotope pattern sheet
            ip = xlfile.wb.create_sheet()
            ip.title = 'Isotope Patterns'
            loc = 1
            for key,val in sorted(sp.items()):
                if sp[key]['affin'] in mskeys:
                    ip.cell(row = 1,column = loc).value = key
                    ip.cell(row = 1,column = loc+1).value = 'intensity'
                    for ind,val in enumerate(sp[key]['spectrum'][0]):
                        ip.cell(row = ind+2,column = loc).value = val # write x value
                        ip.cell(row = ind+2,column = loc+1).value = sp[key]['spectrum'][1][ind] # write y value
                    loc += 3
        
        if 'Function Chromatograms' not in xlfile.wb.get_sheet_names(): #check for preexisting chromatogram sheet
            """
            outputs function chromatograms to a single sheet with the format
            name,space,space,repeat
            units,units,space,repeat
            vals,vals,space,repeat
            """
            fc = xlfile.wb.create_sheet()
            fc.title = 'Function Chromatograms'
            loc = 1
            for key,val in sorted(chroms.items()):
                fc.cell(row=1,column=loc).value = key # name of chromatogram
                fc.cell(row=2,column=loc).value = chroms[key]['xunit'] # x unit
                fc.cell(row=2,column=loc+1).value = chroms[key]['yunit'] # y unit
                for ind,val in enumerate(chroms[key]['x']): # values
                    fc.cell(row = ind+3,column=loc).value = chroms[key]['x'][ind]
                    fc.cell(row = ind+3,column=loc+1).value = chroms[key]['y'][ind]
                loc += 3
        
        uvstuff = False
        for key in sp:
            if sp[key]['affin'] is 'UV':
                uvstuff = True
                break
        if uvstuff is True:
            if 'UV-Vis' not in xlfile.wb.get_sheet_names():
                uvs = xlfile.wb.create_sheet()
                uvs.title = 'UV-Vis'
                uvs['A1'] = 'Time'
                uvs['B1'] = 'TIC'    
                for ind,val in enumerate(rtime['rawUV']): #write time information
                    uvs.cell(row = (ind+2),column = 1).value = rtime['rawUV'][ind] #write time list
                    uvs.cell(row = (ind+2),column = 2).value = TIC['rawUV'][ind]/1000000. #write TIC list
                ind = 2
                for key,val in sorted(sp.items()): # sorts and writes cells
                    if sp[key]['affin'] is 'UV':
                        ind +=1
                        uvs.cell(row = 1,column = ind).value = str(key)
                        for ind2,val2 in enumerate(sp[key]['raw']):
                            uvs.cell(row = (ind2+2),column = ind).value = sp[key]['raw'][ind2]
            
        sys.stdout.write(' DONE\n')        
           

    
    # ----------------------------------------------------------
    # -------------------PROGRAM BEGINS-------------------------
    # ----------------------------------------------------------
    
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
            sp[key]['bounds'] = sp[key]['mol'].bounds(conf=0.999) # generate bounds from molecule object
        if sp[key]['affin'] in mskeys:
            sp[key]['spectrum'] = Spectrum(3,startmz=sp[key]['bounds'][0],endmz=sp[key]['bounds'][1])
            for num in n:
                sp[key]['%s' %(str(num)+'sum')] = []
                sp[key]['%s' %(str(num)+'norm')] = []
    sys.stdout.write(' DONE\n')
    
    rtime = {} # empty dictionaries for time and TIC
    TIC = {}
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
            TIC[modekey] = []
            for col,colval in enumerate(rd.columns):
                for row,rowval in enumerate(colval):
                    if row == 0: #skip first row
                        continue
                    elif colval[0].value == 'Time': # if column is Time, append to that list
                        rtime[modekey].append(rd.cell(row = (row+1), column = (col+1)).value)
                    elif colval[0].value == 'TIC': # if column is TIC, append to that list
                        TIC[modekey].append(rd.cell(row = (row+1), column = (col+1)).value)
                    else: # all other columns
                        sp['%s' %str(colval[0].value)]['raw'].append(rd.cell(row = (row+1), column = (col+1)).value)
                if colval[0].value not in ['Time','TIC']: # define affinity of species
                    sp['%s' %str(colval[0].value)]['affin'] = mode
            sys.stdout.write(' DONE\n')
    
    newpeaks = False
    if rd is not None:
        sptemp = {}
        for key in sp: # checks whether there is a MS species that does not have raw data
            if len(sp[key]['raw']) is 0 and sp[key]['affin'] is not 'UV':
                sptemp[key] = sp[key]
        if len(sptemp.keys()) is not 0:
            newpeaks = True
            sys.stdout.write('Some peaks are not in the raw data, extracting these from raw file.\n')
            mzml = mzML(filename) # load mzML class
            sptemp,TIC,rtime = mzml.pullspeciesdata(sptemp) # pull data
            for key in sptemp:
                sp[key] = sptemp[key]
        
    
    
    if rd is None: # if no raw data is present, process mzML file
        mzml = mzML(filename) # load mzML class
        sp,TIC,rtime = mzml.pullspeciesdata(sp) # pull relevant data from mzML
        chroms = mzml.pullchromdata() # pull chromatograms from mzML
        for key in sp: # compare predicted isotope patterns to the real spectrum and save standard error of the regression
            if sp[key]['formula'] is not None:
                sp[key]['match'] = sp[key]['mol'].compare(sp[key]['spectrum'])
                sp[key]['mol'].plotgaus()
    
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
                    TIC[sumkey] = bindata(num,1,TIC[modekey])
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
                    if sp[key]['affin'] in mskeys: # if species has affinity
                        sp[key][normkey] = []
                        for ind,val in enumerate(sp[key][sumkey]):
                            sp[key][normkey].append(val/(TIC[sumkey+sp[key]['affin']][ind]+0.01)) #+0.01 to avoid div/0 errors
    sys.stdout.write(' DONE\n')
    
    
    
    #import pickle #pickle objects (for troubleshooting)
    #pickle.dump(rtime,open("rtime.p","wb"))
    #pickle.dump(TIC,open("TIC.p","wb"))
    #pickle.dump(chroms,open("chroms.p","wb"))
    #pickle.dump(sp,open("sp.p","wb"))
    
    output() # write data to excel file
    xlfile.updatersimparams(sp) # update summing parameters
    
    sys.stdout.write('\rSaving %s (this may take some time)' %xlsx)
    sys.stdout.flush()
    xlfile.save()
    sys.stdout.write(' DONE\n') 
    
    
    sys.stdout.write('Plotting traces')
    plots() # plots for quick review
    sys.stdout.write(' DONE\n')
    
    stime.printelapsed()           

import sys
#if len(sys.argv) > 1: # if script was initiated from the command line, pull parameters from there
#    import argparse
#    parser=argparse.ArgumentParser(
#        description='''A python script for processing .raw files''',
#        epilog="""Read the .py script for more information""")
#    parser.add_argument('{raw file}', type=str, help='name of the raw file to be processed')
#    parser.add_argument('{excel file}', type=str, help='name of the excel file with the processing parameters')
#    parser.add_argument('{numbers to bin}',type=int,default=1,help='Number of scans to bin (single integer or list')
#    args=parser.parse_args()
#    if len(sys.argv) != 1 and len(sys.argv) != 5:
#        sys.exit('Command line input requires exactly 3 parameters:\nRaw file name\nexcel file name\nnumber of scans to sum (integer or list)')
#    from tome_v02 import strtolist,bindata,PWconvert,loadwb,openpyxlcheck,pullparams
#    from _ScriptTime import ScriptTime
#    from _mzML import mzML
#    from _NoneSpectrum import NoneSpectrum
#    pyrsim(sys.argv[1],sys.argv[2],strtolist(sys.argv[3]))

if __name__ == '__main__':    
    import os
    sys.path.append(os.path.dirname(os.path.realpath(__file__))+'/_classes')
    from tome_v02 import bindata
    from _ScriptTime import ScriptTime
    from _mzML import mzML
    from _Spectrum import Spectrum
    from _Molecule import Molecule
    from _XLSX import XLSX
    pyrsim(filename,xlsx,n)
    sys.stdout.write('fin.')
    sys.stdout.flush()
    import gc
    gc.collect()