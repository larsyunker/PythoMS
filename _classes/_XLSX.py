"""
Class for opening and handling excel files with commonly used data formats
v 1.0

new:
    loads workbook
    accepts incomplete extensions
    pulls spectra
    saves spectra
    pulls rsim parameters for use with summing program
    ---1.0---
    generalized and added rsim output
    generalized and created output for multiple spectra in a single sheet
    added pullmultispectrum function to read multiple spectra from sheet
    consolidated createwb into loadwb
    fixed workbook creation to not return a write-only workbood (this would break cell calls)
    ---1.1---
    added skiplines function to pullspectrum
    added pause for user input if the file is open (no longer requires rerunning entire script after closing)
    added function to convert row and column indicies into excel coordinates
    tweaked pullspectrum to warn users of non-numerical values
    ---1.2---
    user input in save() is now handled by a function
    ---1.3

to add:
    pull rsim raw data
    dynamic creation of workbook (if it's there, load it, if not create)
    pull from pyrsim output sheets
"""

class XLSX(object):
    def __init__(self,bookname,create=False):
        self.bookname = bookname
        self.loadop()
        self.wb,self.bookname = self.loadwb(self.bookname,create=create)

    def __str__(self):
        return 'Loaded excel file "%s"' %self.bookname
    
    def __repr__(self):
        return "%s('%s')" %(self.__class__.__name__,self.bookname)
    
    def checkduplicatesheet(self,sheet):
        """checks for duplicate sheets in the workbook and creates a unique name"""
        i = 1
        while sheet+' ('+`i`+')' in self.wb.get_sheet_names():
            i += 1
        return sheet+' ('+`i`+')'
    
    def correctextension(self,bookname):
        """attempts to correct the extension of the supplied filename"""
        oops = {'.xls':'x','.xl':'sx','.x':'lsx','.':'xlsx','.xlsx':''} # incomplete extensions
        for key in oops:
            if bookname.endswith(key) is True:
                bookname += oops[key]
                return bookname
        return bookname+'.xlsx' # absent extension
    
    def loadop(self):
        """loads openpyxl and checks for lxml"""
        try:
            self.op = __import__('openpyxl')
        except ImportError:
            raise ImportError('openpyxl does not appear to be installed.\nThe XLSX class requires this package to function, please install it.')
        try:
            import lxml
        except ImportError:
            raise ImportError('lxml does not appear to be installed.\nThe XLSX class requires this package to function, please install it.')
    
    def loadwb(self,bookname,create=False):
        """loads specified workbook into class"""
        try:
            wb = self.op.load_workbook(bookname) #try loading specified excel workbook
        except IOError:
            bookname = self.correctextension(bookname) # attempts to correct the extension of the provided workbook
            try:
                wb = self.op.load_workbook(bookname)
            except IOError:
                if create is True:
                    """
                    Due to write-only mode, creating and using that book breaks the cell calls
                    The workbook is therefore created, saved, and reloaded
                    The remove sheet call is to remove the default sheet
                    """
                    wb = self.op.Workbook(bookname,write_only=False) # create workbook
                    wb.save(bookname) # save it
                    wb = self.op.load_workbook(bookname) # load it
                    wb.remove_sheet(wb.worksheets[0]) # remove the old "Sheet"
                else:
                    raise IOError('\nThe excel file "%s" could not be found in the current working directory.' %(self.bookname))
        return wb,bookname
    
    def pullmultispectrum(self,sheetname):
        """reads multispectrum output back into dictionary format"""
        cs = self.wb.get_sheet_by_name(sheetname)
        out = {}
        loc = 1
        while loc < len(cs.rows[0]):
            out[cs.cell(row=1,column=loc).value] = {'xunit':cs.cell(row=1,column=loc+1).value,'yunit':cs.cell(row=1,column=loc+2).value,'x':[],'y':[]}
            ind = 2
            while cs.cell(row=ind,column=loc+1).value is not None:
                out[cs.cell(row=1,column=loc).value]['x'].append(cs.cell(row=ind,column=loc+1).value)
                out[cs.cell(row=1,column=loc).value]['y'].append(cs.cell(row=ind,column=loc+2).value)
                ind += 1
            loc += 4
        return out
        
    def pullspectrum(self,sheet='spectrum',skiplines=0):
        """
        extracts a spectrum from the specified sheet
        skiplines allows that number of lines to be ignored
        
        output: spectrum, xunit, yunit
        """
        def tofloat(value,row,col):
            """attempts to convert to float and raises exception if an error is encountered"""
            try:
                return float(value)
            except ValueError:
                raise ValueError('The value "%s" (cell %s) in "%s" could not be interpreted as a float.\nCheck the value in this cell or change the number of lines skipped' %(value,self.rowandcolumn(row,col),self.bookname))
        skiplines -= 1
        specsheet = self.wb.get_sheet_by_name(sheet)
        spectrum = [[],[]]
        for ind,row in enumerate(specsheet.rows): # for each row append the mz and int values to their respective lists
            if ind > skiplines: # skip specified number of lines
                if ind == skiplines+1: # header row
                    xunit = row[0].value
                    yunit = row[1].value
                    continue
                if row[0].value is not None and row[1].value is not None:
                    spectrum[0].append(tofloat(row[0].value,ind,0)) # append values
                    spectrum[1].append(tofloat(row[1].value,ind,1))
        return spectrum,xunit,yunit
    
    def pullrsimparams(self,sheet='parameters'):
        """
pulls parameters for reconstructed single ion monitoring processing
        
The expected structure for the excel sheet is:
row 1: headers only (not read by this script)

Column #1: name of species (this will be the column heading in the output sheets)
    (not needed if the molecular formula is specified)
Column #2: molecular formula of this species 
    (not required, but can bypass the need for start and end values)
Column #3: what spectrum to find this species in (+,-,UV)
Column #4: start value (m/z or wavelength)
Column #5: end value (m/z or wavelength)

The start value (col#4) is expected to be less than the end value (col#5)
        """
        def othernames(oldsheet):
            """tries to find another common name for the parameters sheet"""
            others = ['Parameters','params','Params']
            for name in others:
                if name in self.wb.get_sheet_names():
                    return name
            raise KeyError('There is no "%s" sheet in "%s".'%(oldsheet,self.bookname))
        
        if sheet not in self.wb.get_sheet_names(): # if the sheet can't be found, try other common names
            sheet = othernames(sheet)
        
        s = self.wb.get_sheet_by_name(sheet) #load sheet in specified excel file
        
        out = {} # output dictionary
        for ind,row in enumerate(s.rows):
            if ind == 0: # skip header row
                continue
            if row[0].value is None:
                name = str(row[1].value)
            else:
                name = str(row[0].value) # species name for dictionary key
            out[name] = {'bounds':[None,None],'affin':None,'formula':None,'match':None} # basic structure
            
            if row[1].value is not None: # set molecular formula if present
                out[name]['formula'] = row[1].value
            
            affin = {'(+)MS':['+','pos','positive','Pos','Positive'], # positive mode valid inputs
            '(-)MS':['-','neg','negative','Neg','Negative'], # negative mode valid inputs
            'UV':['UV','UV-Vis','UVVis','uv','uvvis','uv-vis']} # UV-Vis valid inputs
            try: # set affinity
                if row[2].value in affin['(+)MS']: # sets affinity to positive spectra
                    out[name]['affin'] = '+'
                elif row[2].value in affin['(-)MS']: # sets affinity to negative spectra
                    out[name]['affin'] = '-'
                elif row[2].value in affin['UV']: # sets affinity to negative spectra
                    out[name]['affin'] = 'UV'
                else: # sets affinity to positive (default)
                    out[name]['affin'] = '+'
            except IndexError: # if there is no affinity column
                out[name]['affin'] = '+'
            
            if row[3].value is not None: # set start value if present
                out[name]['bounds'][0] = float(row[3].value)
            
            if row[4].value is not None: # set end value if present
                out[name]['bounds'][1] = float(row[4].value)
            
        for key in out: #checks that start value is less than end value
            if out[key]['bounds'][0] is not None:
                if out[key]['affin'] == 'UV':
                    if out[key]['bounds'][1] is None: # ignore single value bounds
                        continue
                if out[key]['bounds'][1] is None:
                    raise ValueError('\nThere is no end value specified for row "%s".\nStart: %s\n%s'%(key,str(out[key]['bounds'][0]),self.pullrsimparams.__doc__))
                if out[key]['bounds'][0] > out[key]['bounds'][1]:
                    raise ValueError('\nThe end value is larger than the start value for row "%s".\nStart: %s\nEnd: %s\n%s'%(key,str(out[key]['bounds'][0]),str(out[key]['bounds'][1]),self.pullrsimparams.__doc__))
        return out
    
    def removesheets(self,delete):
        """
        removes sheets from the excel file
        delete is a set of strings to be removed
        """
        for sheet in self.wb.get_sheet_names(): # clears sheets that will contain new peak information
            if sheet in delete:
                dels = self.wb.get_sheet_by_name(sheet)
                self.wb.remove_sheet(dels)

    def rowandcolumn(self,row,col):
        """takes an index location of row and column and returns the cell location used by excel"""
        def modrem(val):
            return val//26,val%26
        import string
        alphabet = [x.upper() for x in list(string.ascii_lowercase)] # uppercase it
        col += 1 # offset column to be properly divisible by 26
        mod,rem = modrem(col)
        modl = []
        while mod > 26: # while modulo is greater than the length of the alphabet
            if rem == 0: # if it divided equally
                modl.insert(0,26)
                mod -= 1
            else:
                modl.insert(0,rem)
            mod,rem = modrem(mod)
        if mod == 1 and rem == 0: # exactly 26 == Z
            modl.insert(0,26)
        elif mod == 0: # less than 26
            modl.insert(0,rem)
        else: # other
            modl.insert(0,rem)
            modl.insert(0,mod)
        out = ''
        for i in modl: # build out string
            out += alphabet[i-1]
        return out+str(row+1)    
    
    def save(self,outname=None):
        """commits changes to the workbook"""
        def version_input(string):
            """checks the python version and uses the appropriate version of user input"""
            import sys
            if sys.version.startswith('2.7'):
                return raw_input('%s' %string)
            if sys.version.startswith('3.'):
                return input('%s' %string)
            else:
                raise EnvironmentError('The version_input method encountered an unsupported version of python.')
        
        if outname is None:
            outname = self.bookname
        
        try:
            self.wb.save(outname)
        except IOError:
            version_input('\nThe excel file could not be written. Please close "%s" and press any key to retry save.' %outname)
            try:
                self.wb.save(outname)
            except IOError:
                raise IOError('\nThe excel file "%s" could not be written.' %outname)
    
    def updatersimparams(self,sp,sheet='parameters'):
        """
        updates rsim parameters in the workbook
        sp is a dictionary of species
        """
        if sheet not in self.wb.get_sheet_names(): # if the sheet can't be found, try other common names
            sheet = self.pullrsimparameters.othernames(sheet)
        
        s = self.wb.get_sheet_by_name(sheet) #load sheet in specified excel file
        
        for ind,row in enumerate(s.rows):
            if ind == 0:
                continue
            if row[0].value is None: # if name is not set, use formula
                key = row[1].value
            else:
                key = str(row[0].value)
            if row[2].value is None: # if no affinity defined
                row[2].value = sp[key]['affin']
            if row[3].value is None: # if there are no bounds, update bounds
                row[3].value = sp[key]['bounds'][0]
            if row[4].value is None:
                row[4].value = sp[key]['bounds'][1]
            if sp[key]['match'] is not None:
                if s.cell(row=1,column=6).value is None:
                    s.cell(row=1,column=6).value = 'std err of reg'
                s.cell(row=ind+1,column=6).value = sp[key]['match']
        
    def writemultispectrum(self,xlist,ylist,xunit,yunit,sheetname,specname):
        """
        writes multiple spectra to a single sheet
        can be called multiple times and it will retain the placement
        
        output will be:
        specname|xunit   |yunit   | blank |...repeated
        blank   |xvalues |yvalues | blank |...repeated
        ...     |...     |...     | blank |...repeated
        
        """
        if self.__dict__.has_key('wms') is False: # check for dictionary in self
            self.wms = {sheetname:1}
        if self.wms.has_key(sheetname) is False: # check for preexisting key
            self.wms[sheetname] = 1
        
        if sheetname not in self.wb.get_sheet_names(): # if sheet does not exist
            cs = self.wb.create_sheet()
            cs.title = sheetname
        else:
            cs = self.wb.get_sheet_by_name(sheetname) # load existing sheet
        
        cs.cell(row = 1,column = self.wms[sheetname]).value = specname
        cs.cell(row = 1,column = self.wms[sheetname]+1).value = xunit
        cs.cell(row = 1,column = self.wms[sheetname]+2).value = yunit
        for ind,val in enumerate(xlist):
            cs.cell(row = ind+2,column = self.wms[sheetname]+1).value = val # write x value
            cs.cell(row = ind+2,column = self.wms[sheetname]+2).value = ylist[ind] # write y value
        self.wms[sheetname] += 4 # shift location over 4
    
    def writersim(self,sp,time,key,sheetname,mode,tic=None):
        """
        writes reconstructed single ion monitoring data to sheet
        
        sp is a dictionary of dictionaries for each species
            it is expected that the species' dictionary contains the specified key which is a 1D list
        time is a list of time values
        tic (if specified) is a list of total ion current values
        key is the name of the list within the species' dictionary
        sheetname is what the sheet will be named in the excel file
        mode is the current mode being output (usually either +,-,or UV)
        """
        if sheetname not in self.wb.get_sheet_names():
            cs = self.wb.create_sheet() #create new sheet
            cs.title = sheetname # rename sheet
            cs['A1'] = 'Time'
            offset = 0
            if tic is not None:
                offset += 1 # skip TIC column
                cs['B1'] = 'TIC'
            for ind,val in enumerate(time): #write time information
                cs.cell(row = (ind+2),column = 1).value = time[ind] #write time list
                if tic is not None:
                    cs.cell(row = (ind+2),column = 2).value = tic[ind] #write TIC list
            col = 1 + offset
            for species,dct in sorted(sp.items()):
                if sp[species]['affin'] is mode: # if the species' affinity is the mode
                    col+=1 # +1 from 0 to 1, +1 each to skip Time and TIC columns
                    cs.cell(row = 1,column = col).value = str(species) #write species names
                    for ind,val in enumerate(sp[species][key]):
                        cs.cell(row = (ind+2),column = col).value = sp[species][key][ind]
            
    def writespectrum(self,x,y,sheet='spectrum',xunit='m/z',yunit='counts'):
        """
        writes a provided spectrum to the specified sheet in the workbook
        x and y should be paired lists of values
        """
        if sheet in self.wb.get_sheet_names():
            sheet = self.checkduplicatesheet(sheet)
        ws = self.wb.create_sheet()
        ws.title = sheet
        ws['A1'] = xunit
        ws['B1'] = yunit
        for ind,val in enumerate(x):
            ws.cell(row = ind+2,column = 1).value = x[ind]
            ws.cell(row = ind+2,column = 2).value = y[ind]
        
if __name__ == '__main__':
    name = 'useless delete this'
    xl = XLSX(name,True)