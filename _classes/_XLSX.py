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

to add:
    incorporate soapyoutput
    pull from pyrsim output sheets
"""

class XLSX(object):
    def __init__(self,bookname,create=False,interact=False):
        self.bookname = bookname
        self.loadop()
        if create is True:
            self.wb,self.bookname = self.createwb(self.bookname)
        if create is False:
            self.wb,self.bookname = self.loadwb(self.bookname,interact=interact)

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
    
    def createwb(self,bookname):
        """creates a workbook with the appropriate name"""
        if bookname.endswith('.xlsx') is False: # check file extension
            bookname =self.correctextension(bookname)
        return self.op.Workbook(bookname),bookname
    
    def correctextension(self,bookname):
        """attempts to correct the extension of the supplied filename"""
        oops = {'.xls':'x','.xl':'sx','.x':'lsx','.':'xlsx'} # incomplete extensions
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
    
    def loadwb(self,bookname,interact=False):
        """loads specified workbook into class"""
        try:
            wb = self.op.load_workbook(bookname) #try loading specified excel workbook
        except IOError:
            bookname = self.correctextension(bookname) # attempts to correct the extension of the provided workbook
            try:
                wb = self.op.load_workbook(bookname)
            except IOError:
                if interact is True:
                    if raw_input('The excel file "%s" could not be found in the current working directory, would you like to create it? (Y/N) '%bookname).lower() == 'y':
                        wb = self.createwb(bookname)[0]
                    else:
                        raise IOError('\nThe excel file "%s" could not be found in the current working directory.' %(self.bookname))
                else:
                    raise IOError('\nThe excel file "%s" could not be found in the current working directory.' %(self.bookname))
        return wb,bookname
    
    def pullspectrum(self,sheet='spectrum'):
        """extracts a spectrum from the specified sheet"""
        specsheet = self.wb.get_sheet_by_name(sheet)
        spectrum = [[],[]]
        for ind,row in enumerate(specsheet.rows): # for each row append the mz and int values to their respective lists
            if ind == 0: # header row
                xunit = row[0].value
                yunit = row[1].value
                continue
            spectrum[0].append(row[0].value) # append values
            spectrum[1].append(row[1].value)
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
    
    def save(self):
        """commits changes to the workbook"""
        try:
            self.wb.save(self.bookname)
        except IOError:
            raise IOError('\nThe excel file could not be written. Please close %s .' %self.bookname)
    
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
        ws.cell(row = 1,column = 1).value = xunit
        ws.cell(row = 1,column = 2).value = yunit
        for ind,val in enumerate(x):
            ws.cell(row = ind+2,column = 1).value = x[ind]
            ws.cell(row = ind+2,column = 2).value = y[ind]
        self.save()
        
if __name__ == '__main__':
    name = 'useless delete this'
    xl = XLSX(name)