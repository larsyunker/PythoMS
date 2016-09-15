"""
Class for opening and handling excel files with commonly used data formats
v 1
CHANGELOG:

---1.3
"""

class XLSX(object):
    def __init__(self,bookname,**kwargs):
        self.ks = { # default keyword arguments
        'verbose': True, # toggle verbose
        'create': False, # create new workbook if supplied name is not found in directory
        }
        if set(kwargs.keys()) - set(self.ks.keys()): # check for invalid keyword arguments
            string = ''
            for i in set(kwargs.keys()) - set(self.ks.keys()):
                string += ` i`
            raise KeyError('Unsupported keyword argument(s): %s' %string)
        self.ks.update(kwargs) # update defaules with provided keyword arguments
        
        if self.ks['verbose'] is True:
            self.sys = __import__('sys')
        self.loadop() # check that lxml is present and load openpyxl
        self.wb,self.bookname = self.loadwb(bookname)

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
    
    def get_sheet(self,sheetname):
        """tries to retrieve the specified sheet name, otherwise returns None"""
        try:
            return self.wb.get_sheet_by_name(sheetname)
        except KeyError:
            return None
    
    def loadwb(self,bookname):
        """loads specified workbook into class"""
        if self.ks['verbose'] is True:
            self.sys.stdout.write('\rLoading workbook "%s" into memory' %bookname)
        try:
            wb = self.op.load_workbook(bookname) #try loading specified excel workbook
        except IOError:
            bookname = self.correctextension(bookname) # attempts to correct the extension of the provided workbook
            if self.ks['verbose'] is True:
                self.sys.stdout.write('\rLoading workbook "%s" into memory'%bookname)
            try:
                wb = self.op.load_workbook(bookname)
            except IOError:
                if self.ks['create'] is True:
                    """
                    Due to write-only mode, creating and using that book breaks the cell calls
                    The workbook is therefore created, saved, and reloaded
                    The remove sheet call is to remove the default sheet
                    """
                    if self.ks['verbose'] is True:
                        self.sys.stdout.write('Creating workbook "%s" and loading it into memory' % bookname)
                    wb = self.op.Workbook(bookname,write_only=False) # create workbook
                    wb.save(bookname) # save it
                    wb = self.op.load_workbook(bookname) # load it
                    wb.remove_sheet(wb.worksheets[0]) # remove the old "Sheet"
                else:
                    raise IOError('\nThe excel file "%s" could not be found in the current working directory.' %(self.bookname))
        if self.ks['verbose'] is True:
            self.sys.stdout.write(' DONE\n')
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
    
    def pullrsim(self,sheet):
        """
        pulls rsim data from the specified sheet
        """
        cs = self.wb.get_sheet_by_name(sheet)
        tic = []
        time = []
        data = {}
        for col,colval in enumerate(cs.columns):
            for row,rowval in enumerate(colval):
                if row == 0: #skip first row
                    continue
                elif colval[0].value == 'Time': # if column is Time, append to that list
                    time.append(cs.cell(row = (row+1), column = (col+1)).value)
                elif colval[0].value == 'TIC': # if column is tic, append to that list
                    tic.append(cs.cell(row = (row+1), column = (col+1)).value)
                else: # all other columns
                    if colval[row].value is not None:
                        if data.has_key(str(colval[0].value)) is False:
                            data['%s' %str(colval[0].value)] = {'raw':[]}
                        data['%s' %str(colval[0].value)]['raw'].append(cs.cell(row = (row+1), column = (col+1)).value)
        return data,time,tic  
                
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
        if self.ks['verbose'] is True:
            self.sys.stdout.write('Pulling spectrum from sheet "%s"' % sheet)
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
        if self.ks['verbose'] is True:
            self.sys.stdout.write(' DONE\n')
        return spectrum,xunit,yunit
    
    def pullrsimparams(self,sheet='parameters'):
        """
pulls parameters for reconstructed single ion monitoring processing
        
The expected structure for the excel sheet is:
row 1: headers only (not read by this script)
the headers tell the function what to find in that column

Valid column headers:
Name: name of the species
Formula: molecular formula of the species (if this is specified, it overrides start and end)
Function: what function to find the species in (optional if affinity is specified)
Affinity: what MS mode (+ or -) to find the species (or UV)
start: the start x value to find the species at
end: the end x value to find the species at (if end is specified, the function will integrate between those values; if end is not specified, the function will retrieve the closest value to the start value)

at least one of name, formula, start, or end must be specified for each species
        """
        def othernames(oldsheet):
            """tries to find another common name for the parameters sheet"""
            others = ['Parameters','params','Params']
            for name in others:
                if name in self.wb.get_sheet_names():
                    return name
            raise KeyError('There is no "%s" sheet in "%s".'%(oldsheet,self.bookname))
        
        def lowercase(string):
            """tries to return the lowercase of the supplied value, and avoids NoneType errors of .lower()"""
            if string is None:
                return None
            return string.lower()
        
        if sheet not in self.wb.get_sheet_names(): # if the sheet can't be found, try other common names
            sheet = othernames(sheet)
        
        s = self.wb.get_sheet_by_name(sheet) #load sheet in specified excel file
        
        names = ['name'] # valid name column headers
        formulas = ['formula','form','mf'] # valid molecular formula column headers
        affinities = ['affin','affinity'] # valid affinity column headers
        functions = ['fn','func','function'] # valid function column headers
        levels = ['level','lvl'] # valid level column headers
        smz = ['startmz','start mz','start m/z','startm/z','start'] #valid start m/z column headers
        emz = ['endmz','end mz','end m/z','endm/z','end'] # valid end m/z column headers
        
        self.rsimh = {}
        for ind,col in enumerate(s.columns):
            if lowercase(col[0].value) in names:
                self.rsimh['name'] = ind
            if lowercase(col[0].value) in formulas:
                self.rsimh['formula'] = ind
            if lowercase(col[0].value) in affinities:
                self.rsimh['affin'] = ind
            if lowercase(col[0].value) in functions:
                self.rsimh['function'] = ind
            if lowercase(col[0].value) in levels:
                self.rsimh['level'] = ind
            if lowercase(col[0].value) in smz:
                self.rsimh['bounds'] = [ind,None]
            if lowercase(col[0].value) in emz:
                self.rsimh['bounds'][1] = ind
        
        out = {} # output dictionary
        for ind,row in enumerate(s.rows):
            if ind == 0: # skip header row
                continue
            
            # figure out what to name the species
            name = None
            if self.rsimh.has_key('name') and row[self.rsimh['name']].value is not None: # if there is a name column and the value is not None
                name = row[self.rsimh['name']].value
            elif name is None and self.rsimh.has_key('formula') and row[self.rsimh['formula']].value is not None: # if there is instead a formula
                name = row[self.rsimh['formula']].value
            elif name is None and self.rsimh.has_key('bounds') and row[self.rsimh['bounds'][0]].value is not None: # if there are bounds
                if self.rsimh['bounds'][1] is None:
                    name = row[self.rsimh['bounds'][0]].value
                else:
                    name = '%.1f - %.1f' %(row[self.rsimh['bounds'][0]].value,row[self.rsimh['bounds'][0]].value)
            else:
                raise ValueError('A name could not be determined for row #%d\n%s' %(ind+1,self.pullrsimparams.__doc__))
            out[name] = {} # create name key
            
            for key in self.rsimh: # go through the headers for that row and pull and values provided
                if key == 'name':
                    continue
                if key == 'bounds':
                    if row[self.rsimh[key][0]].value is not None:
                        out[name]['bounds'] = [float(row[self.rsimh[key][0]].value),None]
                    if row[self.rsimh[key][1]].value is not None:
                        out[name]['bounds'][1] = float(row[self.rsimh[key][1]].value)
                    continue
                if row[self.rsimh[key]].value is not None:
                    out[name][key] = row[self.rsimh[key]].value
                
        pos = ['+','pos','positive','Pos','Positive'] # positive mode valid inputs
        neg = ['-','neg','negative','Neg','Negative'] # negative mode valid inputs
        uv = ['UV','UV-Vis','UVVis','uv','uvvis','uv-vis'] # UV-Vis valid inputs
        for name in out:
            if out[name].has_key('affin'): # change affinity to valid ones used by mzML
                if out[name]['affin'] in pos:
                    out[name]['affin'] = '+'
                elif out[name]['affin'] in neg:
                    out[name]['affin'] = '+'
                elif out[name]['affin'] in uv:
                    out[name]['affin'] = '+'
                else:
                    raise ValueError('The affinity "%s" for species "%s" is not valid\n%s' %(out[name]['affin'],name,self.pullrsimparams.__doc__))
            
            if out[name].has_key('function'): # if function is specified, convert to integer
                out[name]['function'] = int(out[name]['function'])
            
            if out[name].has_key('affin') is False and out[name].has_key('function') is False: # if no affinity or function, assume 1st function
                out[name]['function'] = 1
           
            if out[name].has_key('level'): # if level is specified, convert to integer
                out[name]['level'] = int(out[name]['level'])
            
            if out[name].has_key('bounds'): # if bounds are specified, check that end is not less than start
                if out[name]['bounds'][1] is None:
                    pass
                elif out[name]['bounds'][0] > out[name]['bounds'][1]:
                    raise ValueError('\nThe end value is larger than the start value for species "%s".\nStart: %s\nEnd: %s\n%s'%(name,str(out[name]['bounds'][0]),str(out[name]['bounds'][1]),self.pullrsimparams.__doc__))
            
            if not out[name].has_key('formula') and not out[name].has_key('bounds'):
                raise ValueError('The species "%s" does not have a formula or integration bounds' %(name))
            
            if out[name].has_key('formula') is False: # create formula key if not specified (required for pyrsir)
                out[name]['formula'] = None
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
    xl = XLSX(name,create=True)