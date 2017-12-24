'''termcolor'''

class Progress(object):
    def __init__(self,**kwargs):
        """
        A progress output object for use when performing a large number of 
        repititions of the same process and informing the user of progress 
        is desired. 
        
        
        **keyword arguments**
        
        endmsg: *string*
            The message that is printed when the ``fin()`` method is called. 
            Default: DONE
        
        first: *integer*
            The first iteration of the process. Default: 1
        
        fraction: *bool*
            Whether the fractional progress should be written in the printed 
            string. e.g. 1/10. Default: True
        
        hashes: *bool*
            Whether a hash-type progress bar should be written in the printed 
            string. e.g. |#####       | Default: False
        
        hashnum: *integer*
            If *hashes* is True, this is the width of the hash-type progress bar. 
            Default: 20
        
        last: *integer*
            The last iteration of the process. Default: 10
        
        percent: *bool*
            Whether the percent progress should be written in the printed string. 
            e.g. 10.0% Default: True
        
        rng: *bool*
            Whether the iteration range should be written in the printed string. 
            e.g. (1-10) Default: False
        
        string: *string*
            The string prefix that is written (this is usually information about 
            the process. Default: Processing iteration
        
        writeevery: *integer*
            Write an output every n calls. This can be used if an iteration is 
            rapid and printing every output is not particularly useful. 
            Default: 1
        
        
        **Examples**
        
        ::
        
            >>> prog = Progress()
            >>> prog.write(1)
            Processing iteration #1/10 0.0%
            
            >>> prog.write(7)
            Processing iteration #7/10 66.7%
            
            >>> for i in range(1,11):
                    prog.write(i)
                prog.fin()
            Processing iteration #10/10 100.0% DONE
                
        
        """
        self.kw = { # default keyword arguments
        'first':1, # the initial iteration
        'last': 10, # the final iteration
        'string': 'Processing iteration', # the string prefix that is returned
        'fraction': True, # whether to output the fraction of things completed
        'rng': False, # whether to output the range that the iterations span
        'percent': True, # whether to show percent completion
        'hashes': False, # whether to have a hash bar for progress
        'hashnum': 20, # width of the hash progress bar
        'endmsg': 'DONE', # end message
        'writeevery': 1, # write output every n calls
        }
        if set(kwargs.keys()) - set(self.kw.keys()): # check for invalid keyword arguments
            string = ''
            for i in set(kwargs.keys()) - set(self.kw.keys()):
                string += ' %s' % str(i)
            raise KeyError('Unsupported keyword argument(s): %s' %string)
        self.kw.update(kwargs) # update defaules with provided keyword arguments
        import sys
        self.wr = sys.stdout.write
        self.fl = sys.stdout.flush
        self.count = 0
        self.strlen = 0
        if self.kw['hashes'] is True:
            self.spinner = ['|','/','-','\\']
    
    def __str__(self):
        """returns the progress string at the current iteration"""
        return self.write(self.current,True)
    
    def __repr__(self):
        return "%s(%s %d-%d)" %(self.__class__.__name__,self.kw['string'],self.kw['first'],self.kw['last'])
    
    def __getitem__(self,x):
        """
        Prints and returns the progress string at iteration x
        This accomplishes the same thing as write()
        """
        return self.write(x)
    
    def write(self,current,suppress=False):
        """
        writes the progress of the iteration
        suppress 
        """
        self.count += 1 # keep count
        if self.kw['writeevery'] != 1:
            if self.count != self.kw['last'] and self.count % self.kw['writeevery'] != 0: # if the counter does not match the write, bail out
                return None
        self.current = current # saves the current state
        string = '%s' %(self.kw['string']) # begin the string
        if self.kw['fraction'] is True:
            string += ' #%d/%d' %(current-self.kw['first']+1,self.kw['last']-self.kw['first']+1)
        if self.kw['rng'] is True:
            string += ' (%d-%d)' %(self.kw['first'],self.kw['last'])
        if self.kw['percent'] is True:
            try:
                self.perc = (float(current)-self.kw['first'])/float(self.kw['last']-self.kw['first'])*100.
                string += ' %.1f%%' %(self.perc)
            except ZeroDivisionError:
                string += ' err%%'
        if self.kw['hashes'] is True:
            string += self.hashes()
        if len(string) < self.strlen: # create space filler if output has somehow shrunk below the length of the previous output
            string += ' '* (self.strlen-len(string))
        self.strlen = len(string)
        if suppress is True: # does not write string to terminal, instead returns the generated progress string
            return string
        try:
            self.wr('\r%s' %string)
        except ValueError: # a catch for I/O errors that sometimes pop up for large iteration processes
            pass
        return string
    
    def hashes(self):
        """generates the hash-type progress bar if called for"""
        if self.kw['percent'] is False: # if the percent has not been calculated
            self.perc = (float(self.current)-self.kw['first'])/float(self.kw['last']-self.kw['first'])*100.
        out = ' |'
        num = int(self.perc/100.*self.kw['hashnum'])
        out += '#'*num # add completed
        if num < self.kw['hashnum']: # add spinner
            out += self.spinner[self.count%4]
        out += ' '*(self.kw['hashnum']-len(out)+2) # add still to go
        out += '|'
        return out
    
    def updatestring(self,string):
        """updates the output string to the specified string"""
        self.kw['string'] = string
    
    def fin(self,msg=None):
        """
        writes the completion message of the object and starts a new line
        if msg is specified, that message will be written instead of the object's 
        completion message
        """
        if msg is None:
            self.wr(' %s\n' %self.kw['endmsg'])
        else: # if a custom exit message was supplied
            self.wr(' %s\n' %msg)
        self.fl() # flush output
    


if __name__ == '__main__': # for testing and troubleshooting 
    first = 1
    finish = 255
    prog = Progress(
    first = first,
    last = finish,
    #hashes = True,
    )
    import time
    for i in range(first,finish+1):
        prog.write(i)
        time.sleep(0.01)
    prog.fin()