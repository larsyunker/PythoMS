"""
ScriptTime class
records timepoints in a python script

CHANGELOG
---1.1---
- added plots method so the user can check whether the function call time increases over time
- modified profile printing to return 'n/a' if the number of calls is less than 3
---1.2

to add:
    use time.time() in unix and time.clock() in windows
"""
    

class ScriptTime(object):
    def __init__(self, profile=False):
        """
        Class for storing timepoints in a python script
        profile toggles whether the profiling functionality of the class will operate
        """
        self.time = __import__('time')
        self.sys = __import__('sys')
        self.time_zero = self.time.clock()
        self.start = self.time.localtime()
        self.profile = profile # toggle for profiling functions
        self.profiles = {}
        
    def __str__(self):
        """The string that is returned when printed"""
        return "The initiation time of this class was {}".format(self.time.strftime('%I:%M:%S %p',self.start))
    
    def __repr__(self):
        """The representation that is returned"""
        return "{}({})".format(self.__class__.__name__,self.time.strftime('%I:%M:%S %p',self.start))
        
    def clearprofiles(self):
        """clears the profile data"""
        self.profiles = {}
        
    def formattime(self,t):
        """
        formats a time value in seconds to the appropriate string
        roughly based on formattime 
        """
        def rediv(b,lst):
            """elongates a list based on the quotient and remainder of the value and supplied list"""
            return list(divmod(b[0],lst))+b[1:]
        
        def hour(lst):
            """formats to the hour"""
            return "%d:%02d:%02d.%03d" % (lst[0],lst[1],lst[2],round(lst[3]))
        
        def minute(lst):
            """formats to the minute"""
            return "%02d:%02d.%01d" %(lst[1],lst[2],round(lst[3]))

        def second(lst):
            """formats to the second"""
            return "%d.%03d s" %(lst[2],round(lst[3]))
        
        def lessthansecond(lst):
            """used to return a str with appropriate units if the time is less than a second"""
            keys = {3:'ms',6:'us',9:'ns',12:'ps'} # units for given orders of magnitude
            return "%.1f %s" %(round(lst[-1],1),keys[order])
        
        order = 3
        out = reduce(rediv,[[t*10**order,],1000,60,60])
        if out[0] == 0: # if less than an hour
            if out[1] == 0: # if less than a minute
                if out[2] == 0: # if less than a second
                    if out[3] // 1 == 0: # if less than 1 millisecond
                        while out[-1] // 1 == 0: # determine the appropriate order of magnitude
                            order += 3
                            out = reduce(rediv,[[t*10**order,],1000,60,60])
                            if order > 9: # if order is greater than handled (typically systems can only track to the microsecond)
                                return '<1 ns'
                        return lessthansecond(out)
                    else:
                        return lessthansecond(out)
                else:
                    return second(out)
            else:
                return minute(out)
        else:
            return hour(out)
    
    def periter(self,num):
        """
        calculated elapsed time per unit iteration (designed for timing scripts)
        """
        self.sys.stdout.write('Average time per iteration: %s\n' %(self.formattime(self.elap/float(num))))
    
    def plots(self):
        """plots the duration of each function call against the function call number"""
        import pylab as pl
        import scipy as sp
        for key in self.profiles:
            pl.plot(sp.arange(self.profiles[key][0])+1,self.profiles[key][1],label=key)
        pl.legend()
        pl.show()
    
    def printelapsed(self):
        """prints the elapsed time of the object"""
        if self.__dict__.has_key('end_time') is False:
            self.triggerend()
        self.sys.stdout.write('Elapsed time: %s\n' %(self.formattime(self.elap)))
    
    def printend(self):
        """prints the end time and the elapsed time of the object"""
        if self.__dict__.has_key('end_time') is False:
            self.triggerend()
        self.sys.stdout.write('End time: %s (elapsed: %s)\n' %(self.time.strftime('%I:%M:%S %p',self.end),self.formattime(self.elap)))
    
    def printprofiles(self):
        """prints the data for the profiled functions"""
        if self.__dict__.has_key('m') is False:
            self.m = __import__('math')
        self.sys.stdout.write('\nFunction profile data:\n')
        self.sys.stdout.write('%15s  %6s  %13s  %13s  %13s  %13s\n' %('function','called','avg','stdev','max','min'))
        for fname, data in self.profiles.items():
            avg = sum(data[1])/len(data[1])
            if len(data[1]) < 3: # if the number of calls is less than 3, standard deviation is not applicable
                stdev = 'n/a'
            else:
                stdev = self.formattime(self.m.sqrt(sum((i-avg)**2 for i in data[1])/(len(data[1])-1)))
            self.sys.stdout.write('%15s  %6d  %13s  %13s  %13s  %13s\n' %(fname,data[0], self.formattime(avg), stdev, self.formattime(max(data[1])), self.formattime(min(data[1]))))
            #self.sys.stdout.write('Function %s called %d times. ' % (fname, data[0]))
            #self.sys.stdout.write('Execution time max: %s, min: %s, average: %s, stdev: %s\n' % (self.formattime(max(data[1])), self.formattime(min(data[1])), self.formattime(avg),self.formattime(self.m.sqrt(sum((i-avg)**2 for i in data[1])/(len(data[1])-1))) ))
    
    def printstart(self):
        """prints the start (trigger) time of the object"""
        self.sys.stdout.write('Start time: %s\n' %(self.time.strftime('%I:%M:%S %p',self.start)))
    
    def profilefn(self,fn):
        """generates a profiled version of the supplied function"""
        #from functools import wraps # unsure why these lines are present
        #@wraps(fn)
        def with_profiling(*args, **kwargs):
            """decorates function with profiling commands"""
            start_time = self.time.clock() # time that the function was called
            ret = fn(*args, **kwargs) # calls the function
    
            elapsed_time = self.time.clock() - start_time # end time of the function
            if fn.__name__ not in self.profiles: # generates a dictionary key based on the function name if not present
                self.profiles[fn.__name__] = [0, []] # [number of times called, [list of durations]]
            self.profiles[fn.__name__][0] += 1
            self.profiles[fn.__name__][1].append(elapsed_time)
            return ret # returns the calculated call of the function
        if self.profile is True:
            return with_profiling # returns the decorated function
        else:
            return fn
    
    def triggerend(self):
        """triggers endpoint and calculates elapsed time since start"""
        self.time_end = self.time.clock()
        self.end = self.time.localtime()
        self.elap = self.time_end - self.time_zero
    
if __name__ == '__main__':
    st = ScriptTime(profile=True)
    t = 0.236592
    print st.formattime(t)
    #st.printstart()
    #
    #@st.profilefn
    #def first(x):
    #    return x+100
    #
    #@st.profilefn
    #def second(x):
    #    return x**2
    #i = 10
    #for i in range(1000):
    #    i = first(i)
    #    i = second(i)
    #st.printend()
    #st.printprofiles()
    
    
    