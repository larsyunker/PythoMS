"""
ScriptTime class
records timepoints in a python script

new:
    functional
    added periter to print the average time per supplied number of iterations
    created formattime to handle times less than 1 ms
    removed secondstostr and replaced all calls with formattime
    added function profiling
    added toggle for profiling
    changed from time.time() to time.clock() which seems to give much higher resolution
    ---1.0---
    removed timepoint function (now redundant with profiling capability)
    updated print profile data function to be more detailed and easier to read
    ---1.1---
    ---1.2

to add:
    use time.time() in unix and time.clock() in windows
"""
import sys
import time
import numpy as np
import datetime


class ScriptTime(object):
    def __init__(self, profile=False):
        """
        Class for storing timepoints in a python script
        profile toggles whether the profiling functionality of the class will operate
        """
        self._start_seconds = time.time()  # start time (seconds since epoch)
        self._end_seconds = None  # end time (seconds since epoch)
        self._start_clock = time.localtime()  # start clock time
        self._end_clock = None  # end clock time
        self.profile = profile  # toggle for profiling functions
        self.profiles = {}

    def __str__(self):
        """The string that is returned when printed"""
        return f'{self.__class__.__name__} initiated at {self.start_time}'

    def __repr__(self):
        """The representation that is returned"""
        return f'{self.__class__.__name__}({self.start_time})'

    @property
    def start_time(self):
        return time.strftime("%I:%M:%S %p", self._start_clock)

    @property
    def end_time(self):
        if self._end_clock is None:
            return None
        return time.strftime("%I:%M:%S %p", self._end_clock)

    @property
    def elapsed_time(self):
        if self._end_seconds is None:
            return time.time() - self._start_seconds
        return self._end_seconds - self._start_seconds

    def clearprofiles(self):
        """clears the profile data"""
        self.profiles = {}

    def formattime(self, t):
        """
        Formats a time value in seconds to the appropriate string. This is now included for legacy support.
        roughly based on formattime 
        """
        return str(datetime.timedelta(t))

    def periter(self, num):
        """
        calculated elapsed time per unit iteration (designed for timing scripts)
        """
        sys.stdout.write('Average time per iteration: %s\n' % (self.formattime(self.elapsed_time / float(num))))

    def printelapsed(self):
        """prints the elapsed time of the object"""
        if self._end_seconds is None:
            self.triggerend()
        sys.stdout.write(f'Elapsed time: {datetime.timedelta(seconds=self.elapsed_time)}\n')

    def printend(self):
        """prints the end time and the elapsed time of the object"""
        if self._end_seconds is None:
            self.triggerend()
        sys.stdout.write(f'End time: {self.end_time} (elapsed: {datetime.timedelta(seconds=self.elapsed_time)})\n')

    def printprofiles(self):
        """prints the data for the profiled functions"""
        sys.stdout.write('\nFunction profile data:\n')
        sys.stdout.write(
            '%15s  %6s  %13s  %13s  %13s  %13s\n' % ('function', 'called', 'avg', 'standard_deviation', 'max', 'min'))
        for fname, data in self.profiles.items():
            avg = sum(data[1]) / len(data[1])
            sys.stdout.write('%15s  %6d  %13s  %13s  %13s  %13s\n' % (
                fname,
                data[0],
                self.formattime(avg),
                self.formattime(np.sqrt(sum((i - avg) ** 2 for i in data[1]) / (len(data[1]) - 1))),
                self.formattime(max(data[1])),
                self.formattime(min(data[1]))
            ))

    def printstart(self):
        """prints the start (trigger) time of the object"""
        sys.stdout.write('Start time: %s\n' % (time.strftime('%I:%M:%S %p', self._start_clock)))

    def profilefn(self, fn):
        """generates a profiled version of the supplied function"""

        # from functools import wraps # unsure why these lines are present
        # @wraps(fn)
        def with_profiling(*args, **kwargs):
            """decorates function with profiling commands"""
            start_time = time.clock()  # time that the function was called
            ret = fn(*args, **kwargs)  # calls the function

            elapsed_time = time.clock() - start_time  # end time of the function
            if fn.__name__ not in self.profiles:  # generates a dictionary key based on the function name if not present
                self.profiles[fn.__name__] = [0, []]  # [number of times called, [list of durations]]
            self.profiles[fn.__name__][0] += 1
            self.profiles[fn.__name__][1].append(elapsed_time)
            return ret  # returns the calculated call of the function

        if self.profile is True:
            return with_profiling  # returns the decorated function
        else:
            return fn

    def triggerend(self):
        """triggers endpoint and calculates elapsed time since start"""
        self._end_seconds = time.time()
        self._end_clock = time.localtime()


if __name__ == '__main__':
    st = ScriptTime()
    time.sleep(2.)
    st.triggerend()
    st.printend()
