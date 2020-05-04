import sys
import warnings

warnings.warn(
    'the pythoms.progress module as been deprecated, please switch to tqdm',
    DeprecationWarning,
    stacklevel=2,
)


class Progress(object):
    def __init__(self,
                 first: int = 1,  # the initial iteration
                 last: int = 10,  # the final iteration
                 string: str = 'Processing iteration',  # the string prefix that is returned
                 fraction: bool = True,  # whether to output the fraction of things completed
                 rng: bool = False,  # whether to output the range that the iterations span
                 percent: bool = True,  # whether to show percent completion
                 hash: bool = False,  # whether to have a hash bar for progress
                 hashnum: int = 20,  # width of the hash progress bar
                 endmsg: str = 'DONE',  # end message
                 writeevery: int = 1,  # write output every n calls
                 ):
        """
        A progress output object for use when performing a large number of
        repititions of the same process and informing the user of progress
        is desired.

        :param first: The first iteration of the process.
        :param last: The last iteration of the process.
        :param string: The string prefix that is written (this is usually information about
            the process.
        :param fraction: Whether the fractional progress should be written in the printed string. e.g. 1/10.
        :param rng: Whether the iteration range should be written in the printed string.
            e.g. (1-10)
        :param percent: Whether the percent progress should be written in the printed string.
            e.g. 10.0%
        :param hash: Whether a hash-type progress bar should be written in the printed
            string. e.g. |#####       | Default: False
        :param hashnum: If *hashes* is True, this is the width of the hash-type progress bar.
        :param endmsg: The message that is printed when the `fin()` method is called.
        :param writeevery: Write an output every n calls. This can be used if an iteration is
            rapid and printing every output is not particularly useful.
        
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
        self.first = first
        self.last = last
        self.string = string
        self.fraction = fraction
        self.rng = rng
        self.percent = percent
        self.hash = hash
        self.hashnum = hashnum
        self.endmsg = endmsg
        self.writeevery = writeevery
        self.wr = sys.stdout.write
        self.fl = sys.stdout.flush
        self.count = 0
        self.strlen = 0
        self.current = 0  # current state
        self.spinner = ['|', '/', '-', '\\']

    def __str__(self):
        """returns the progress string at the current iteration"""
        return self.write(self.current, True)

    def __repr__(self):
        return f'{self.__class__.__name__}({self.string} {self.first}-{self.last})'

    def __getitem__(self, x):
        """
        Prints and returns the progress string at iteration x
        This accomplishes the same thing as write()
        """
        return self.write(x)

    @property
    def perc(self):
        try:
            return round(
                (float(self.current) - self.first)
                / float(self.last - self.first)
                * 100.,
                1
            )
        except ZeroDivisionError:
            return 0.

    def write(self, current, suppress=False):
        """
        Writes the progress of the iteration

        :param current: current iteration
        :param suppress: suppress output
        :return: formatted output string
        """
        self.count += 1  # keep count
        if self.writeevery != 1:
            # if the counter does not match the write, bail out
            if self.count != self.last and self.count % self.writeevery != 0:
                return None
        self.current = current  # saves the current state
        string = f'{self.string}'  # begin the string
        if self.fraction is True:
            string += f' #{current - self.first + 1}/{self.last - self.first + 1}'
        if self.rng is True:
            string += f' ({self.first}-{self.last})'
        if self.percent is True:
            string += f' {self.perc}%'
        if self.hash is True:
            string += self.hashes()
        # create space filler if output has somehow shrunk below the length of the previous output
        if len(string) < self.strlen:
            string += ' ' * (self.strlen - len(string))
        self.strlen = len(string)
        # does not write string to terminal, instead returns the generated progress string
        if suppress is True:
            return string
        try:
            self.wr(f'\r{string}')
        # a catch for I/O errors that sometimes pop up for large iteration processes
        except ValueError:
            pass
        return string

    def hashes(self):
        """generates the hash-type progress bar if called for"""
        out = ' |'
        num = int(self.perc / 100. * self.hashnum)
        out += '#' * num  # add completed
        if num < self.hashnum:  # add spinner
            out += self.spinner[self.count % 4]
        out += ' ' * (self.hashnum - len(out) + 2)  # add still to go
        out += '|'
        return out

    def fin(self, msg=None):
        """
        writes the completion message of the object and starts a new line
        if msg is specified, that message will be written instead of the object's 
        completion message
        """
        if msg is None:
           msg = self.endmsg
        self.wr(f' {msg}\n')
        self.fl()  # flush output


if __name__ == '__main__':  # for testing and troubleshooting
    first = 1
    finish = 255
    prog = Progress(
        first=first,
        last=finish,
        # hashes = True,
    )
    import time

    for i in range(first, finish + 1):
        prog.write(i)
        time.sleep(0.01)
    prog.fin()
