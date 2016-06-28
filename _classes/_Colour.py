class Colour(object):
    """
    script for taking a colour input and returning an RGB tuple usable by pyplot v002

    select a colour to use:
        A) supply a letter combination (ex. 'b' for blue or 'db' for dark blue)
        B) or supply a R,G,B value in (###,###,###) format [each R/G/B value should be an integer in 0-255]
        C) or supply a colour hash code (ex. '#881b1b' for (136,27,27))
    b: blue
    r: red
    g: olive green
    p: purple
    o: orange
    a: aqua
    k: black
    
    d_: dark colour (ex. 'db' = dark blue)
    l_: light colour (ex. 'lb' = light blue)
    """
    def __init__(self,c):
        self.c = c
        self.ty = type(self.c)
        if self.ty == str:
            if self.c[0] == '#':
                self.tup = self.hex_to_rgb(c)
            else:
                self.tup = self.colourdict(c)
        elif self.ty == tuple or self.ty == list:
            for val in self.c:
                if type(val) != int:
                    raise ValueError('One of the values in your colour tuple is not an integer ("%s").\nPlease correct the value.\n%s' %(val,self.__class__.__doc__))
                elif val > 255 or val < 0:
                    raise ValueError('One of the values in your colour tuple is outside of the valid RGB range of 0-255 ("%s").\n%s' %(val,self.__class__.__doc__))
            self.tup = self.c
        self.mpl = self.mplcolour(self.tup) # generate matplotlib colour tuple
    
    def __str__(self):
        return 'Colour {}'.format(self.c)
    
    def __repr__(self):
        return '{}({})'.format(self.__class__.__name__,self.c)
    
    def colourdict(self,string):
        """
        converts string key to predefined colour
        """
        coldict = {'b': (79,129,189),
        'r': (192,80,77),
        'g': (155,187,89),
        'p': (128,100,162),
        'a': (75,172,198),
        'o': (247,150,70),
        'lb': (149,179,215),
        'lr': (217,150,148),
        'lg': (195,214,155),
        'lp': (179,162,199),
        'lo': (147,205,221),
        'la': (250,192,144),
        'db': (54,96,146),
        'dr': (99,37,35),
        'dg': (79,98,40),
        'dp': (64,49,82),
        'do': (33,89,104),
        'da': (152,72,7),
        'k': (0,0,0)} 
        if coldict.has_key(string) is True:
            return coldict[string]
        else:
            raise ValueError('\nThe provided colour (%s) could not be interpreted.\n%s'%(str(string),self.__class__.__doc__))
    
    def hex_to_rgb(self,string):
        """
        converts hex colour to tuple
        borrowed from http://stackoverflow.com/questions/214359/converting-hex-color-to-rgb-and-vice-versa
        """
        string = string.lstrip('#')
        lv = len(string)
        return tuple(int(string[i:i + lv // 3], 16) for i in range(0, lv, lv // 3))
        #import struct
        #return struct.unpack('BBB',string.decode('hex'))
    
    def mplcolour(self,c):
        """
        converts 0-255 integers to 0-1 floats used by matplotlib
        """
        out = []
        for val in c: # normalize output (necessary if using matplotlib)
            out.append(round(val/255.,3))
        return tuple(out)
    
    def rgb_to_hex(self,tup):
        """
        converts RGB tuple to hex code
        """
        return '#%02x%02x%02x' % tup
        
        
if __name__ == '__main__':
    col = Colour('#881b1b')
    print col                
    
    
    