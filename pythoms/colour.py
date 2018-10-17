"""
IGNORE:
takes an input string, tuple, or list and attempts to interpret it as a colour

new:
    several methods were rewritten to be less restrictive
    updated docstring
    added cmyk support (both to and from)
    added printdetails method for easy colour conversion
    made it so hex input does not need a '#'
    ---2.1---
    ---2.2
IGNORE
"""


class Colour(object):
    def __init__(self, c, rgb_scale=255, cmyk_scale=100):
        """
        Takes an input colour string and stores several colour properties of that colour

        **Parameters**

        c: *string* or *tuple*
            Input string (predefined colour or hex code) or R,G,B tuple.

        rgb_scale: *integer*
            The maximum value of the RGB scale for this object. Default is 255, implying
            RGB values between 0-255.

        cmyk_scale: *integer*
            The maximum value of the CMYK scale. Currently unused


        **Examples**

        ::

            >>> col = Colour((89,89,89))
            >>> col.print_details()
            Colour (89, 89, 89)
            RGB:        89, 89, 89
            RGB (0-1):  0.349, 0.349, 0.349
            CMYK:       0.0, 0.0, 0.0, 65.1
            hex:        #595959

            >>> col.mpl # matplotlib RGB format (floats between 0 and 1)
            (0.34901960784313724, 0.34901960784313724, 0.34901960784313724)

            >>> col2 = Colour('#A8E5B0')
            >>> col2.rgb
            (168, 229, 176)
            >>> col2.cmyk
            (26.637554585152838, 0.0, 23.144104803493455, 10.196078431372547)


        **Predefined Colours**

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
        self.c = c
        self.rgb_scale = float(rgb_scale)
        self.cmyk_scale = float(cmyk_scale)
        self.ty = type(self.c)
        if self.ty == str:  # interpret as string
            self.interpretstring(c)
        elif self.ty == tuple or self.ty == list:  # interpret as list or tuple
            if len(self.c) == 3:
                self.validatergb(self.c)
                self.rgb = self.c
                self.cmyk = self.rgb_to_cmyk(self.rgb)
            elif len(self.c) == 4:
                self.validatecmyk(self.c)
                self.cmyk = self.c
                self.rgb = self.cmyk_to_rgb(self.cmyk)
            self.hex = self.rgb_to_hex(self.rgb)
        self.mpl = self.mplcolour(self.rgb)  # generate matplotlib colour tuple
        self.rgb_to_others()  # calculate other values

    def __str__(self):
        return 'Colour {}'.format(self.c)

    def __repr__(self):
        return '{}({})'.format(self.__class__.__name__, self.c)

    def hex_to_rgb(self, string):
        """
        converts hex colour to tuple
        snagged from http://stackoverflow.com/questions/214359/converting-hex-color-to-rgb-and-vice-versa
        """
        string = string.lstrip('#')  # remove hash if present
        lv = len(string)
        if lv % 3 == 0:
            out = []
            for i in range(0, lv, lv // 3):
                out.append(int(string[i:i + lv // 3], 16))  # append integer in base 16
            return tuple(out)
            # return tuple(int(string[i:i + lv // 3], 16) for i in range(0, lv, lv // 3)) #list-comprehension version
        else:
            raise ValueError('The provided colour string "%s" could not be interpreted as a colour.\n%s' % (
            string, self.__class__.__doc__))

    def cmyk_to_rgb(self, tup):
        """converts cmyk tuple to rgb"""
        # set values and normalize
        c, m, y, k = tup
        c /= self.cmyk_scale
        m /= self.cmyk_scale
        y /= self.cmyk_scale
        k /= self.cmyk_scale

        # convert to RGB and scale
        r = int(round(self.rgb_scale * (1.0 - c) * (1.0 - k), 0))
        g = int(round(self.rgb_scale * (1.0 - m) * (1.0 - k), 0))
        b = int(round(self.rgb_scale * (1.0 - y) * (1.0 - k), 0))
        return r, g, b

    def interpretstring(self, string):
        """interprets a string as a colour"""
        colourdict = {
            'b': (79, 129, 189),
            'r': (192, 80, 77),
            'g': (155, 187, 89),
            'p': (128, 100, 162),
            'a': (75, 172, 198),
            'o': (247, 150, 70),
            'lb': (149, 179, 215),
            'lr': (217, 150, 148),
            'lg': (195, 214, 155),
            'lp': (179, 162, 199),
            'lo': (147, 205, 221),
            'la': (250, 192, 144),
            'db': (54, 96, 146),
            'dr': (99, 37, 35),
            'dg': (79, 98, 40),
            'dp': (64, 49, 82),
            'do': (33, 89, 104),
            'da': (152, 72, 7),
            'k': (0, 0, 0),
        }
        if string in colourdict:
            self.hex = self.rgb_to_hex(colourdict[string])
            self.rgb = colourdict[string]
        else:
            self.rgb = self.hex_to_rgb(string)
            self.hex = string
        self.cmyk = self.rgb_to_cmyk(self.rgb)

    def mplcolour(self, c):
        """converts 0-255 integers to 0-1 floats used by matplotlib"""
        out = []
        for val in c:  # normalize output to 255
            out.append(val / self.rgb_scale)
        return tuple(out)

    def printdetails(self):
        """prints the details of the object"""
        import sys
        sys.stdout.write('Colour %s\n' % str(self.c))
        sys.stdout.write('RGB:        %d, %d, %d\n' % self.rgb)
        sys.stdout.write('RGB (0-1):  %.3f, %.3f, %.3f\n' % self.mpl)
        sys.stdout.write('CMYK:       %.1f, %.1f, %.1f, %.1f\n' % self.cmyk)
        # sys.stdout.write('HLS:        %.3f, %.3f, %.3f\n' %self.hls)
        # sys.stdout.write('HSV:        %.3f, %.3f, %.3f\n' %self.hsv)
        # sys.stdout.write('YIQ:        %.3f, %.3f, %.3f\n' %self.yiq)
        sys.stdout.write('hex:        %s' % self.hex)

    def rgb_to_cmyk(self, tup):
        """converts rgb tuple to cmyk"""
        if sum(tup) == 0:  # black
            return 0, 0, 0, self.cmyk_scale

        # set values and normalize
        r, g, b = tup
        r /= self.rgb_scale
        g /= self.rgb_scale
        b /= self.rgb_scale

        # extract CMYK values
        k = 1 - max(r, g, b)
        c = (1 - r - k) / (1 - k)
        m = (1 - g - k) / (1 - k)
        y = (1 - b - k) / (1 - k)

        # scale and return
        return c * self.cmyk_scale, m * self.cmyk_scale, y * self.cmyk_scale, k * self.cmyk_scale

    def rgb_to_hex(self, tup):
        """converts RGB tuple to hex code"""
        return '#%02x%02x%02x' % tup

    def rgb_to_others(self):
        """
        converts RGB tuple to:
            HLS (Hue Lightness Saturation)
            HSV (Hue Saturation Value)
            YIQ"""
        import colorsys
        r, g, b = [i / self.rgb_scale for i in self.rgb]
        self.hls = colorsys.rgb_to_hls(r, g, b)
        self.hsv = colorsys.rgb_to_hsv(r, g, b)
        self.yiq = colorsys.rgb_to_yiq(r, g, b)

    def validatecmyk(self, tup):
        """checks that values in a CMYK tuple are valid"""
        for val in tup:
            if type(val) != int and type(val) != float:
                raise ValueError('One of the values in the CMYK colour tuple is not a number (%s).\n%s' % (
                str(val), self.__class__.__doc__))
            if val > self.cmyk_scale or val < 0:
                raise ValueError(
                    'One of the values in the CMYK colour tuple (%s) is outside of the valid range (0-%d)\n%s' % (
                    str(val), self.cmyk_scale, self.__class__.__doc__))

    def validatergb(self, tup):
        """checks that values in an RGB tuple are valid"""
        for val in tup:
            if type(val) != int:
                raise ValueError(
                    'One of the values in the RGB colour tuple is not an integer ("%s").\nPlease correct the value.\n%s' % (
                    str(val), self.__class__.__doc__))
            elif val > self.rgb_scale or val < 0:
                raise ValueError(
                    'One of the values in the RGB colour tuple (%s) is outside of the valid RGB range of 0-%d.\n%s' % (
                    str(val), self.rgb_scale, self.__class__.__doc__))


if __name__ == '__main__':
    col = Colour((89, 89, 89))
    col.printdetails()
