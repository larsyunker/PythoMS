"""
A script for plotting all spectra in an excel file as mass spectra and saving them
This allows for spectra without the same dimension to be compared

new:
    functional. creates simple figures with little customization
    ---1.0---
    now uses plotms
    ---1.1---
"""

import sys
from pythoms.xlsx import XLSX
from pythoms.tome import normalize, plot_mass_spectrum

# start m/z for the spectrum
startmz = 50

# end m/z for the spectrum
endmz = 700

# colour for the plotted spectrum (see _Colour for more details)
colour = 'k'

# the excel file name to find the x,y values for the spectra
filename = 'Dy'

# the type of the spectrum (either 'continuum' or 'centroid')
spectype = 'continuum'

#########################################################################################
# there is no need to modify the rest of the script
if __name__ == '__main__':
    xlfile = XLSX(filename)  # load excel file

    for sheet in xlfile.wb.sheetnames:  # for each sheet
        sys.stdout.write('Plotting spectrum "%s"' % sheet)
        sys.stdout.flush()
        spec = xlfile.pullspectrum(sheet, 8)[0]  # pull spectrum
        spec[1] = normalize(spec[1], 100.)  # normalize to 100
        outname = filename + ' - ' + sheet + ' (' + str(startmz) + '-' + str(endmz) + ')'
        plot_mass_spectrum(
            spec,
            outname=outname,
            spectype=spectype,
            mz=[startmz, endmz],
            padding=[0.1, 0.95, 0.13, 0.96],
            verbose=False,
            speccolour=colour
        )
        sys.stdout.write(' DONE\n')
    sys.stdout.write('fin.\n')
