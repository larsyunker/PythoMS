"""
This script generates then iterates through a mzML file, binning a specified number of spectra
CHANGELOG:
    ---4.2
"""


def sumspectra(filename, start=None, end=None, excel=None):
    """
    Sums spectra from raw file and outputs to excel file 
    
    input: filename, *kwargs
    filename:
        name of raw file
    sr:
        scan range to sum
        default 'all'
        specify with [start scan,end scan]
    """
    from PythoMS._classes._mzML import mzML
    from PythoMS._classes._XLSX import XLSX
    from PythoMS._classes._ScriptTime import ScriptTime

    st = ScriptTime()
    st.printstart()
    mzml = mzML(filename)  # create mzML object
    if start is None:
        start = mzml.functions[1]['sr'][0] + 1
    if end is None:
        end = mzml.functions[1]['sr'][1] + 1
    x, y = mzml.sum_scans(start=start, end=end)
    xlfile = XLSX(filename, create=True)
    xlfile.writespectrum(x, y, 'summed spectra (scans %d-%d)' % (start, end))
    xlfile.save()
    st.printend()


if __name__ == '__main__':
    sumspectra(
        'LY-2015-01-20 02',  # raw filename to use
        start=None,  # start scan number
        end=None,  # end scan number
    )
