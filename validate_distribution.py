"""
a script for validating the functionality of the commonly used scripts
"""


def test_pyrsir():
    sys.stdout.write('Testing PyRSIR...')
    from PyRSIR import pyrsir
    import shutil
    cwd = os.getcwd()
    old = cwd + '\\pyrsir_validation.xlsx'
    new = cwd + '\\pyrsir_validation (backup).xlsx'
    shutil.copy(old, new)
    pyrsir('MultiTest', 'pyrsir_validation', 3, plot=False, verbose=False)
    shutil.copy(new, old)
    os.remove(new)
    sys.stdout.write(' PASS\n')


def test_molecule():
    sys.stdout.write('Testing Molecule class...')
    from _classes._Molecule import Molecule
    mol1 = Molecule('L2PdAr+I',
                    # decpl=4,
                    # dropmethod='consolidate',
                    # threshold=0.01,
                    )
    if mol1.sf != 'C61H51IP3Pd':
        raise ValueError('Bad string formula generation')
    if mol1.em != 1109.1303776667514:
        raise ValueError('Bad exact mass calculation')
    if mol1.barip != [[1105.1304418,
                       1106.1338223494815,
                       1107.1290318632552,
                       1108.130515182168,
                       1109.1303776667514,
                       1110.1328896203245,
                       1111.1302404110806,
                       1112.132637641164,
                       1113.1319333869242,
                       1114.1341549662084,
                       1115.137151901555,
                       1116.1403684682987,
                       1117.1436663423106,
                       1118.1469889913205,
                       1119.1502750634809,
                       1120.1534694],
                      [2.2855073137216264,
                       1.5212910351820117,
                       25.459435969443284,
                       66.75651043139929,
                       100.0,
                       52.83757222099844,
                       75.16317471176205,
                       42.72403820863219,
                       39.685904925538146,
                       20.31615697845775,
                       6.1694175764333625,
                       1.2842197722332993,
                       0.20022034300192526,
                       0.024670825781130442,
                       0.0024698214870115874,
                       0.0001983611246228382]]:
        raise ValueError('Bad bar isotope pattern generation')
    mol1 - 'PPh3'  # test subtraction
    mol1 + 'PPh3'  # test addition
    mol2 = Molecule('N(Et)2(CH2(13C)H2(2H))2')
    mol1 + mol2  # test class addition
    mol1.gaussian_isotope_pattern(mol1.barip)
    sys.stdout.write(' PASS\n')


def test_mzml():
    sys.stdout.write('Testing mzML class...')
    from _classes._mzML import mzML
    mzml = mzML('MultiTest', verbose=False)
    if mzml.functions.keys() != [1, 3, 4]:
        raise ValueError('Did not pull the correct functions')

    @mzml.foreachchrom
    def testperchrom(chromatogram):
        attr = mzml.attributes(chromatogram)
        return attr['id']

    if testperchrom() != [u'TIC', u'SRM SIC Q1=200 Q3=100 function=2 offset=0']:
        raise ValueError('For each chromatogram or attributes function failed')

    @mzml.foreachscan
    def testperspec(spectrum):
        p = mzml.cvparam(spectrum)
        return p["MS:1000016"]

    if testperspec() != [0.0171000008, 0.135733336, 0.254333347, 0.372983336, 0.491699994, 0.0510833338, 0.169750005,
                         0.288383335, 0.407000005, 0.525833309, 0.0847499967, 0.20341666, 0.322033346, 0.440683335]:
        raise ValueError('For each scan or cvparam function failed')
    if sum(mzml.sum_scans()[1]) != 162804754:
        raise ValueError('sum_scans function failed')
    if sum((mzml[2])[1]) != 6742121:
        raise ValueError('scan indexing failed')
    if sum((mzml[0.01])[1]) != 56270834:
        raise ValueError('time indexing failed')
    sys.stdout.write(' PASS\n')


def test_xlsx():
    sys.stdout.write('Testing XLSX class...')
    from _classes._XLSX import XLSX
    xlfile = XLSX('xlsx_validation', verbose=False)
    spec, xunit, yunit = xlfile.pullspectrum('example MS spectrum')
    multispec = xlfile.pullmultispectrum('example multi-spectrum')
    rsimparams = xlfile.pullrsimparams()
    xlout = XLSX('xlsxtestout.xlsx', create=True, verbose=False)
    xlout.writespectrum(spec[0], spec[1], 'test single spectrum out', xunit, yunit)
    for key, val in sorted(multispec.items()):
        xlout.writemultispectrum(multispec[key]['x'], multispec[key]['y'], multispec[key]['xunit'],
                                 multispec[key]['yunit'], 'Function Chromatograms', key)
    xlout.save()
    os.remove(os.path.dirname(os.path.realpath(__file__)) + '\\validation_files\\xlsxtestout.xlsx')
    sys.stdout.write(' PASS\n')


def test_spectrum():
    sys.stdout.write('Testing Spectrum class...')
    from _classes._Spectrum import Spectrum
    spec = Spectrum(3)
    spec.addvalue(479.1, 1000)
    spec2 = Spectrum(3)
    spec2.addvalue(443.1, 1000)
    spec += spec2
    spec3 = Spectrum(3, start=50, end=2500)
    spec3.addvalue(2150.954, 1000)
    spec += spec3
    output = spec.trim(True)
    expected = [[50.0, 443.1, 479.1, 2150.954], [0, 1000, 1000, 1000]]
    if output != expected:
        raise ValueError('The output of the spectrum sequence did not match what was expected')
    sys.stdout.write(' PASS\n')


if __name__ == '__main__':
    import sys, os

    if os.path.dirname(os.path.realpath(__file__)) not in sys.path:
        sys.path.append(os.path.dirname(os.path.realpath(__file__)))
        sys.path.append(os.path.dirname(os.path.realpath(__file__)) + '\\validation_files')
    test_molecule()
    test_mzml()
    test_spectrum()
    test_xlsx()
    test_pyrsir()
