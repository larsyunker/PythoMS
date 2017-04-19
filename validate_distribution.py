"""
a script for validating the functionality of the commonly used scripts
"""

def test_pyrsir():
    sys.stdout.write('Testing PyRSIR...')
    from PyRSIR import pyrsir
    import shutil
    cwd = os.getcwd()
    old = cwd+'\\pyrsir_validation.xlsx'
    new = cwd+'\\pyrsir_validation (backup).xlsx'
    shutil.copy(old,new)
    pyrsir('MultiTest','pyrsir_validation',3,plot=False,verbose=False)
    shutil.copy(new,old)
    os.remove(new)
    sys.stdout.write(' PASS\n')    

def test_molecule():
    sys.stdout.write('Testing Molecule class...')
    from _classes._Molecule import Molecule
    mol1 = Molecule('L2PdAr+I')
    if mol1.sf != 'C61H51IP3Pd':
        raise ValueError('Bad string formula generation')
    if mol1.em != 1109.12832052557:
        raise ValueError('Bad exact mass calculation')
    if mol1.barip != [
    [1156.21,1157.21,1158.2001922415977,1159.2075274860124,1160.2033155288823,1161.2020504146174,1162.2003014864583,1163.200072566986,1164.2066774602326,1165.208663672193,1166.209393895099,1167.21,1168.21],
    [2.2945353033055147,1.5138412943167003,25.55112366668965,66.86927422750642,100.0,52.44501730463592,75.10171505680198,42.34934508332113,39.483066353566024,20.077321620243318,6.008127007653248,1.1883833361405927,0.1612963991856637]
    ]:
        raise ValueError('Bad bar isotope pattern generation')
    mol1 - 'PPh3' # test subtraction
    mol1 + 'PPh3' # test addition
    mol2 = Molecule('N(Et)2(CH2(13C)H2(2H))2')
    mol1 + mol2 # test class addition
    mol1.gaussianisotopepattern()
    sys.stdout.write(' PASS\n')
    
def test_mzml():
    sys.stdout.write('Testing mzML class...')
    from _classes._mzML import mzML
    mzml = mzML('MultiTest',verbose=False)
    if mzml.functions.keys() != [1,3,4]:
        raise ValueError('Did not pull the correct functions')
    @mzml._foreachchrom
    def testperchrom(chromatogram):
        attr = mzml.attributes(chromatogram)
        return attr['id']
    if testperchrom() != [u'TIC', u'SRM SIC Q1=200 Q3=100 function=2 offset=0']:
        raise ValueError('For each chromatogram or attributes function failed')
    @mzml._foreachscan
    def testperspec(spectrum):
        p = mzml.cvparam(spectrum)
        return p["MS:1000016"]
    if testperspec() != [0.0171000008, 0.135733336, 0.254333347, 0.372983336, 0.491699994, 0.0510833338, 0.169750005, 0.288383335, 0.407000005, 0.525833309, 0.0847499967, 0.20341666, 0.322033346, 0.440683335]:
        raise ValueError('For each scan or cvparam function failed')
    if sum(mzml.sum_scans()[1]) != 162806964:
        raise ValueError('sum_scans function failed')
    if sum((mzml[2])[1]) != 6742121:
        raise ValueError('scan indexing failed')
    if sum((mzml[0.01])[1]) != 56270834:
        raise ValueError('time indexing failed')
    sys.stdout.write(' PASS\n')
    
def test_xlsx():
    sys.stdout.write('Testing XLSX class...')
    from _classes._XLSX import XLSX
    xlfile = XLSX('xlsx_validation')
    spec,xunit,yunit = xlfile.pullspectrum('example MS spectrum')
    multispec = xlfile.pullmultispectrum('example multi-spectrum')
    rsimparams = xlfile.pullrsimparams()
    xlout = XLSX('xlsxtestout.xlsx',create=True)
    xlout.writespectrum(spec[0],spec[1],'test single spectrum out',xunit,yunit)
    for key,val in sorted(multispec.items()):
        xlout.writemultispectrum(multispec[key]['x'],multispec[key]['y'],multispec[key]['xunit'],multispec[key]['yunit'],'Function Chromatograms',key)
    xlout.save()
    os.remove(os.path.dirname(os.path.realpath(__file__))+'\\validation_files\\xlsxtestout.xlsx')
    sys.stdout.write(' PASS\n')
    
def test_spectrum():
    sys.stdout.write('Testing Spectrum class...')
    from _classes._Spectrum import Spectrum
    spec = Spectrum(3)
    spec.addvalue(479.1,1000)
    spec2 = Spectrum(3)
    spec2.addvalue(443.1,1000)
    spec += spec2
    spec3 = Spectrum(3,50,2500)
    spec3.addvalue(2150.954,1000)
    spec += spec3
    output = spec.trim(True)
    expected = [[50.0, 443.1, 479.1, 2150.954], [0, 1000, 1000, 1000]]
    if output != expected:
        raise ValueError('The output of the spectrum sequence did not match what was expected')
    sys.stdout.write(' PASS\n')
        
if __name__ == '__main__':
    import sys,os
    if os.path.dirname(os.path.realpath(__file__)) not in sys.path:
        sys.path.append(os.path.dirname(os.path.realpath(__file__)))
        sys.path.append(os.path.dirname(os.path.realpath(__file__))+'\\validation_files')
    test_molecule()
    test_mzml()
    test_spectrum()
    test_xlsx()
    test_pyrsir()