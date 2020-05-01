"""
a script for validating the functionality of the commonly used scripts
"""
import unittest
import os
import pathlib
import shutil
from random import random
from pythoms.mzml import mzML, branch_attributes
from pythoms.molecule import Molecule, IPMolecule, VALID_DROPMETHODS, VALID_IPMETHODS, element_intensity_list
from pythoms.spectrum import Spectrum
from pythoms.xlsx import XLSX
from pythoms.psims import CVParameterSet
from PyRSIR import pyrsir


validation_path = pathlib.Path(os.getcwd()) / 'validation_files'


class TestPyRSIR(unittest.TestCase):
    def setUp(self) -> None:
        # target xlsx and mzml pairings
        self.targets = [
            [
                'multitest_pyrsir_validation.xlsx',
                'MultiTest'
            ],
            [
                'LY-2015-09-15 06 pyrsir example.xlsx',
                'LY-2015-09-15 06',
            ]
        ]
        # create backups
        for xlfile, _ in self.targets:
            shutil.copy(
                validation_path / xlfile,
                validation_path / f'{xlfile}.bak'
            )

    def test(self):
        for xlfile, msfile in self.targets:
            pyrsir(
                validation_path / msfile,
                validation_path / xlfile,
                3,
                plot=False,
                verbose=False,
            )

    def tearDown(self) -> None:
        for xlfile, _ in self.targets:
            shutil.copy(
                validation_path / f'{xlfile}.bak',
                validation_path / xlfile,
            )
            os.remove(validation_path / f'{xlfile}.bak')


class TestMolecule(unittest.TestCase):
    def setUp(self):
        self.mol = Molecule('L2PdAr+I')
        self.ipmol = IPMolecule(
            'L2PdAr+I',
            ipmethod='multiplicative',
            dropmethod='threshold',
            threshold=0.01,
        )

    def test_molecule(self):
        self.assertEqual(
            self.mol.molecular_formula,
            'C61H51IP3Pd'
        )
        self.assertEqual(
            self.mol.composition,
            {'C': 61, 'H': 51, 'P': 3, 'Pd': 1, 'I': 1}
        )
        self.assertEqual(
            self.mol.molecular_weight,
            1110.300954404405
        )
        self.assertEqual(
            self.mol.monoisotopic_mass,
            1109.12832
        )

    def test_ipmolecule_methods(self):
        for ipmethod in VALID_IPMETHODS:
            for dropmethod in VALID_DROPMETHODS:
                mol = IPMolecule(
                    'Pd2C10H5',
                    ipmethod=ipmethod,
                    dropmethod=dropmethod,
                )
                test = mol.gaussian_isotope_pattern  # test gaussian isotope pattern generation

    def test_ipmolecule(self):
        self.assertEqual(
            self.ipmol.estimated_exact_mass,
            1109.1303706381723,
        )
        self.assertEqual(
            self.ipmol.barip,
            [[1105.130443, 1106.133823749481, 1107.1290292337153, 1108.1305157201678, 1109.1303706381723,
              1110.1328590930914, 1111.1301978511672, 1112.1325950611867, 1113.1318575059308, 1114.134086933976,
              1115.1370272665604, 1116.140052, 1117.143407],
             [2.287794397621507, 1.5228133756325326, 25.476059354316945, 66.8193866193291, 100.0, 52.65050639843156,
              74.88108058795096, 42.5730473226288, 39.36707265932168, 20.17253048748261, 5.990476280101723,
              1.1848920932846654, 0.16082254122736006]]
        )
        self.ipmol - 'PPh3'  # test subtraction
        self.ipmol + 'PPh3'  # test addition
        mol2 = IPMolecule('N(Et)2(CH2(13C)H2(2H))2')
        self.ipmol + mol2  # test class addition


class TestmzML(unittest.TestCase):
    def test_mzml(self):
        mzml = mzML(
            validation_path / 'MultiTest',
            verbose=False
        )
        self.assertEqual(  # check that the correct function keys were pulled
            mzml.functions.keys(),
            {1, 3, 4},
        )

        @mzml.foreachchrom
        def testperchrom(chromatogram):
            attr = branch_attributes(chromatogram)
            return attr['id']

        self.assertEqual(  # test chromatogram decorator
            testperchrom(),
            [u'TIC', u'SRM SIC Q1=200 Q3=100 function=2 offset=0']
        )

        @mzml.foreachscan
        def testperspec(spectrum):
            p = CVParameterSet.create_from_branch(spectrum)
            return p["MS:1000016"].value

        self.assertEqual(  # test spectrum decorator
            testperspec(),
            [0.0171000008, 0.135733336, 0.254333347, 0.372983336, 0.491699994, 0.0510833338, 0.169750005,
             0.288383335, 0.407000005, 0.525833309, 0.0847499967, 0.20341666, 0.322033346, 0.440683335]
        )

        self.assertEqual(  # test intensity summing
            sum(mzml.sum_scans()[1]),
            162804754.0
        )

        self.assertEqual(  # test scan indexing
            sum((mzml[2])[1]),
            6742121
        )

        self.assertEqual(  # test time indexing
            sum((mzml[0.01])[1]),
            56270834
        )


class TestXLSX(unittest.TestCase):
    def test_xlsx(self):
        xlfile = XLSX(
            validation_path / 'xlsx_validation',
            verbose=False
        )
        spec, xunit, yunit = xlfile.pullspectrum('example MS spectrum')
        multispec = xlfile.pullmultispectrum('example multi-spectrum')
        rsimparams = xlfile.pullrsimparams()
        xlout = XLSX(
            validation_path / 'xlsxtestout.xlsx',
            create=True,
            verbose=False
        )
        xlout.writespectrum(spec[0], spec[1], 'test single spectrum out', xunit, yunit)
        for key, val in sorted(multispec.items()):
            xlout.writemultispectrum(
                multispec[key]['x'],
                multispec[key]['y'],
                multispec[key]['xunit'],
                multispec[key]['yunit'],
                'Function Chromatograms',
                key
            )
        xlout.save()
        os.remove(
            validation_path / 'xlsxtestout.xlsx'
        )


class TestSpectrum(unittest.TestCase):
    def test_spectrum(self):
        spec = Spectrum(3)
        spec.add_value(479.1, 1000)
        self.assertEqual(
            spec.trim(),
            [[479.1], [1000]]
        )

        spec2 = Spectrum(3)
        spec2.add_value(443.1, 1000)
        self.assertEqual(
            spec2.trim(),
            [[443.1], [1000]]
        )
        spec += spec2
        self.assertEqual(
            spec.trim(),
            [[443.1, 479.1], [1000, 1000]]
        )
        spec3 = Spectrum(3, start=50, end=2500)
        spec3.add_value(2150.9544, 1000)
        self.assertEqual(
            spec3.trim(),
            [[2150.954], [1000]]
        )
        spec += spec3

        self.assertEqual(
            spec.trim(True),
            [[50.0, 443.1, 479.1, 2150.954, 2500], [0.0, 1000, 1000, 1000, 0.0]]
        )
        spec.end = 2100.
        self.assertEqual(
            spec.trim(),
            [[443.1, 479.1], [1000, 1000]]
        )

    def test_element(self):
        mol = Spectrum(
            3,
            start=0.,
            end=100.,
            filler=0.,
        )
        mol.add_spectrum(  # start with a Cl
            *element_intensity_list('Cl')
        )
        mol.add_element(  # add another Cl
            *element_intensity_list('Cl')
        )
        self.assertEqual(
            mol.trim(),
            [[69.938, 71.935, 73.932], [0.5739577600000001, 0.36728448, 0.05875776]]
        )
        mol.add_element(
            *element_intensity_list('Pd')
        )
        self.assertEqual(
            mol.trim(),
            [[171.843, 173.84, 173.841, 174.842, 175.837, 175.838, 175.841, 176.839, 177.835, 177.838, 177.841, 178.836,
              179.835, 179.838, 179.843, 181.835, 181.84, 183.837],
             [0.005854369152000001, 0.0037463016960000003, 0.06393889446400002, 0.128164767808, 0.000599329152,
              0.040915491072, 0.15686265580800002, 0.082014624384, 0.006545614464, 0.100378848384, 0.15186922329600003,
              0.013120607808, 0.016058495808, 0.097183473408, 0.06726784947200001, 0.015547303296, 0.043045741056,
              0.0068864094719999994]]

        )
        mol.charge = 2
        self.assertEqual(
            mol.trim()[0],
            [85.922, 86.92, 86.921, 87.421, 87.919, 87.919, 87.921, 88.42, 88.918, 88.919, 88.921, 89.418, 89.918,
             89.919, 89.922, 90.918, 90.92, 91.919]

        )
        del mol.charge
        self.assertEqual(
            mol.trim()[0],
            [171.843, 173.84, 173.841, 174.842, 175.837, 175.838, 175.841, 176.839, 177.835, 177.838, 177.841, 178.836,
             179.835, 179.838, 179.843, 181.835, 181.84, 183.837]
        )

    def test_indexing(self):
        """tests calculated indexing for filled Spectrum objects"""
        spec = Spectrum(3, empty=False)
        for i in range(1000):
            num = random()
            mz = num * spec.end
            try:
                index = spec.index(mz)
            except ValueError:
                continue
            self.assertEqual(
                round(mz, 3),
                round(spec.x[index], 3)
            )


if __name__ == '__main__':
    unittest.main(verbosity=2)
