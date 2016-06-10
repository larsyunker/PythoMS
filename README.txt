mass-spec-python-tools is a collection of scripts to aid in the processing and interpretation of mass spectrometric data.

The project was created by Lars Yunker at the University of Victoria, Victoria, B.C.

Currently available tools:
PyRSIM: takes supplied raw and parameters files and generates reconstructred single ion monitoring traces
isotope pattern overlay: takes a supplied mass spectrum and overlays the predicted isotope pattern onto it
video frame renderer: a tool for generating a series of images showing mass spectrum and traces which can be combined into a video
y-axis zoom figure: renders a series of images which zoom into the y-axis to illustrate the dynamic range of mass spectrometers

Tools and classes:
tome_v02: a collection of scripts which are required by several of the tools
_Colour: interprets a provided colour and converts it into several other formats (mostly this is used for pyplot colour conversion)
_Molecule: takes a supplied molecular formula and calculates several physico-chemical properties of that formula
_ScriptTime: a class for timing python scripts
_Spectrum: generates a spectrum which can be added to (very useful for combining spectra that have high-precision x values)
_XLSX: a class for handling and writing excel files (uses openpyxl)
_mzML: a class for loading and interpreting mzML files

_crc_mass: a dictionary of exact masses and natural abundances used by the Molecule class (obtained from the CRC Handbook of Chemistry and Physics 2015)
_nist_mass: as above but obtained from the NIST database
_formabbrvs: a dictionary of common abbreviations in chemical formulas used by the Molecule class


In beta:
ED-ESI plot: this will eventually take a MS/MS plot and generate an energy dependant mass spectrum plot