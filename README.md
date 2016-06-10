# mass-spec-python-tools: mass spec made easier

## What is it?
This is a collection of scripts to aid in the processing and interpretation of mass spectrometric data. 

The project was created by Lars Yunker at the University of Victoria, Victoria, B.C.

## Requirements and Installation:
The scripts were written for Python 2.7 but support for 3.x is in development.

## Getting started
* Download the entire repository to a folder and load scripts as needed. 
* Inexperienced python users can edit the input parameters in the files directly. 
* Experienced python users can import many of the classes and scripts directly (e.g. in iPython). Not all scripts are completely callable yet (work in progres). 

For specific instructions on using the PyRSIM script, see [this tutorial video](https://www.youtube.com/watch?v=zc8i54EiCGY)

## Errors
If you encounter an error, please submit an Issue in Github with as much information as possible
Things that are helpful to include:
* the raw file you were trying to parse (zip it first)
* the exact parameters you were using
* any additional files you were supplying to the script
* the error output

## Currently available tools:
### PyRSIM
This script takes supplied raw and parameters files and generates reconstructred single ion monitoring traces
### isotope pattern overlay
Takes a supplied mass spectrum and overlays the predicted isotope pattern onto it
### video frame renderer
A tool for generating a series of images showing mass spectrum and traces which can be combined into a video
### y-axis zoom figure
Renders a series of images which zoom into the y-axis to illustrate the dynamic range of mass spectrometers

## Classes and packages
### tome_v02
A collection of scripts which are required by several of the tools
### _Colour
Interprets a provided colour and converts it into several other formats (mostly this is used for pyplot colour conversion)
### _Molecule
Takes a supplied molecular formula and calculates several physico-chemical properties of that formula
Requires several other dictionaries to function:
* _crc_mass: a dictionary of exact masses and natural abundances used by the Molecule class (obtained from the CRC Handbook of Chemistry and Physics 2015)
* _nist_mass: as above but obtained from the NIST database
* _formabbrvs: a dictionary of common abbreviations in chemical formulas used by the Molecule class
### _ScriptTime
A class for timing python scripts
### _Spectrum
Generates a spectrum which can be added to (very useful for combining spectra that have high-precision x values)
### _XLSX
A class for handling and writing excel files (uses openpyxl)
### _mzML
A class for loading and interpreting mzML files

## Still under development:
### ED-ESI plot
This will eventually take a MS/MS plot and generate an energy dependant mass spectrum plot

## License
these tools are not yet licensed