# mzml and psims modules

The `mzml` module contains a variety of methods and classes related to interpreting mzML files. Many of the
methods are wrappers for the `xml.dom.minidom` package, and are specific to mzML files.

The `psims` module contains several methods and classes for interpreting and interacting with
HUPO-PSI-MS controlled variable definitions. These are used in mzML files, and rather than requiring
the user to look up any undefined values, these tools were created for easy access to controlled
variables and their definitions.

# mzml module

todo


# psims module

This module was designed to all users to easily retrieve a Python-object version of a HUPO-PSI-MS
obo file. The module will first search the current working directory for an *.obo file (any located
obo files will first be validated to ensure that they are HUPO-PSI-MS obo cv files), and
if it fails to find one, the module will attempt to download the correct obo directly from the
HUPO-PSI Github repository.

The `CVParam` class is designed to define and describe a single controlled variable (CV) parameter
and its attribute. The `CVParameterSet` class holds a collection of `CVParam` objects, and the
`CVParameterDefinitions` is a special subclass of `CVParameterSet` which is loaded directly from an
obo file.

```
>>> from pythoms.psims import cv_param_def  # automatic instance of the CVParameterDefinition class
>>> param = cv_param_def['MS:1000073']  # retrieve a parameter using item retrieval
>>> param  # the retrieved parameter is a CVParam instance
CVParam(MS:1000073)
>>> param.name  # name of the parameter
'electrospray ionization'
>>> param.definition  # definition
'"A process in which ionized species in the gas phase are produced from an analyte-containing solution via
highly charged fine droplets, by means of spraying the solution from a narrow-bore needle tip at atmospheric
pressure in the presence of a high electric field. When a pressurized gas is used to aid in the formation of
a stable spray, the term pneumatically assisted electrospray ionization is used. The term ion spray is not
recommended." [PSI:MS]'
```

The `CVParameterSet` is used extensively in the `mzML` class to allow for straightforward access
to attributes, names, and values defined in the mzML files. The `CVParam` class attribute retrievals are
structured to first look whether an attribute value was defined, and if it is not, the attribute's value
will be pulled from the definitions. In this way, even an undefined value may still be accessed because
it is defined in the obo file.

```
>>> from pythoms.psims import CVParam
>>> test_param = CVParam(id='MS:1000073')  # define a CV paramter with only an ID (no specified attributes)
>>> test_param.name  # the obo-defined attributes are accessed
'electrospray ionization'
```
