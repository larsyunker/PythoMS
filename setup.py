from setuptools import setup, find_packages
import pythoms

NAME = 'pythoms'
AUTHOR = 'Lars Yunker'

PACKAGES = find_packages()
KEYWORDS = ', '.join([
    'mass spectrometry',
    'mass spec',
    'mzML',
    'isotope pattern',
    'HUPO PSI-MS',
])

with open('README.MD') as f:
    long_description = f.read()

setup(
    name=NAME,
    version=pythoms.__version__,
    description='A Python library to aid in the processing and interpretation of mass spectrometric data.',
    long_description=long_description,
    long_description_content_type='text/markdown',
    author=AUTHOR,
    url='https://github.com/larsyunker/PythoMS',
    packages=PACKAGES,
    license='MIT',
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Development Status :: 5 - Production/Stable',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Chemistry',
        'Topic :: Scientific/Engineering :: Information Analysis',
        'Operating System :: OS Independent',
        'Natural Language :: English'
    ],
    install_requires=[
        'numpy>=1.14.2',
        'openpyxl>=2.5.2',
        'matplotlib>=2.1.2',
        'scipy>=1.1.0',
        'sympy>=1.1.1',
        'obonet==0.2.5',  # they changed attribute names without deprecationwarnings, so only this version is verified
        'isospecpy>=2.0.2',
    ],
    keywords=KEYWORDS,
)
