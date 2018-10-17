from pythoms.molecule import Molecule, IPMolecule, molecular_weight_error
import pickle
import os


kwargs = {'string': 'C61H51IP3Pd', 'dropmethod': 'consolidate'}

mol = IPMolecule(
    verbose=True,
    **kwargs,
)

# mol2 = IPMolecule(
#
#     ** kwargs,
#     ipmethod='combinatorics',
# )

mol3 = IPMolecule(
    **kwargs,
    ipmethod='hybrid',
)

