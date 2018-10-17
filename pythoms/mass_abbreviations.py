"""
Dictionary of common abbreviations used in chemical formulas
A new 'element' is triggered by an upper case letter in the MolecularFormula class
Additional entries can be added by using a unique trigger starting with an upper case letter and followed by lower case
letters.
Each new entry should be a dictionary giving the molecular formula
"""

abbrvs = {
    # common chemical abbreviations
    'Ac': {'C': 2, 'H': 3, 'O': 1},  # acetyl
    'Ar+': {'C': 25, 'H': 21, 'P': 1},  # phosphonium tagged aryl group
    'Bn': {'C': 7, 'H': 7},  # benzyl
    'Bu': {'C': 4, 'H': 8},  # butyl
    'Cod': {'C': 8, 'H': 12},  # cyclooctadiene
    'Cp': {'C': 5, 'H': 5},  # cyclopenadienyl
    'Cp*': {'C': 10, 'H': 15},  # pentamethylcyclopentadienyl
    'Cy': {'C': 6, 'H': 11},  # cyclohexyl
    'Et': {'C': 2, 'H': 4},  # ethyl
    'Dba': {'C': 17, 'H': 14, 'O': 1},  # dibenzylideneacetone
    'Dbu': {'C': 9, 'H': 16, 'N': 2},  # 1,8-diazabicycloundec-7-ene
    'Dmf': {'C': 3, 'H': 7, 'N': 1, 'O': 1},  # dimethylformamide
    'Dppm': {'C': 25, 'H': 22, 'P': 2},  # bis(diphenylphosphino)methane
    'Dppe': {'C': 26, 'H': 24, 'P': 2},  # bis(diphenylphosphino)ethane
    'Dppp': {'C': 27, 'H': 26, 'P': 2},  # bis(diphenylphosphino)propane
    'L': {'C': 18, 'H': 15, 'P': 1},  # triphenylphosphine
    'Im+': {'C': 11, 'H': 12, 'N': 2},  # imidazolium tagged aryl group
    'Me': {'C': 1, 'H': 3},  # methyl
    'Mes': {'C': 9, 'H': 12},  # mesitylene
    'Ph': {'C': 6, 'H': 5},  # phenyl
    'Pr': {'C': 3, 'H': 6},  # propyl
    'Py': {'C': 5, 'H': 5, 'N': 1},  # pyridine
    'Tf': {'C': 1, 'F': 3, 'O': 2, 'S': 1},  # trifluoromethanesulfonyl
    'Thf': {'C': 4, 'H': 8, 'O': 1},  # tetrahydrofuran
    'Tol': {'C': 7, 'H': 7},  # tolyl
    'Tr': {'C': 19, 'H': 15},  # triphenylmethyl

    # amino acids (assumes in polypeptide form; subtracting 1 O and 2 H from each)
    'Ala': {'C': 3, 'H': 5, 'N': 1, 'O': 1},  # alanine
    'Arg': {'C': 6, 'H': 11, 'N': 4, 'O': 1},  # arginine
    'Asn': {'C': 4, 'H': 6, 'N': 2, 'O': 2},  # asparagine
    'Asp': {'C': 4, 'H': 5, 'N': 1, 'O': 3},  # aspartic acid
    'Cys': {'C': 3, 'H': 5, 'N': 1, 'O': 1, 'S': 1},  # cysteine
    'Gln': {'C': 5, 'H': 8, 'N': 2, 'O': 1},  # glutamine
    'Glu': {'C': 5, 'H': 7, 'N': 1, 'O': 3},  # glutamic acid
    'Gly': {'C': 2, 'H': 3, 'N': 1, 'O': 1},  # glycine
    'His': {'C': 6, 'H': 7, 'N': 3, 'O': 1},  # histidine
    'Ile': {'C': 6, 'H': 11, 'N': 1, 'O': 1},  # isoleucine
    'Leu': {'C': 6, 'H': 11, 'N': 1, 'O': 1},  # leucine
    'Lys': {'C': 6, 'H': 12, 'N': 2, 'O': 1},  # lysine
    'Met': {'C': 5, 'H': 9, 'N': 1, 'O': 1, 'S': 1},  # methionine
    'Phe': {'C': 9, 'H': 9, 'N': 1, 'O': 1},  # phenylalanine
    'Pro': {'C': 5, 'H': 7, 'N': 1, 'O': 1},  # proline
    'Ser': {'C': 3, 'H': 5, 'N': 1, 'O': 2},  # serine
    'Thr': {'C': 4, 'H': 7, 'N': 1, 'O': 2},  # threonine
    'Trp': {'C': 11, 'H': 10, 'N': 2, 'O': 1},  # tryptophan
    'Tyr': {'C': 9, 'H': 9, 'N': 1, 'O': 2},  # tryosine
    'Val': {'C': 5, 'H': 9, 'N': 1, 'O': 1},  # valine
}
