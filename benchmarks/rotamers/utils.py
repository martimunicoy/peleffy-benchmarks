"""
This module contains a variety of helpful tools for tests.
"""

import re
from collections import Counter


# ---------------------
# Operations with lists
# ---------------------
def commond(x, y):
    """
    Function that gets the commond elements of two lists regardless
    the order.
    """
    # Validate both lists, confirm either are empty
    if len(x) == 0 and len(y) > 0:
        return y  # All y values are unique if x is empty
    elif len(y) == 0 and len(x) > 0:
        return x  # All x values are unique if y is empty

    first_set = set(map(tuple, x))
    secnd_set = set(map(tuple, y))
    commond = []
    for first_element in first_set:
        for second_element in secnd_set:
            if Counter(first_element) == Counter(second_element):
                commond.append(first_element)
    commond_set = set(map(tuple, commond))
    return commond_set


def differentiate(x, y):
    """
    Retrieve a unique of list of elements that do not exist in both x and y.
    Capable of parsing one-dimensional (flat) and two-dimensional
    (lists of lists) lists.

    :param x: list #1
    :param y: list #2
    :return: list of unique values
    """
    commond_set = commond(x, y)
    first_set = set(x)
    secnd_set = set(y)
    diff1 = list(first_set - commond_set)
    diff2 = list(secnd_set - commond_set)
    for diff2_element in diff2:
        for commond_element in commond_set:
            if Counter(diff2_element) == Counter(commond_element):
                diff2.remove(diff2_element)
    for diff1_element in diff1:
        for commond_element in commond_set:
            if Counter(diff1_element) == Counter(commond_element):
                diff1.remove(diff1_element)

    return diff1, diff2


# -----
# Bonds
# -----
def get_bonds(rot_lib):
    """
    Reads from a file in .pdb/.mae format and returs a list with the
    bonds in the rotamer library
    """
    bonds = []
    for rot in rot_lib:
        rot_list = list(rot.split(" "))
        if len(rot_list) == 8:
            bond = (rot_list[5], rot_list[6])
            bonds.append(bond)
    return bonds


# ----------
# Resolution
# ----------
def get_resolution(rot_lib):
    """
    Reads from a file in .pdb/.mae format and returs a list with the
    resolution for each bond in the rotamer library
    """
    resolution_list = []
    bonds = []
    for rot in rot_lib:
        rot = re.sub('_', '', rot)
        rot_list = list(rot.split(" "))
        if len(rot_list) == 8:
            resolution = rot_list[4]
            bond = (rot_list[5], rot_list[6])
            bonds.append(bond)
            resolution_element = (resolution, bond)
            resolution_list.append(resolution_element)
    return resolution_list, bonds


def get_resolution_elements(pdb_file, mae_file):
    """
    Gets the resolution from both rotamer libraries (PlopRotTemp and
    offpele).
    """
    resolution_offpele, bonds_offpele = get_resolution(pdb_file)
    resolution_PlopRotTemp, bonds_PlopRotTemp = get_resolution(mae_file)
    return resolution_offpele, resolution_PlopRotTemp


def compare_resolution(resolution1, resolution2):
    """
    Compares the resolution from two lists and returns which bonds have
    different resolution in the two rotamers libraries.
    """
    resolution1_set = set(map(tuple, resolution1))
    resolution2_set = set(map(tuple, resolution2))
    results_list = []
    for element_resolution1 in resolution1_set:
        for element_resolution2 in resolution2_set:
            if ((Counter(element_resolution1[1])
                    == Counter(element_resolution2[1]))
                    and element_resolution1[0] != element_resolution2[0]):
                result = element_resolution1[0], element_resolution2[0], \
                    element_resolution1[1]
                results_list.append(result)
    return results_list


# ------
# Groups
# ------
def get_number_group(rot_lib):
    """
    To do
    """
    text_string = rot_lib.read()
    pattern = 'newgrp'
    return text_string.count(pattern) + 1


def get_groups(rot_lib):
    """
    To do
    """
    groups = []
    group = []
    for rot in rot_lib:
        rot_list = list(rot.split(" "))
        if 'newgrp' in rot_list:
            groups.append(group)
            group = []
        else:
            if len(rot_list) == 8:
                bond = (rot_list[5], rot_list[6])
                group.append(bond)
    groups.append(group)
    return groups


def compare_bonds_groups(groups_offpele, groups_PlopRotTemp):
    """
    To do
    """
    results = []
    result = []
    for group1 in groups_offpele:
        for group2 in groups_PlopRotTemp:
            result.append(list(commond(group1, group2)))
            results.append(result)
            result = []
    print('RESULTS', results)
    print('Group1:', len(groups_offpele[0]), 'Group2:', len(groups_PlopRotTemp[0]), 'Common:', len(results[0][0]))


def compare_groups(groups_offpele, groups_PlopRotTemp):
    """
    To do
    """
    num_offpele = len(groups_offpele)
    num_PlopRotTemp = len(groups_PlopRotTemp)
    result = []
    if num_offpele == num_PlopRotTemp:
        result.append('The number of groups is the same.')
        print('The number of groups is the same.')
        compare_bonds_groups(groups_offpele, groups_PlopRotTemp)
    else:
        result.append('The number of groups is not the same.')
        print('The number of groups is not the same.')


# -----------------
# Writing in a file
# -----------------
def write_diff_bonds(file, diff1, diff2):
    """
    Writes into a fill the differences between the to rotamers
    libraries.
    """
    file.write('Bonds:' + '\n')
    for i in diff1:
        file.write('+' + str(i) + '\n')


def write_diff_resolution(file, resolution_differences):
    """
    To do
    """
    file.write('Resolution:' + '\n')
    for difference in resolution_differences:
        file.write(str(difference) + '\n')


# ------------
# Reading file
# ------------
def extract_rotamer_info_from(path_to_templates, relative_path_to_pdbs):
    """
    It searches and returns a collection of data useful for the analysis
    of rotamer libraries.

    Parameters
    ----------
    path_to_templates : str
        The pattern to all the templates to analyze
    relative_path_to_pbds : str
        The relative path from templates from which pdbs can be retrieved

    Returns
    -------
    mol_ids : list[int]
        The ordered list of molecule ids
    mol_names : list[str]
        The ordered list of molecule names
    groups : list[list[list]]
        The ordered list of bonds defined in each rotamer's library group
        (or branch)
    locations : list[dict]
        The ordered list of dictionaries that link a PDB atom name (key)
        with its location in molecule (value). The location is represented
        by a single char which can be either 'M', if it is in the core,
        or 'S', otherwise
    pdb_path : list[str]
        The ordered list of paths to the pdb
    """
    from glob import glob
    from pathlib import Path

    mol_ids = list()
    mol_names = list()
    groups = list()
    locations = list()
    pdbs = list()

    for template_path in glob(path_to_templates):
        template_path = Path(template_path)
        with open(template_path, 'r') as template:
            try:
                molecule_id = int(template_path.parent.name)
                molecule_name = template_path.name.split('.')[0]
            except ValueError:
                continue

            groups_to_append = get_groups(template)

        parameters_path = template_path.parent.joinpath(
            molecule_name.lower() + 'z')

        try:
            with open(parameters_path, 'r') as template:
                locations.append(get_atom_locations(template))
        except FileNotFoundError:
            continue

        pdb_path = template_path.parent.joinpath(
            relative_path_to_pdbs, '{}.pdb'.format(molecule_id))

        if not pdb_path.is_file():
            continue

        groups.append(groups_to_append)

        mol_ids.append(molecule_id)
        mol_names.append(molecule_name)

        pdbs.append(str(pdb_path))

    assert len(mol_ids) == len(mol_names) \
        and len(mol_ids) == len(groups) \
        and len(mol_ids) == len(locations) \
        and len(mol_ids) == len(pdbs), 'Extracted data has wrong size'

    return mol_ids, mol_names, groups, locations, pdbs


# -------------
# Impact parser
# -------------
def get_atom_locations(parameters_path):
    """
    It links each PDB atom name with its location inside the molecule.
    Whether it is in the core (M) or in a branch (S).

    Parameters
    ----------
    parameters_path : _io.TextIOWrapper
        The TextIOWrapper obtained from an Impact file.

    Returns
    -------
    locations : dict
        The dictionary that links a PDB atom name (key) with its location
        in molecule (value). The location is represented by a single
        char which can be either 'M', if it is in the core, or 'S',
        otherwise
    """
    locations = dict()

    for line in parameters_path:
        # Ignore header lines
        if line.startswith('*'):
            continue

        # The section we are interesed lies before NBON
        if line.startswith('NBON'):
            break

        fields = line.split()

        # We are looking for lines with 9 columns delimited by white spaces
        if len(fields) != 9:
            continue

        # PDB atom names are in the fifth column
        PDB_atom_name = fields[4]

        # Locations are in the third column
        location = fields[2]

        locations[PDB_atom_name] = location

    return locations


# -------------
# Draw rotamers
# -------------
def draw_rotamers(pdb_path, mol_name, groups, locations):
    """
    To do
    """

    COLORS = [(82 / 255, 215 / 255, 255 / 255), (255 / 255, 154 / 255, 71 / 255),
              (161 / 255, 255 / 255, 102 / 255), (255 / 255, 173 / 255, 209 / 255),
              (154 / 255, 92 / 255, 255 / 255), (66 / 255, 255 / 255, 167 / 255),
              (251 / 255, 255 / 255, 17 / 255)]

    from rdkit import Chem
    from rdkit.Chem.Draw import rdMolDraw2D

    mol = Chem.MolFromPDBFile(str(pdb_path), removeHs=False)

    # Hydrogens do not need to be displayed
    mol = Chem.RemoveHs(mol)

    bond_indexes = list()
    bond_color_dict = dict()

    for bond in mol.GetBonds():
        atom1_name = bond.GetBeginAtom().GetPDBResidueInfo().GetName().replace(' ', '_')
        atom2_name = bond.GetEndAtom().GetPDBResidueInfo().GetName().replace(' ', '_')

        bond_id = (atom1_name, atom2_name)
        bond_id_inverse = (atom2_name, atom1_name)

        for color_index, group in enumerate(groups):
            if bond_id in group or bond_id_inverse in group:
                bond_indexes.append(bond.GetIdx())
                try:
                    bond_color_dict[bond.GetIdx()] = COLORS[color_index]
                except IndexError:
                    bond_color_dict[bond.GetIdx()] = (99 / 255,
                                                      122 / 255,
                                                      126 / 255)
                break

    atom_indexes = list()
    radii_dict = dict()
    atom_color_dict = dict()

    for atom in mol.GetAtoms():
        atom_name = atom.GetPDBResidueInfo().GetName().replace(' ', '_')
        if locations[atom_name] == 'M':
            atom_indexes.append(atom.GetIdx())
            radii_dict[atom.GetIdx()] = 0.6
            atom_color_dict[atom.GetIdx()] = (255 / 255, 243 / 255, 133 / 255)

    Chem.rdDepictor.Compute2DCoords(mol)
    draw = rdMolDraw2D.MolDraw2DSVG(500, 500)
    draw.SetLineWidth(4)
    rdMolDraw2D.PrepareAndDrawMolecule(draw, mol, mol_name,
                                       highlightAtoms=atom_indexes,
                                       highlightAtomRadii=radii_dict,
                                       highlightAtomColors=atom_color_dict,
                                       highlightBonds=bond_indexes,
                                       highlightBondColors=bond_color_dict)
    draw.FinishDrawing()

    from IPython.display import SVG

    image = SVG(draw.GetDrawingText())

    return image
