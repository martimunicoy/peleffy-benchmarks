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
    # Validate both lists, confirm either are empty
    if len(x) == 0 and len(y) > 0:
        return y  # All y values are unique if x is empty
    elif len(y) == 0 and len(x) > 0:
        return x  # All x values are unique if y is empty

    try:
        # Immutable and Unique - Convert list of tuples into set of tuples
        first_set = set(map(tuple, x))
        secnd_set = set(map(tuple, y))

    # Dealing with a 1D dataset (list of items)
    except TypeError:
        # Unique values only
        first_set = set(x)
        secnd_set = set(y)
    commond = []
    for first_element in first_set:
        for second_element in secnd_set:
            if Counter(first_element) == Counter(second_element):
                commond.append(first_element)
    return list(set(first_set) - set(commond)), list(set(secnd_set) - set(commond))


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
        rot = re.sub('_', '', rot)
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
            if (Counter(element_resolution1[1]) == Counter(element_resolution2[1])) and element_resolution1[0] != element_resolution2[0]:
                result = element_resolution1[0], element_resolution2[0], element_resolution1[1]
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
    i = 0  # unused variable 'i'
    for rot in rot_lib:
        rot = re.sub('_', '', rot)
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
