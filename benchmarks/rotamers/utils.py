"""
This module contains a variety of helpful tools for tests.
"""

import re
from collections import Counter

# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#                                                            Operations with lists
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

def commond(x,y): 
    """Function that gets the commond elements of two lists regardless the order."""
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
    commond_set = set(map(tuple,commond))
    return commond_set


def differentiate(x, y):
    """Function that returns the different elements in two lists"""
    commond_set = commond(x,y)
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

# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#                                                            Bonds
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

def get_bonds(rot_lib):
    """Reads from a file in .pdb/.mae format and returs a list with the bonds in the rotamer library"""
    bonds = []
    for rot in rot_lib: 
        rot = re.sub('_', '', rot)
        rot_list = list(rot.split(" "))
        if len(rot_list) == 8:
             bond = (rot_list[5],rot_list[6])
             bonds.append(bond)
    return bonds

# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#                                                            Resolution
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
def get_resolution(rot_lib):
    """Reads from a file in .pdb/.mae format and returs a list with the resolution for each bond in the rotamer library"""
    resolution_list = []
    bonds = []
    for rot in rot_lib: 
        rot = re.sub('_', '', rot)
        rot_list = list(rot.split(" "))
        if len(rot_list) == 8:
            resolution = rot_list[4]
            bond = (rot_list[5],rot_list[6])
            bonds.append(bond)
            resolution_element = (resolution,bond)
            resolution_list.append(resolution_element)
    return resolution_list, bonds



def get_resolution_elements(pdb_file, mae_file):
    """Gets the resolution for """
    resolution_offpele, bonds_offpele = get_resolution(pdb_file)
    resolution_PlopRotTemp, bonds_PlopRotTemp = get_resolution(mae_file)
    return resolution_offpele, resolution_PlopRotTemp

def compare_resolution(resolution1, resolution2):
    """Compares the resolution from two lists and returns which bonds ave different resolution in the two rotamers libraries."""
    resolution1_set = set(map(tuple,resolution1))
    resolution2_set = set(map(tuple,resolution2))
    results_list = []
    for element_resolution1 in resolution1_set:
        for element_resolution2 in resolution2_set: 
            if (Counter(element_resolution1[1]) == Counter(element_resolution2[1])) and element_resolution1[0] != element_resolution2[0]:
                result = element_resolution1[0],element_resolution2[0], element_resolution1[1]
                results_list.append(result)
    return results_list

# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#                                                            Writing in a file
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
def write_diff_bonds(file,diff1,diff2):
    """Writes into a fill the differences between the to rotamers libraries."""
    file.write('Bonds:' + '\n')
    for i in diff1:
        file.write('+' + str(i) + '\n')
    for j in diff2:
        file.write('-' + str(j) + '\n')

def write_diff_resolution(file, resolution_differences):
    file.write('Resolution:' + '\n')
    for difference in resolution_differences: 
        file.write(str(difference) + '\n')
