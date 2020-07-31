"""
This module contains a variety of helpful tools for tests.
"""

def get_bonds(rot_lib):
	""" Obtains a list of all the pairs of atoms with the rotable bonds.
	"""
	bonds = []
    for rot in rot_lib: 
        #print(line1)
        rot = re.sub('_', '', rot)
        rot_list = list(rot.split(" "))
        if len(rot_list) == 8:
             bond = (rot_list[5],rot_list[6])
             bonds.append(bond)
    print(bonds)
    return bonds

def differentiate(x, y):
    """
    Retrieve a unique of list of elements that do not exist in both x and y.
    Capable of parsing one-dimensional (flat) and two-dimensional (lists of lists) lists.

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
            print(first_element, second_element, first_element == second_element)
            if Counter(first_element) == Counter(second_element):
                commond.append(first_element)
    return list(set(first_set) - set(commond)), list(Counter(set(secnd_set)) - Coundet(set(commond)))

def write_diff(diff1,diff2):
	f.write(diff1)
	f.write(diff2)
