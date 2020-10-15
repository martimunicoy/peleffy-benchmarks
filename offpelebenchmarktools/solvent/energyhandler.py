"""
It compares the output energies after the minimization using the
OpenForceField toolkit with PELE considering a solvent (OBC) and in vacuum,
and it computes the Hydration Free Energy.
"""


def parse_file(file):
    """
    Parsing the date from the input file.
    """
    import re

    with open(file, 'r') as f:
        data = list()
        group = dict()
        for key, value in re.findall(r'(.*):\s*([\dE+-.]+)', f.read()):
            if key in group:
                data.append(group)
                group = dict()
            group[key] = value
        data.append(group)
    return data


def get_energy(file):
    """
    Given a file, it does a data parsing of the file and creates a dictionary
    of the energies after minimization. Returns the total energy computed
    after the minimization.
    """
    data = parse_file(file)
    energies_minimized = data[-1]
    label = list(energies_minimized.keys())[-1]
    return (energies_minimized.get(label))


def compute_energies(vacuum_file, OBC_file):
    """
    Given two files, it returns the vacuum energy, the OBC energy and the
    difference energy.
    """
    vacuum_energy = get_energy(vacuum_file)
    OBC_energy = get_energy(OBC_file)
    hydration_energy = float(OBC_energy) - float(vacuum_energy)
    return vacuum_energy, OBC_energy, hydration_energy
