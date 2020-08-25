
"""

It compares the output energies after the minimization using the OpenForceField toolkit with PELE 
considering a solvent (OBC) and in vacuum, and it computes the Hydration Free Energy. 

Example: 
-------------

python /home/lauramalo/tests/solvent/DiffEnergies.py VACUUM_minimization.out OBC_minimization.out
"""


import argparse
import os.path as os
import re


def parse_args():
    """
        Parse command line arguments
        :returns: object -- Object containing command line options
    """
    parser = argparse.ArgumentParser(description="Compute The Hydration Free Energy")
    parser.add_argument("vacuum_input_file", type=str, help="Output file of minimization in VACUUM.")
    parser.add_argument("OBC_input_file", type=str, help="Output file of minimization in OBC.")
    parser.add_argument("-o", "--output_path", type=str, default="", help="Path where to write the results")
    parser.add_argument("-of", "--output_file", type=str, default="hydration_free_energy.csv", help="File with the energy difference.")
    args = parser.parse_args()
    return args.vacuum_input_file,args.OBC_input_file, args.output_path, args.output_file


def parse_file(file):
    """
    Parsing the date from the input file.
    """
    with open(file, 'r') as f:
        data  = list()
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
    Given a file, first it does a data parsing of the file and creates a dictionary of the nergies after minimization. 
    Returns the total energy computed after the minimization.
    """
    data = parse_file(file)
    energies_minimized = data[-1]
    label = list(energies_minimized.keys())[-1]
    return (energies_minimized.get(label))

def compute_energies(vacuum_file, OBC_file):
    """
    Given two files, computes the minimized energy for the two cases and returns the
    difference written in the out_file.
    """
    vacuum_energy = get_energy(vacuum_file)
    OBC_energy = get_energy(OBC_file)
    hydration_energy = float(OBC_energy) - float(vacuum_energy)
    return vacuum_energy, OBC_energy, hydration_energy
    

def write_energies(f, vacuum_energy, OBC_energy, hydration_energy):
    f.write('Vacuum_energy:' + vacuum_energy + '\n')
    f.write('OBC energy:' +  OBC_energy + '\n')
    f.write('Hydration free energy:' + str(hydration_energy) + '\n')


def compare_energies(f, vacuum_input_file, OBC_input_file):
    vacuum_energy, OBC_energy, hydration_energy = compute_energies(vacuum_input_file, OBC_input_file)
    write_energies(f, vacuum_energy, OBC_energy, hydration_energy)


def main(vacuum_input_file, OBC_input_file, out_path, out_file):

    args = parse_args()
    if out_path != "": out_file = os.join(out_path,out_file)
    f = open(out_file, 'w')
    compare_energies(f, vacuum_input_file, OBC_input_file)
    f.close()


if __name__ == "__main__":
    vacuum_input_file, OBC_input_file, out_path, out_file = parse_args()
    main(vacuum_input_file, OBC_input_file, out_path, out_file)
