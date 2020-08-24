"""
It computes the Hydration Free Energy for a molecule from a PDB file. 

Example: 
---------

>>>  python hydrationfreeenergy.py ligand.pdb 
"""

import argparse
import os.path 
import re
import os 
import MoleculeMinimized as MM
import DiffEnergies as DiffEnergies


def parse_args():
    """
        Parse command line arguments
        :returns: object -- Object containing command line options
    """
    parser = argparse.ArgumentParser(description="Compute The Hydration Free Energy")
    parser.add_argument("pdb_input_file", type=str, help="PDB for the ligand to minimize and compute the Hydration Free Energy.")
    parser.add_argument("-o", "--output_path", type=str, default="", help="Path where to write the results")
    parser.add_argument("-of", "--output_file", type=str, default="hydration_free_energy.csv", help="File with the energy difference.")
    args = parser.parse_args()
    return args.pdb_input_file, args.output_path, args.output_file

def get_outputs(path, input_file):
    """
    It gets the two PDB files from the minimization for the molecule in vacuum and OBC.
    """

    vacuum_input_file = os.path.join(path, 'VACUUM_minimization.out')
    OBC_input_file = os.path.join(path, 'OBC_minimization.out')
    return vacuum_input_file, OBC_input_file


def main(pdb_input_file, out_path, out_file):

    args = parse_args()
    if out_path != "": out_file = os.path.join(out_path,out_file)
    #by default, if it need to be a parameter this should be modified
    PELE_version = '/home/municoy/builds/PELE/PELE-repo_serial/PELE-1.6'
    #minimize the molecule
    new_ligand = MM.MoleculeMinimized(pdb_input_file, PELE_version)
    new_ligand.minimize(pdb_input_file, PELE_version)
    #computes the Hydration Free Energy
    path = os.getcwd()
    vacuum_input_file, OBC_input_file = get_outputs(path, pdb_input_file)
    name = os.path.splitext(pdb_input_file)[0]
    f = open(os.path.join(os.getcwd(), out_file), 'w')
    DiffEnergies.compare_energies(f, vacuum_input_file, OBC_input_file)
    f.close()

if __name__ == "__main__":
    pdb_input_file, out_path, out_file = parse_args()
    main(pdb_input_file, out_path, out_file)