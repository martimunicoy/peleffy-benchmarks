
"""
Exemple of how to use the class MoleculeMinimized to minimize using PELE a ligand [2.pdb] using the PELE executable version
[/home/municoy/builds/PELE/PELE-repo_serial/PELE-1.6]
python test_minimization.py 2.pdb /home/municoy/builds/PELE/PELE-repo_serial/PELE-1.6
"""
import MoleculeMinimized as MM
import os
import argparse


def parse_args():
    """
        Parse command line arguments
        :returns: object -- Object containing command line options
    """
    parser = argparse.ArgumentParser(description="Perform PELE minimizations")
    parser.add_argument("input_file", type=str, help="Input .pdb file.")
    parser.add_argument("PELE_version", type=str, help = "Exacutable PELE.")
    args = parser.parse_args()
    return args.input_file, args.PELE_version


def main(input_file, PELE_version):
	new_molecule = MM.MoleculeMinimized(input_file, PELE_version)
	new_molecule.minimize(input_file = input_file, PELE_version = PELE_version)

if __name__ == "__main__":
    input_file, PELE_version = parse_args()
    main(input_file, PELE_version)

