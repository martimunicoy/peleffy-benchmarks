"""

Example:
-------------
python AlignRMSD.py /home/lauramalo/repos/offpele-benchmarks/benchmarks/data/SMIRNOFF_coverage_set_1/pdb/5.pdb /home/lauramalo/repos/offpele-benchmarks/benchmarks/data/SMIRNOFF_coverage_set_1/pdb/pdbs_to_smarts.json /home/lauramalo/repos/offpele-benchmarks/benchmarks/data/SMIRNOFF_coverage_set_1/QM/ids_to_smarts.json

"""


import os 
import argparse
import re
import json
from pathlib import Path
import subprocess

from rdkit import Chem 
import rdkit.Chem.rdmolfiles

import link_structures as LinkStructures
import MoleculeMinimized as MM



def parse_args():
    """
        Parse command line arguments
        :returns: object -- Object containing command line options
    """
    parser = argparse.ArgumentParser(description="Align and compute RMSD")
    parser.add_argument("pdb_input_file", type=str, help="Original pdb file.")    
    parser.add_argument("pdb_to_smarts", type=str, help="JSON file connecting pdbs to SMARTS")
    parser.add_argument("ids_to_smarts", type=str, help="JSON file connecting ids to SMARTS.")
    parser.add_argument("-o", "--output_path", type=str, default="./output/", help="Path where to write the results")
    parser.add_argument("-of", "--output_file", type=str, default="hydration_free_energy.csv", help="File with the energy difference.")
    args = parser.parse_args()
    return args.pdb_input_file, args.pdb_to_smarts, args.ids_to_smarts, args.output_path, args.output_file



def get_xyz_from_pdb(pdbfile):
    """
        It generates a coordinates file (.xyz) from a PDB file. The coordinates file will be saved as molecule.xyz.
        The function returns the path of the new generated coordinates file. 
    """
    p = Path(pdbfile)
    path = p.parents[0]
    mol = rdkit.Chem.rdmolfiles.MolFromPDBFile(pdbfile, removeHs = False)
    rdkit.Chem.rdmolfiles.MolToXYZFile(mol, os.path.join(path, 'molecule.xyz'))
    return(os.path.join(path, 'molecule.xyz'))

def get_comformations(inputpdb, pdb_to_smarts, ids_to_smarts):
    """
        It gets all the corformations from the /QM/ database corresponding to the inputpdb molecule.  
    """
    p = Path(inputpdb)
    name_mol = p.name
    return LinkStructures.main(name_mol, pdb_to_smarts, ids_to_smarts)

def get_molecule_minimized(mol_pdb, solvent):
    """
        It minimized the molecule using PELE's minimization. 
        Depending on the solvent parameter it returns the path to the PDB of the minimized molecule. 
    """
    new_molecule = MM.MoleculeMinimized(mol_pdb, '/home/municoy/builds/PELE/PELE-repo_serial/PELE-1.6')
    new_molecule.minimize(input_file = mol_pdb, PELE_version = '/home/municoy/builds/PELE/PELE-repo_serial/PELE-1.6')
    p = Path(mol_pdb)
    file, folder = p.name, p.parents[0]
    name_folder = os.path.splitext(file)[0] 
    if solvent == False:
        path = os.path.join(os.getcwd(),name_folder,'VACUUM_minimized.pdb')
    if solvent == True:
        path = os.path.join(os.getcwd(),name_folder,'OBC_minimized.pdb')
    return path

def compute_rmsd(mol1,mol2):
    """
        It aligns the to molecules and computes RMSD. Returns the RMSD value. 
        Documentation: https://github.com/charnley/rmsd
        (This implementation could be done importing the package rmsd.)
    """
    rmsd = os.popen("python calculate_rmsd --reorder %s %s" %(mol1,mol2)).readlines()
    return re.sub('\n', '', rmsd[0])


def main(pdb_input_file,pdb_to_smarts, ids_to_smarts, out_path, out_file):

    args = parse_args()
    if out_path != "": out_file = os.path.join(out_path,out_file)

    # USes PELE's minimization.
    mol1 = get_molecule_minimized(pdb_input_file, solvent = False)

    # Saves a coordinates file (.xyz) for the minimazed molecule.
    mol1_path = get_xyz_from_pdb(mol1)

    # Gets all the matching componds for the input molecule
    mol2_list = get_comformations(pdb_input_file, pdb_to_smarts, ids_to_smarts)

    # Performs RMSD with all the conformations and obtains the best one.
    results = []
    name_list = []
    for mol2_path in mol2_list:
        rmsd = compute_rmsd(mol1_path,mol2_path)
        results.append(rmsd)
        comp = Path(mol2_path)
        name_list.append(comp.name)
    dic = dict(zip(name_list, results))
    min_val = min(dic.values())
    a = ([k for k,v in dic.items() if v == min_val])

    #Writes out the obtained results. 
    print('The minimum RMSD is %s for the %s conformation.' %(str(min_val), a[0]))


if __name__ == "__main__":
    pdb_input_file, pdb_to_smarts, ids_to_smarts, out_path, out_file = parse_args()
    main(pdb_input_file, pdb_to_smarts, ids_to_smarts, out_path, out_file)
