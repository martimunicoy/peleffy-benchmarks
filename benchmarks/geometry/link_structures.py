"""
Script that links PDB ligant structure which all the corresponding QM conformations of the minimized ligand in the
SMIRNOFF database.
"""
import argparse
import os
import re
import json
from pathlib import Path

def parse_args():
    """
        Parse command line arguments
        :returns: object -- Object containing command line options
    """
    parser = argparse.ArgumentParser(description="Link each minimized structure with the corresponding QM conformations")
    parser.add_argument("pdb_input_file", type=str, help="Original pdb file to be linked.")
    parser.add_argument("pdb_to_smarts", type=str, help="JSON file connecting pdbs to SMARTS")
    parser.add_argument("ids_to_smarts", type=str, help="JSON file connecting ids to SMARTS.")
    args = parser.parse_args()
    return args.pdb_input_file,args.pdb_to_smarts, args.ids_to_smarts


def parse_json_file(file):
    """
    Parsing the date from a json file into a dictionary.
    """
    with open(file, 'r') as json_file:
    	data = json.load(json_file)
    return data

def get_ids(mydict,smarts_label):
	"""
	Given a dictionary and a value for an item returns all the keys with that value.
	"""
	items = mydict.items()
	ids_list = []
	for ids,smarts in items:
		if smarts == smarts_label: ids_list.append(ids)
	return ids_list

def get_path(ids):
	"""
	Returns a list of the paths of the corresponding conformations in /QM/ of the given ligand.
	"""
	PATH_FOLDER = 'QM'
	path_list = []
	for id in ids:
		path = os.path.join('/home/lauramalo/repos/peleffy-benchmarks/benchmarks/data/SMIRNOFF_coverage_set_1',PATH_FOLDER, id + '.xyz')
		path_list.append(path)
	return path_list

def main(pdb_file, pdb_to_smarts, ids_to_smarts):
	"""
	Command line:
	----------
	>>> python link_structures.py [ligand.pdb] [pdbs_to_smarts.json] [ids_to_smarts.json]

	Example:
	----------
	>>> python link_structures.py 2.pdb pdbs_to_smarts.json ids_to_smarts.json

	PATH of the json files:
	----------
	pdbs_to_smarts.json -> /home/lauramalo/repos/peleffy-benchmarks/benchmarks/data/SMIRNOFF_coverage_set_1/pdb/pdbs_to_smarts.json

	ids_to_smarts.json -> /home/lauramalo/repos/peleffy-benchmarks/benchmarks/data/SMIRNOFF_coverage_set_1/ids/ids_to_smarts.json

	"""
	args = parse_args()
	p = Path(pdb_file)
	label = p.name
	label = label.replace('.pdb','')

	dict_pdb = parse_json_file(pdb_to_smarts)
	dict_ids = parse_json_file(ids_to_smarts)
	ids = get_ids(dict_ids,dict_pdb.get(label))
	return get_path(ids)

if __name__ == "__main__":
    pdb_file, pdb_to_smarts, ids_to_smarts = parse_args()
    main(pdb_file, pdb_to_smarts, ids_to_smarts)
