
import argparse
import os.path as os
import re
import json


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
    Parsing the date from s json file into a dictionary.
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
	PATH_FOLDER = 'QM'
	path_list = []
	for id in ids:
		path = os.join('..',PATH_FOLDER, id + '.xyz')
		path_list.append(path)
	print(path_list)
	return path_list

def main(pdb_file, pdb_to_smarts, ids_to_smarts):
	args = parse_args()
	label = pdb_file.replace('.pdb','')
	dict_pdb = parse_json_file(pdb_to_smarts)
	dict_ids = parse_json_file(ids_to_smarts)
	ids = get_ids(dict_ids,dict_pdb.get(label))
	return get_path(ids)

if __name__ == "__main__":
    pdb_file, pdb_to_smarts, ids_to_smarts = parse_args()
    main(pdb_file, pdb_to_smarts, ids_to_smarts)
