"""
It builds up the PDB files from the SMARTS tags in the json file.
It also creates a new json file with the pairs of PDB indexes and
SMARTS tags.
"""

import os
import json
from rdkit import Chem
from rdkit.Chem import AllChem


INPUT_FOLDER = "QM"
OUTPUT_FOLDER = "pdb"
INPUT_JSON_NAME = "ids_to_smarts.json"
OUTPUT_JSON_NAME = "pdbs_to_smarts.json"

with open('../{}/{}'.format(INPUT_FOLDER, INPUT_JSON_NAME), 'r') as f:
    tags = json.load(f)

unique_tags = set(tags.values())

os.makedirs("../{}".format(OUTPUT_FOLDER), exist_ok=True)

pdb_to_smarts = dict()

for index, tag in enumerate(unique_tags):
    mol = Chem.MolFromSmiles(tag)
    if mol:
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol)
        Chem.MolToPDBFile(mol, '../{}/{}.pdb'.format(OUTPUT_FOLDER, index))
        pdb_to_smarts[index] = tag

json_output = json.dumps(pdb_to_smarts)

with open("../{}/{}".format(OUTPUT_FOLDER, OUTPUT_JSON_NAME), "w") as f:
    f.write(json_output)

