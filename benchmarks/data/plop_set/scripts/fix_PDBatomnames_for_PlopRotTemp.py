"""
It fixes the PDB atom names found in a PDB file to fix the parsing errors
from PlopRotTemp. Please, note that in order to fix these errors, PDB
atom names must have, if possible, a whitespace in the first position.
"""

import os
from glob import glob


INPUT_FOLDER = "original"
OUTPUT_FOLDER = "fixed_original"


os.makedirs("../{}".format(OUTPUT_FOLDER))

for pdb_file in glob("../{}/*.pdb".format(INPUT_FOLDER)):
    pdb_name = os.path.basename(pdb_file)
    fixed_pdb = ""
    with open(pdb_file) as f:
        for line in f:
            if line.startswith("HETATM"):
                PDB_atom_name = line[12:16]
                if PDB_atom_name[0] != " " and PDB_atom_name[3] == " ":
                    fixed_PDB_atom_name = " " + PDB_atom_name[0:3]
                    line = line[:12] + fixed_PDB_atom_name + line[16:]

            fixed_pdb += line

    with open("../{}/{}".format(OUTPUT_FOLDER, pdb_name), 'w') as f:
        f.write(fixed_pdb)

