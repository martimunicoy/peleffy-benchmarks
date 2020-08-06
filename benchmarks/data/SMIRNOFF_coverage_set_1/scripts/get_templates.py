"""
It creates the PELE templates of all ligands defined as PDB files.
"""

import os
from glob import glob
from offpele.main import run_offpele


INPUT_FOLDER = "pdb"
FORCEFIELD = "openff_unconstrained-1.2.0.offxml"
RESOLUTION = 30
CHARGES_METHOD = "gasteiger"
OUTPUT_FOLDER = "templates"


for pdb_file in glob("../{}/*.pdb".format(INPUT_FOLDER)):
    pdb_name = os.path.basename(pdb_file)
    pdb_name = os.path.splitext(pdb_name)[0]
    os.makedirs("../{}/{}/".format(OUTPUT_FOLDER, pdb_name), exist_ok=True)
    try:
        run_offpele(pdb_file, forcefield=FORCEFIELD, resolution=RESOLUTION,
                    charges_method=CHARGES_METHOD,
                    output="../{}/{}/".format(OUTPUT_FOLDER, pdb_name))
    except:
        print(' - Skipping {}.pdb '.format(pdb_name) +
              'because of an error when generating the template')

