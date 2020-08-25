import os
import pandas as pd
from offpele.topology import Molecule, RotamerLibrary
from offpele.template import Impact
from offpele.solvent import OBC2


PATH_TO_FREESOLV_DATABASE = '../data/FreeSolv/FreeSolv0.52.txt'
OFF_FORCEFIELD = 'openff_unconstrained-1.2.0.offxml'
CHARGES_METHOD = 'gasteiger'

freesolv_db = pd.read_csv(PATH_TO_FREESOLV_DATABASE, delimiter=';',
                          skipinitialspace=True, skiprows=[0, 1, 2], header=0,
                          names=['compound id', 'SMILES', 'iupac name',
                                 'experimental value',
                                 'experimental uncertainty',
                                 'calculated value (GAFF)',
                                 'calculated uncertainty',
                                 'experimental reference',
                                 'calculated reference',
                                 'notes'])

compound_ids = freesolv_db['compound id'].to_list()
smiles_tags = freesolv_db['SMILES'].to_list()

for cid, tag in zip(compound_ids, smiles_tags):
    try:
        # Create dir
        os.makedirs('out/{}'.format(cid), exist_ok=True)

        # Create representation of a particular molecule
        molecule = Molecule(smiles=tag)
        molecule.set_name('LIG')

        # Generate its rotamer library
        rotamer_library = RotamerLibrary(molecule)
        rotamer_library.to_file('out/{}/{}.rot.assign'.format(
            cid, molecule.name.upper()))

        # Generate its parameters and template file
        molecule.parameterize(OFF_FORCEFIELD, charges_method=CHARGES_METHOD)
        impact = Impact(molecule)
        impact.write('out/{}/{}z'.format(cid, molecule.name.lower()))

        # Generate its solvent parameters
        solvent = OBC2(molecule)
        solvent.to_json_file('out/{}/ligandParams.txt'.format(cid))

        # Write molecule's PDB
        print('out/{}/{}.pdb'.format(cid, molecule.name.lower()))
        molecule.to_pdb_file('out/{}/{}.pdb'.format(cid, molecule.name.lower()))

    except Exception:
        print('Skipping:', cid, '-', tag)
        continue
