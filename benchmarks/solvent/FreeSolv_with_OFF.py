import os
import pandas as pd
from offpele.topology import Molecule, RotamerLibrary
from offpele.template import Impact
from offpele.solvent import OBC2
from offpele.main import handle_output_paths
import DiffEnergies as DiffEnergies


PATH_TO_FREESOLV_DATABASE = '/home/lauramalo/repos/offpele-benchmarks/benchmarks/data/FreeSolv/FreeSolv0.52.txt'
OFF_FORCEFIELD = 'openff_unconstrained-1.2.0.offxml'
CHARGES_METHOD = 'gasteiger'
VACUUM_CF = '/home/lauramalo/repos/offpele-benchmarks/benchmarks/solvent/Conf/VACUUM_minimization.conf'
OBC_CF = '/home/lauramalo/repos/offpele-benchmarks/benchmarks/solvent/Conf/OBC_minimization.conf'
PELE_EXEC = '/home/municoy/builds/PELE/PELE-repo_serial/PELE-1.6'
PELE_SRC = '/home/municoy/repos/PELE-repo/'


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
experimental_v = freesolv_db['experimental value'].to_list()

energies = list()

for cid, tag, exp_v in zip(compound_ids, smiles_tags, experimental_v):
    try: 
        # Create dir
        os.makedirs('out/{}'.format(cid), exist_ok=True)

        # Create representation of a particular molecule
        molecule = Molecule(smiles=tag)
        molecule.set_name('UNL')

        # Handle paths
        rotamer_library_output_path, impact_output_path, solvent_output_path = \
            handle_output_paths(molecule=molecule, output='out/{}'.format(cid), as_datalocal=True)

        # Generate its rotamer library
        rotamer_library = RotamerLibrary(molecule)
        rotamer_library.to_file(rotamer_library_output_path)

        # Generate its parameters and template file
        molecule.parameterize(OFF_FORCEFIELD, charges_method=CHARGES_METHOD)
        impact = Impact(molecule)
        impact.write(impact_output_path)

        # Generate its solvent parameters
        solvent = OBC2(molecule)
        solvent.to_json_file(solvent_output_path)

        # Write molecule's PDB
        molecule.to_pdb_file('out/{}/ligand.pdb'.format(cid))

        # Link folders
        os.symlink('{}Data'.format(PELE_SRC), os.path.join(os.getcwd(), 'out/{}'.format(cid), 'Data'))
        os.symlink('{}Documents'.format(PELE_SRC), os.path.join(os.getcwd(), 'out/{}'.format(cid), 'Documents'))

        # Minimization
        os.chdir("out/{}".format(cid))
        os.system(" %s %s > VACUUM_minimization.out" % (PELE_EXEC, VACUUM_CF))
        os.system(" %s %s > OBC_minimization.out" % (PELE_EXEC, OBC_CF))
        os.chdir("../..")

        # Calculate energetic difference
        difference = DiffEnergies.compute_energies("out/{}/VACUUM_minimization.out".format(cid), "out/{}/OBC_minimization.out".format(cid))[2]
        energies.append(tuple((cid, difference, exp_v)))
    except Exception:
        print('Skipping:', cid, '-', tag)
        continue


# Use pandas to save energies list of tuples to file
# CID    energetic diff    experimental value
df = pd.DataFrame(energies, columns = ['CID','Energetic Difference', 'Experimental value'])
df.to_csv('out/results.txt')



