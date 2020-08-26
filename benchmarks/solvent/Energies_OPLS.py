
import os
import pandas as pd
import DiffEnergies as DiffEnergies


PATH_TO_FREESOLV_DATABASE = '/home/lauramalo/repos/offpele-benchmarks/benchmarks/data/FreeSolv/FreeSolv0.52.txt'
PATH_OPLS = '/home/lauramalo/tests/FreeSolv/OPLS_out'

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
experimental_v = freesolv_db['experimental value'].to_list()

energies = list()

for cid, exp_v in zip(compound_ids, experimental_v):
    try:
        difference = DiffEnergies.compute_energies(os.path.join(PATH_OPLS,cid,'VACUUM_minimization.out'), os.path.join(PATH_OPLS,cid,'OBC_minimization.out'))[2]
        energies.append(tuple((cid, difference, exp_v)))

    except Exception:
        print('Skipping:', cid)
        continue

df = pd.DataFrame(energies, columns = ['CID','Energetic Difference', 'Experimental value'])
df.to_csv('OPLS_out/results.txt')


