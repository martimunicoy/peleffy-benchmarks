"""
It runs two PELE minimizations, at vacuum and in solvent, and
outputs the energetic difference (hydration free energy).
"""


import os
import shutil
from peleffy.topology import Molecule
from peleffybenchmarktools.utils.pele import PELEMinimization
from tqdm import tqdm
from functools import partial
from multiprocessing import Pool


def parallel_run(output_path, solvent, charge_method,
                 pele_exec, pele_src, pele_license,
                 forcefield_name, forcefield, entry):
    """Parallel runner."""

    cid, tag, exp_v = entry

    try:
        molecule = Molecule(smiles=tag, name=cid, tag='LIG')
        if forcefield_name is not None:
            molecule.parameterize(forcefield_name,
                                  charge_method=charge_method)
        elif forcefield is not None:
            molecule.set_forcefield(forcefield)
            molecule.parameterize(charge_method=charge_method)

        if ((molecule.forcefield.type == 'OPLS2005')
            or (molecule.forcefield.type == 'OpenFF + OPLS2005'
                and molecule.forcefield._nonbonding == 'opls2005')):
            if solvent == 'OBC':
                # Generate OBC parameters for OPLS2005
                from peleffy.template import Impact

                os.makedirs(os.path.join(output_path, cid), exist_ok=True)

                impact = Impact(molecule)
                impact.write(os.path.join(output_path, cid, 'ligz'))

                os.makedirs(os.path.join(output_path, cid, 'DataLocal/OBC/'),
                            exist_ok=True)
                os.system('python2 {}scripts/solventOBCParamsGenerator.py '.format(pele_src)
                          + os.path.join(output_path, cid, 'ligz') + ' >> /dev/null')
                shutil.copy('{}Data/OBC/solventParamsHCTOBC.txt'.format(pele_src),
                            os.path.join(output_path, cid, 'DataLocal/OBC/'))
                os.system('cat ' + os.path.join(output_path, cid, 'ligz_OBCParams.txt >> '
                                                + os.path.join(output_path, cid,
                                                               'DataLocal/OBC/solventParamsHCTOBC.txt')))

                os.remove(os.path.join(output_path, cid, 'ligz'))
                os.remove(os.path.join(output_path, cid, 'ligz_OBCParams.txt'))

            forcefield_tag = 'OPLS2005'
        else:
            forcefield_tag = 'OpenForceField'

        # Minimization
        pele_vacuum_min = PELEMinimization(
            pele_exec, pele_src, pele_license,
            solvent_type='VACUUM',
            output_path=output_path,
            forcefield=forcefield_tag)
        pele_vacuum_out = pele_vacuum_min.run(molecule,
                                              output_file='vacuum_out.txt')

        pele_obc_min = PELEMinimization(
            pele_exec, pele_src, pele_license,
            solvent_type=solvent,
            output_path=output_path,
            forcefield=forcefield_tag)
        pele_obc_out = pele_obc_min.run(molecule,
                                        output_file='solvent_out.txt')

        # Calculate energetic difference
        difference = compute_energies(pele_vacuum_out, pele_obc_out)[2]

        return tuple((cid, difference, exp_v))

    except Exception as e:
        print('Exception found with compound {}: '.format(cid)
              + str(e))


def runner(output_path, compound_ids, smiles_tags,
           experimental_v, solvent, charge_method,
           pele_exec, pele_src, pele_license,
           n_proc=1, forcefield_name=None, forcefield=None):
    """Main runner."""
    if forcefield_name is None and forcefield is None:
        raise Exception('Either a forcefield_name or a forcefield object '
                        + 'is required')

    energies = list()

    parallel_runner = partial(parallel_run, output_path,
                              solvent, charge_method,
                              pele_exec, pele_src, pele_license,
                              forcefield_name, forcefield)

    with Pool(n_proc) as p:
        listed_data = list(tqdm(p.imap(parallel_runner,
                                       zip(compound_ids, smiles_tags,
                                           experimental_v)),
                                total=len(compound_ids)))

    for entry in listed_data:
        if entry is None:
            continue
        energies.append(entry)

    return energies


def parse_file(file):
    """
    Parsing the date from the input file.
    """
    import re

    with open(file, 'r') as f:
        data = list()
        group = dict()
        for key, value in re.findall(r'(.*):\s*([\dE+-.]+)', f.read()):
            if key in group:
                data.append(group)
                group = dict()
            group[key] = value
        data.append(group)
    return data


def get_energy(file):
    """
    Given a file, it does a data parsing of the file and creates a dictionary
    of the energies after minimization. Returns the total energy computed
    after the minimization.
    """
    data = parse_file(file)
    try:
        energies_minimized = data[-1]
        label = list(energies_minimized.keys())[-1]
    except IndexError:
        raise Exception('An invalid PELE output file has been found')
    return (energies_minimized.get(label))


def compute_energies(vacuum_file, OBC_file):
    """
    Given two files, it returns the vacuum energy, the OBC energy and the
    difference energy.
    """
    vacuum_energy = get_energy(vacuum_file)
    OBC_energy = get_energy(OBC_file)
    hydration_energy = float(OBC_energy) - float(vacuum_energy)
    return vacuum_energy, OBC_energy, hydration_energy
