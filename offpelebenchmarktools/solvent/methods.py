import os
from offpele.topology import Molecule, RotamerLibrary
from offpele.template import Impact
from offpele.solvent import OBC2
from offpele.main import handle_output_paths
from time import sleep
import shutil
from offpelebenchmarktools.solvent.energyhandler import compute_energies
from offpelebenchmarktools.utils.pele import PELEMinimization
from tqdm import tqdm
from functools import partial
from multiprocessing import Pool


def parallel_OFF_run(output_path, solvent, off_forcefield,
                     charges_method, pele_exec, pele_src, pele_license,
                     entry):
    cid, tag, exp_v = entry

    try:
        molecule = Molecule(smiles=tag, name=cid, tag='LIG')
        molecule.parameterize(off_forcefield,
                              charges_method=charges_method)

        # Minimization
        pele_vacuum_min = PELEMinimization(
            pele_exec, pele_src, pele_license,
            solvent_type='VACUUM',
            output_path=output_path)
        pele_vacuum_out = pele_vacuum_min.run(molecule,
                                              output_file='vacuum_out.txt')

        pele_obc_min = PELEMinimization(
            pele_exec, pele_src, pele_license,
            solvent_type=solvent,
            output_path=output_path)
        pele_obc_out = pele_obc_min.run(molecule,
                                        output_file='solvent_out.txt')

        # Calculate energetic difference
        difference = compute_energies(pele_vacuum_out, pele_obc_out)[2]

        return tuple(cid, difference, exp_v)

    except Exception as e:
        print('Exception found with compound {}: '.format(cid)
              + str(e))


def method_OFF(output_path, compound_ids, smiles_tags,
               experimental_v, solvent, off_forcefield,
               charges_method, pele_exec, pele_src, pele_license,
               n_proc=1):
    energies = list()

    parallel_runner = partial(parallel_OFF_run, output_path,
                              solvent, off_forcefield,
                              charges_method, pele_exec,
                              pele_src, pele_license)

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


def method_OPLS(output_path, compound_ids, smiles_tags, experimental_v,
                solvent, PELE_EXEC, SCHRODINGER_SRC, PLOPROTTEMP_SRC,
                pele_src, IMPACT_TEMPLATE_PATH, ROTAMER_LIBRARY_PATH, VACUUM_CF, OBC_CF, VDGBNP_CF):
    energies = list()
    diff_list = list()
    exp_list = list()
    for cid, tag, exp_v in zip(compound_ids, smiles_tags, experimental_v):
        try:
            # Create dir
            os.makedirs(os.path.join(output_path, cid), exist_ok=True)

            # Create representation of a particular molecule
            molecule = Molecule(smiles=tag, tag='UNL')

            # Write molecule's PDB
            molecule.to_pdb_file(os.path.join(output_path, cid, 'ligand.pdb'))

            # Generate MAE file and run PlopRotTemp
            os.chdir(os.path.join(output_path, cid))
            os.system('{}utilities/prepwizard '.format(SCHRODINGER_SRC)
                      + 'ligand.pdb ligand.mae -noepik '
                      + '-noprotassign -noccd -noimpref')
            for _ in range(0, 10):
                if os.path.isfile('ligand.mae'):
                    break
                else:
                    sleep(1)
            else:
                print('Time out')
                raise RuntimeError

            os.system('{}utilities/python {} ligand.mae'.format(SCHRODINGER_SRC,
                                                                PLOPROTTEMP_SRC))
            # Fix atom names from PlopRotTemp
            with open('ligand.pdb', 'r') as f:
                new_pdb_block = ''
                for line in f:
                    if line.startswith('HETATM'):
                        name = line[12:16]
                        if name[3] == ' ' and name[0] != ' ':
                            name = ' ' + name[0:3]

                        new_pdb_block += line[:12] + name + line[16:]
                    else:
                        new_pdb_block += line

            # Link folders
            try:
                os.symlink('{}Data'.format(pele_src), os.path.join(os.getcwd(), 'Data'))
                os.symlink('{}Documents'.format(pele_src), os.path.join(os.getcwd(), 'Documents'))
            except FileExistsError:
                pass

            # Generate OBC parameters
            os.makedirs('DataLocal/OBC/', exist_ok=True)
            os.system('python2 {}scripts/solventOBCParamsGenerator.py '.format(pele_src)
                      + 'unlz >> /dev/null')
            shutil.copy('{}Data/OBC/solventParamsHCTOBC.txt'.format(pele_src), 'DataLocal/OBC')
            os.system('cat unlz_OBCParams.txt >> DataLocal/OBC/solventParamsHCTOBC.txt')

            with open('ligand.pdb', 'w') as f:
                f.write(new_pdb_block)

            os.makedirs(IMPACT_TEMPLATE_PATH, exist_ok=True)
            os.makedirs(ROTAMER_LIBRARY_PATH, exist_ok=True)
            shutil.copy('unlz', IMPACT_TEMPLATE_PATH)
            shutil.copy('UNL.rot.assign', ROTAMER_LIBRARY_PATH)
            os.chdir('../..')

            # Generate solvent parameters, minimization and compute energetic difference
            if solvent == 'OBC':

                # Minimization
                os.chdir(os.path.join(output_path, cid))
                os.system(" %s %s > VACUUM_minimization.out" % (PELE_EXEC, VACUUM_CF))
                os.system(" %s %s > OBC_minimization.out" % (PELE_EXEC, OBC_CF))
                os.chdir("../..")

                # Calculate energetic difference
                difference = compute_energies(os.path.join(output_path, cid, 'VACUUM_minimization.out'), os.path.join(output_path, cid, 'OBC_minimization.out'))[2]
                energies.append(tuple((cid, difference, exp_v)))

            if solvent == 'VDGBNP':

                print('Start minimization')
                print(os.getcwd())
                # Minimization
                os.chdir(os.path.join(output_path, cid))
                os.system(" %s %s > VACUUM_minimization.out" % (PELE_EXEC, VACUUM_CF))
                os.system(" %s %s > VDGBNP_minimization.out" % (PELE_EXEC, VDGBNP_CF))
                os.chdir("../..")
                print('Finish minimization')

                # Calculate energetic difference
                difference = compute_energies(os.path.join(output_path, cid, 'VACUUM_minimization.out'), os.path.join(output_path, cid, 'VDGBNP_minimization.out'))[2]
                energies.append(tuple((cid, difference, exp_v)))
                diff_list.append(difference)
                exp_list.append(exp_v)
        except Exception:
            print('Skipping:', cid, '-', tag)
            continue
    return energies, diff_list, exp_list


def method_OFFOPLS(output_path, compound_ids, smiles_tags, experimental_v,
                   solvent, non_bonds, bonds_angles, off_forcefield,
                   charges_method, pele_exec, pele_src, pele_license):
    energies = list()
    diff_list = list()
    exp_list = list()
    for cid, tag, exp_v in zip(compound_ids, smiles_tags, experimental_v):
        try:
            # Create dir
            os.makedirs(os.path.join(output_path, cid), exist_ok=True)

            # Create representation of a particular molecule
            molecule = Molecule(smiles=tag, tag='UNL')

            # Handle paths
            rotamer_library_output_path, impact_output_path, \
                solvent_output_path = handle_output_paths(
                    molecule=molecule,
                    output=os.path.join(output_path, cid),
                    as_datalocal=True)

            # Generate its rotamer library
            rotamer_library = RotamerLibrary(molecule)
            rotamer_library.to_file(rotamer_library_output_path)

            # Generate its parameters and template file
            molecule.parameterize(OFF_FORCEFIELD,
                                  charges_method=charges_method)

            if non_bonds:
                molecule.add_OPLS_nonbonding_params()
            if bonds_angles:
                molecule.add_OPLS_bonds_and_angles()

            impact = Impact(molecule)
            impact.write(impact_output_path)

            # Generate its solvent parameters
            if solvent == 'VDGBNP':
                solvent = OBC2(molecule)
                solvent.to_json_file(solvent_output_path)

            # Write molecule's PDB
            molecule.to_pdb_file(os.path.join(output_path, cid,
                                              'ligand.pdb'))

            os.chdir(os.path.join(output_path, cid))
            os.makedirs(IMPACT_TEMPLATE_PATH, exist_ok=True)
            os.makedirs(ROTAMER_LIBRARY_PATH, exist_ok=True)
            shutil.copy('DataLocal/Templates/OFF/Parsley/HeteroAtoms/unlz',
                        IMPACT_TEMPLATE_PATH)
            os.chdir('../..')

            # Link folders
            try:
                os.symlink('{}Data'.format(pele_src),
                           os.path.join(output_path, cid, 'Data'))
                os.symlink('{}Documents'.format(pele_src),
                           os.path.join(output_path, cid, 'Documents'))
            except FileExistsError:
                pass

            # Minimization
            os.chdir(os.path.join(output_path, cid))
            os.system(" %s %s > VACUUM_minimization.out" % (pele_exec,
                                                            VACUUM_CF))
            os.system(" %s %s > OBC_minimization.out" % (pele_exec,
                                                         OBC_CF))
            os.chdir("../..")

            # Calculate energetic difference
            difference = compute_energies(os.path.join(
                output_path, cid, 'VACUUM_minimization.out'),
                os.path.join(output_path, cid, 'OBC_minimization.out'))[2]
            energies.append(tuple((cid, difference, exp_v)))
            diff_list.append(difference)
            exp_list.append(exp_v)

        except Exception:
            continue
    return energies, diff_list, exp_list
