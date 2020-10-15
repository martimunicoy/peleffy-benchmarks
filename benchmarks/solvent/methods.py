import os
import argparse 
import pandas as pd
from offpele.topology import Molecule, RotamerLibrary
from offpele.template import Impact
from offpele.solvent import OBC2
from offpele.main import handle_output_paths
import DiffEnergies as DiffEnergies
from time import sleep
import shutil


def method_OFF(OUTPUT_PATH, compound_ids, smiles_tags, experimental_v, OFF_FORCEFIELD, CHARGES_METHOD, PELE_EXEC, VACUUM_CF, OBC_CF,PELE_SRC):
    energies = list()
    diff_list = list()
    exp_list = list()
    for cid, tag, exp_v in zip(compound_ids, smiles_tags, experimental_v):
        try: 
            # Create dir
            os.makedirs(os.path.join(OUTPUT_PATH, cid), exist_ok=True)

            # Create representation of a particular molecule
            molecule = Molecule(smiles=tag, tag='UNL')

            # Handle paths
            rotamer_library_output_path, impact_output_path, solvent_output_path = \
                handle_output_paths(molecule=molecule, output=os.path.join(OUTPUT_PATH, cid), as_datalocal=True)

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
            molecule.to_pdb_file(os.path.join(OUTPUT_PATH, cid, 'ligand.pdb'))

            # Link folders
            try:
                os.symlink('{}Data'.format(PELE_SRC), os.path.join(os.getcwd(), OUTPUT_PATH, cid, 'Data'))
                os.symlink('{}Documents'.format(PELE_SRC), os.path.join(os.getcwd(), OUTPUT_PATH, cid, 'Documents'))
            except FileExistsError:
                pass

            # Minimization
            os.chdir(os.path.join(OUTPUT_PATH, cid))
            os.system(" %s %s > VACUUM_minimization.out" % (PELE_EXEC, VACUUM_CF))
            os.system(" %s %s > OBC_minimization.out" % (PELE_EXEC, OBC_CF))
            os.chdir("../..")

            # Calculate energetic difference
            difference = DiffEnergies.compute_energies(os.path.join(OUTPUT_PATH, cid,'VACUUM_minimization.out'), os.path.join(OUTPUT_PATH, cid,'OBC_minimization.out'))[2]
            energies.append(tuple((cid, difference, exp_v)))
            diff_list.append(difference)
            exp_list.append(exp_v)
        except Exception:
            continue
    return energies, diff_list, exp_list


def method_OPLS(OUTPUT_PATH, compound_ids, smiles_tags, experimental_v,solvent, PELE_EXEC, SCHRODINGER_SRC, PLOPROTTEMP_SRC, PELE_SRC, IMPACT_TEMPLATE_PATH, ROTAMER_LIBRARY_PATH, VACUUM_CF, OBC_CF,VDGBNP_CF):
    energies = list()
    diff_list = list()
    exp_list = list()
    for cid, tag, exp_v in zip(compound_ids, smiles_tags, experimental_v):
        try:
            # Create dir
            os.makedirs(os.path.join(OUTPUT_PATH, cid), exist_ok=True)

            # Create representation of a particular molecule
            molecule = Molecule(smiles=tag, tag='UNL')

            # Write molecule's PDB
            molecule.to_pdb_file(os.path.join(OUTPUT_PATH, cid, 'ligand.pdb'))

            # Generate MAE file and run PlopRotTemp
            os.chdir(os.path.join(OUTPUT_PATH, cid))
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
                os.symlink('{}Data'.format(PELE_SRC), os.path.join(os.getcwd(),  'Data'))
                os.symlink('{}Documents'.format(PELE_SRC), os.path.join(os.getcwd(), 'Documents'))
            except FileExistsError:
                pass
            
            # Generate OBC parameters
            os.makedirs('DataLocal/OBC/', exist_ok=True)
            os.system('python2 {}scripts/solventOBCParamsGenerator.py '.format(PELE_SRC)
                      + 'unlz >> /dev/null')
            shutil.copy('{}Data/OBC/solventParamsHCTOBC.txt'.format(PELE_SRC), 'DataLocal/OBC')
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
                os.chdir(os.path.join(OUTPUT_PATH, cid))
                os.system(" %s %s > VACUUM_minimization.out" % (PELE_EXEC, VACUUM_CF))
                os.system(" %s %s > OBC_minimization.out" % (PELE_EXEC, OBC_CF))
                os.chdir("../..")

                # Calculate energetic difference
                difference = DiffEnergies.compute_energies(os.path.join(OUTPUT_PATH,cid,'VACUUM_minimization.out'), os.path.join(OUTPUT_PATH,cid,'OBC_minimization.out'))[2]
                energies.append(tuple((cid, difference, exp_v)))
            
            if solvent == 'VDGBNP':

                print('Start minimization')
                print(os.getcwd())
                # Minimization
                os.chdir(os.path.join(OUTPUT_PATH, cid))
       	        os.system(" %s %s > VACUUM_minimization.out" % (PELE_EXEC, VACUUM_CF))
                os.system(" %s %s > VDGBNP_minimization.out" % (PELE_EXEC, VDGBNP_CF))
                os.chdir("../..")
                print('Finish minimization')

                # Calculate energetic difference
                difference = DiffEnergies.compute_energies(os.path.join(OUTPUT_PATH, cid, 'VACUUM_minimization.out'), os.path.join(OUTPUT_PATH, cid, 'VDGBNP_minimization.out'))[2]
                energies.append(tuple((cid, difference, exp_v)))
                diff_list.append(difference)
                exp_list.append(exp_v)
        except Exception :
            print('Skipping:', cid, '-', tag)
            continue
    return energies, diff_list, exp_list


def method_OFFOPLS(OUTPUT_PATH, compound_ids, smiles_tags, experimental_v, solvent, non_bonds, bonds_angles, OFF_FORCEFIELD, PELE_EXEC,
    CHARGES_METHOD,  PELE_SRC, VACUUM_CF, OBC_CF, IMPACT_TEMPLATE_PATH, ROTAMER_LIBRARY_PATH):
    energies = list()
    diff_list = list()
    exp_list = list()
    for cid, tag, exp_v in zip(compound_ids, smiles_tags, experimental_v):
        try: 
            # Create dir
            os.makedirs(os.path.join(OUTPUT_PATH, cid), exist_ok=True)

            # Create representation of a particular molecule
            molecule = Molecule(smiles=tag, tag='UNL')

            # Handle paths
            rotamer_library_output_path, impact_output_path, solvent_output_path = \
                handle_output_paths(molecule=molecule, output=os.path.join(OUTPUT_PATH, cid), as_datalocal=True)

            # Generate its rotamer library
            rotamer_library = RotamerLibrary(molecule)
            rotamer_library.to_file(rotamer_library_output_path)

            # Generate its parameters and template file
            molecule.parameterize(OFF_FORCEFIELD, charges_method=CHARGES_METHOD)

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
            molecule.to_pdb_file(os.path.join(OUTPUT_PATH, cid, 'ligand.pdb'))

            os.chdir(os.path.join(OUTPUT_PATH, cid))
            os.makedirs(IMPACT_TEMPLATE_PATH, exist_ok=True)
            os.makedirs(ROTAMER_LIBRARY_PATH, exist_ok=True)
            shutil.copy('DataLocal/Templates/OFF/Parsley/HeteroAtoms/unlz', IMPACT_TEMPLATE_PATH)
            os.chdir('../..')

            # Link folders
            try:
                os.symlink('{}Data'.format(PELE_SRC), os.path.join( OUTPUT_PATH, cid, 'Data'))
                os.symlink('{}Documents'.format(PELE_SRC), os.path.join( OUTPUT_PATH, cid, 'Documents'))
            except FileExistsError:
                pass

            # Minimization
            os.chdir(os.path.join(OUTPUT_PATH, cid))
            os.system(" %s %s > VACUUM_minimization.out" % (PELE_EXEC, VACUUM_CF))
            os.system(" %s %s > OBC_minimization.out" % (PELE_EXEC, OBC_CF))
            os.chdir("../..")

            # Calculate energetic difference
            difference = DiffEnergies.compute_energies(os.path.join(OUTPUT_PATH, cid,'VACUUM_minimization.out'), os.path.join(OUTPUT_PATH, cid,'OBC_minimization.out'))[2]
            energies.append(tuple((cid, difference, exp_v)))
            diff_list.append(difference)
            exp_list.append(exp_v)

        except Exception:
            continue
    return energies, diff_list, exp_list
