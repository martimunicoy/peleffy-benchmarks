"""
It prepares the system files and runs a PELE minimization.
"""
import os
LOCAL_DIR = os.path.dirname(os.path.abspath(__file__))


class MultiMinimizer(object):
    """
    It handles multiple calls to the Minimizer methods.
    """

    def __init__(self, PELE_exec, PELE_src, n_proc=1):
        """
        It initializes a MultiMinimizer object.

        Parameters
        ----------
        PELE_exec : str
            Path to the PELE executable
        PELE_src : str
            Path to PELE source folder
        n_proc : int
            The number of processors to employ to gather and parse data
        """
        # Supress INFO messages from peleffy
        from peleffy.utils import Logger
        log = Logger()
        log.set_level('WARNING')

        self._PELE_exec = PELE_exec
        self._PELE_src = PELE_src
        self._output_path = None
        self.n_proc = n_proc

    def run(self, pdb_paths, output_path='output'):
        """
        It runs a bunch of minimizations.

        Parameters
        ----------
        pdb_paths : list[str]
            List of PDB paths that contain the molecules to minimize
        output_path : str
            The output path where results will be saved
        """
        from multiprocessing import Pool
        from tqdm import tqdm

        self._output_path = output_path

        with Pool(self.n_proc) as p:
            list(tqdm(p.imap(self._parallel_minimizer,
                             enumerate(pdb_paths)),
                      total=len(pdb_paths)))

    def _parallel_minimizer(self, mol_ids):
        """
        It runs a minimization in parallel.

        Parameters
        ----------
        mol_ids : tuple([int, str])
            The identifiers ot he molecule to minimize
        """
        mol_id, mol_path = mol_ids
        minimizer = Minimizer(PELE_exec=self._PELE_exec,
                              PELE_src=self._PELE_src)

        try:
            minimizer.minimize(mol_path, mol_id,
                               output_path=self._output_path)
        except Exception as e:
            print('Exception found for molecule '
                  + '{} in {}'.format(mol_id, mol_path)
                  + ': ' + str(e))


class Minimizer(object):
    """
    It contains all the tools to minimize a molecule with PELE and the
    OpenForceField toolkit for PELE.
    """

    CONTROL_FILES = {
        'vacuum': os.path.join(LOCAL_DIR, 'data/VACUUM_minimization.conf'),
        'OBC': os.path.join(LOCAL_DIR, 'data/OBC_minimization.conf')}

    def __init__(self, PELE_exec, PELE_src):
        """
        It initializes a Minimizer object.

        Parameters
        ----------
        PELE_exec : str
            Path to the PELE executable
        PELE_src : str
            Path to PELE source folder
        """
        self._PELE_exec = PELE_exec
        self._PELE_src = PELE_src
        self._output_path = None

    def minimize(self, smiles, mol_id, solvent='vacuum',
                 forcefield='openff_unconstrained-1.2.0.offxml',
                 charge_method='am1bcc', output_path=None):
        """
        Given an input PDB file, it runs a minimization with PELE.

        Parameters
        ----------
        smiles : str
            The smiles tag representing the molecule to minimize
        mol_id : str
            Unique id to identify the molecule to minimize
        solvent : str
            The solvent name. One of ['vacuum', 'OBC']. Default is 'vacuum'
        forcefield : str
            The Open Force Field force field to generate the parameters
            with
        charge_method : str
            The charge method to calculate the partial charges with
        output_path : str
            The output path where results will be saved
        """

        if solvent not in self.CONTROL_FILES:
            raise ValueError('Invalid solvent:', solvent,
                             'It must be one of', self.CONTROL_FILES.keys())

        output_path = self._get_output_path(output_path, mol_id)

        self._create_directory(output_path)

        self._link_folders(output_path)

        self._generate_parameters(smiles, mol_id, output_path,
                                  forcefield=forcefield,
                                  charge_method=charge_method)

        self._run_PELE_minimization(solvent, output_path)

    def _get_output_path(self, output_path, mol_id):
        """
        It sets the output path.

        Parameters
        ----------
        output_path : str
            The output path where results will be saved
        mol_id : str
            Unique id to identify the molecule to minimize

        Returns
        -------
        output_path : str
            The output path that has been assigned
        """
        import os

        if output_path is None:
            output_path = os.path.join('output', str(mol_id))
        else:
            output_path = os.path.join(output_path, str(mol_id))

        return output_path

    def _create_directory(self, output_path):
        """
        It creates an output directory where all the results will be saved
        and copies the ligand's pdb inside.

        Parameters
        ----------
        output_path : str
            The output path where results will be saved
        """
        import os

        # It makes the output directory
        os.makedirs(output_path, exist_ok=True)

    def _link_folders(self, output_path):
        """
        It links the necessary folders to the output folder.

        Parameters
        ----------
        output_path : str
            The output path where results will be saved
        """
        import os

        # Link to Data
        link_path = os.path.join(os.getcwd(), output_path, 'Data')
        if os.path.isdir(link_path):
            os.remove(link_path)
        os.symlink(os.path.join(self._PELE_src, 'Data'),
                   link_path)

        # Link to Documents
        link_path = os.path.join(os.getcwd(), output_path, 'Documents')
        if os.path.isdir(link_path):
            os.remove(link_path)
        os.symlink(os.path.join(self._PELE_src, 'Documents'),
                   link_path)

    def _generate_parameters(self, smiles, mol_id, output_path,
                             forcefield='openff_unconstrained-1.2.0.offxml',
                             charge_method='am1bcc'):
        """
        It generates the parameters of the molecule (from the input_file)
        as DataLocal in the output folder.

        Parameters
        ----------
        smiles : str
            The smiles tag representing the molecule to minimize
        mol_id : str
            Unique id to identify the molecule to minimize
        output_path : str
            The output path where parameters will be saved
        forcefield : str
            The Open Force Field force field to generate the parameters
            with
        charge_method : str
            The charge method to calculate the partial charges with
        """
        import peleffy
        from peleffy.topology import Molecule
        from peleffy.template import Impact
        from peleffy.solvent import OBC2
        from peleffy.main import handle_output_paths
        import os

        # Create representation of a particular molecule
        molecule = Molecule(smiles=smiles, name=mol_id, tag='UNL')

        # Save molecule to PDB file
        molecule.to_pdb_file(os.path.join(output_path, 'ligand.pdb'))

        # Saving paths
        rotamer_library_output_path, impact_output_path, \
            solvent_output_path = handle_output_paths(molecule=molecule,
                                                      output=output_path,
                                                      as_datalocal=True)

        # Generate its rotamer library
        rotamer_library = peleffy.topology.RotamerLibrary(molecule)
        rotamer_library.to_file(rotamer_library_output_path)

        # Generate its parameters and template file
        molecule.parameterize(forcefield, charge_method=charge_method)
        impact = Impact(molecule)
        impact.write(impact_output_path)

        # Generate its solvent parameters
        solvent = OBC2(molecule)
        solvent.to_json_file(solvent_output_path)

    def _run_PELE_minimization(self, solvent, output_path):
        """
        It runs a PELE minimization.

        Parameters
        ----------
        solvent : str
            The solvent name. One of ['vacuum', 'OBC']. Default is 'vacuum'
        output_path : str
            The output path where parameters will be saved
        """
        import os

        # Minimization
        previous_dir = os.getcwd()
        os.chdir(os.path.join(os.getcwd(), output_path))
        os.system("{} {} > {}_minimization.out".format(
            self._PELE_exec, self.CONTROL_FILES[solvent], solvent))
        os.chdir(previous_dir)


class MinimizationBenchmark(object):
    """
    It defines a MinimizationBenchmark object.
    """

    def __init__(self, dataset_name, out_folder, n_proc=1):
        """
        It initializes a MinimizationBenchmark object.

        Parameters:
        ----------
        dataset_name: str
            The name of the collection you want to extract from QCPortal
        out_folder: str
            The name for the output folder
        n_proc : int
            The number of processors to employ

        Examples:
        ----------

        >>> benchmark = MinimizationBenchmark(
                dataset = 'Kinase Inhibitors: WBO Distributions',
                out_folder = 'SET1')
        >>> benchmark.run()

        """
        # Supress INFO messages from peleffy
        from peleffy.utils import Logger
        log = Logger()
        log.set_level('WARNING')

        self.dataset_name = dataset_name
        self.out_folder = out_folder
        self.n_proc = n_proc

    def _get_molecule_minimized(self, pdb_path):
        """
        It minimized the molecule using PELE's minimization.

        Parameters
        ----------
        pdb_path : str
            The path to the PDB of the molecule to minimize
        """
        from peleffy.topology import Molecule
        from peleffybenchmarktools.utils.pele import PELEMinimization

        # Load and parameterize the molecule
        try:
            mol = Molecule(pdb_path)
            mol.parameterize('openff_unconstrained-1.2.1.offxml',
                             charge_method='gasteiger')

            # Runs a PELE Minimization
            pele_minimization = PELEMinimization(
                PELE_exec='/home/municoy/builds/PELE/PELE-repo_serial/PELE-1.6',
                PELE_src='/home/municoy/repos/PELE-repo/',
                PELE_license='/home/municoy/builds/PELE/license')
            pele_minimization.run(mol)

        except Exception as e:
            print('Skipping minimization for molecule {}: '.format(pdb_path)
                  + str(e))

    def _parse_output_file(self, file):
        """
        Parsing the date from an output file after a PELEJob.
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

    def _has_nonstandard_dihedral(self, smiles_tag):
        """
        It parameterizes the molecule belonging to the smiles tag supplied
        to determine whether it contains any non-standard dihedral or not.
        A non-standard dihedral is a dihedral with a phase constant
        different from 0 and 180 degrees.

        Parameters
        ----------
        smiles_tag : str
            The smiles tag to construct the molecule with

        Returns
        -------
        answer : bool
            True if smiles' Molecule has a non-standard dihedral, False
            otherwise
        """
        from peleffy.topology import Molecule
        from simtk import unit
        from contextlib import suppress

        with suppress(Exception):
            mol = Molecule(smiles=smiles_tag)
            mol.parameterize('openff_unconstrained-1.2.1.offxml',
                             charge_method='gasteiger')
            for p in mol.propers:
                if p.phase not in (unit.Quantity(0, unit.degree),
                                   unit.Quantity(180, unit.degree)):
                    return True
        return False

    def _get_dataset_structures(self, filter_dihedrals):
        """
        It gets the Dataset from the QCPortal and saves in a folder
        the optimized molecules as PDB.

        Parameters
        ----------
        filter_dihedrals : bool
            Whether to filter entries to keep only non-standard dihedrals
            or use them all
        """
        from peleffybenchmarktools.utils import QCPortal
        import qcportal as ptl
        from functools import partial
        from multiprocessing import Pool
        from tqdm import tqdm

        # Get the optimization dataset
        client = ptl.FractalClient()
        ds = client.get_collection('OptimizationDataset', self.dataset_name)
        nmols = len(list(ds.data.records.keys()))

        # Initializes the QCPortal class
        qc_portal = QCPortal()

        # Handles output paths
        set_folder = os.path.join(os.getcwd(), self.out_folder, 'QM')
        os.makedirs(set_folder, exist_ok=True)

        filtered_indexes = list()
        if filter_dihedrals:
            print(' - Filtering non-standard dihedrals')
            smiles_tags = \
                [ds.get_entry(ds.df.index[index]).attributes.get(
                    'canonical_smiles') for index in range(nmols)]
            with Pool(self.n_proc) as p:
                contains_nonstandard_dihedral = \
                    list(tqdm(p.imap(self._has_nonstandard_dihedral,
                                     smiles_tags),
                              total=len(smiles_tags)))

            for idx, nonstandard_dihedral_inside in \
                    enumerate(contains_nonstandard_dihedral):
                if nonstandard_dihedral_inside:
                    filtered_indexes.append(idx)

        else:
            filtered_indexes = [idx for idx in range(nmols)]

        parallel_function = partial(qc_portal._parallel_struct_getter,
                                    ds, set_folder)

        # Build optimized molecules and save them as a PDB
        print(' - Generating optimized molecules')
        with Pool(self.n_proc) as p:
            list(tqdm(p.imap(parallel_function,
                             filtered_indexes),
                      total=len(filtered_indexes)))

    def run(self, filter_dihedrals=False):
        """
        It generates the folders for the optimized with QM structures
        to PDB files and these PDBs minnimized with PELE.

        Parameters
        ----------
        filter_dihedrals : bool
            Whether to filter entries to keep only non-standard dihedrals
            or use them all
        """
        import shutil
        import glob
        from multiprocessing import Pool
        from tqdm import tqdm

        # Delete output folder, if it already exists

        # Obtain the dataset
        self._get_dataset_structures(filter_dihedrals=filter_dihedrals)

        # Gets all the PDB files from the set
        pdb_files = glob.glob(
            os.path.join(os.path.join(self.out_folder, 'QM'), "*pdb"))

        # Loads the molecules from the Dataset and runs a PELEMinimization
        print(' - Minimizing molecules')
        with Pool(self.n_proc) as p:
            list(tqdm(p.imap(self._get_molecule_minimized,
                             pdb_files),
                      total=len(pdb_files)))

        # Moves the output folder(created by the PELEMinimization) to the desired output folder
        shutil.move(os.path.join(os.getcwd(), 'output'), os.path.join(os.getcwd(), self.out_folder, 'PELE'))

    def compute_RMSD(self):
        """
            For a collection of structures,
            it saves a CSV file with a dictionary of the RMSD comparison between PELE and QM minimized structures.
            It generats an histogram of the computed RMSD values.
        """
        import mdtraj as md
        import glob
        import re
        import pandas as pd
        import matplotlib.pyplot as plt

        # Gets all the structures from the output folder
        pdb_files = glob.glob(os.path.join(os.path.join(self.out_folder, 'QM'), "*pdb"))

        # Computes the RMSD between PELE and QM minimized structures
        d = {}
        rmsd_results = []
        for pdb_file in pdb_files:
            _, index = os.path.split(pdb_file)
            index = re.sub('.pdb', '', index)
            try:
                molQM = md.load_pdb(os.path.join(self.out_folder, 'QM/{}.pdb'.format(index)))
                molPELE = md.load_pdb(os.path.join(self.out_folder, 'PELE/{}/minimized.pdb'.format(index)))
                rsmd = md.rmsd(molQM, molPELE)[0]
                d.update({index: rsmd})
                rmsd_results.append(rsmd)
            except Exception:
                pass

        # Writes out a CSV file with the dictionary of the RMSD results.
        df = pd.DataFrame(d.items(), columns=['Ligand ID', 'RMSD'])
        df.to_csv(os.path.join(self.out_folder, 'rmsd.csv'))

        # Plots an histogram of the computed RMSD values
        plt.figure(figsize=(10, 7))
        plt.hist(rmsd_results, bins=10, rwidth=0.8, align='mid', range=(0, 1), color='gray')
        plt.xlabel('RMSD')
        plt.ylabel('Frequency')
        plt.title('Structural Histogram')
        plt.savefig(os.path.join(self.out_folder, 'rmsd.png'))

    def energetic_difference(self):
        """
            For a ollection of structures, gets the energy before and after the minimization and computes its difference.
            It generates a file with the energies and an energetic histogram.
        """
        import glob
        import re
        import pandas as pd
        import maplotlib.pyplot as plt

        # Gets all the structures from the output folder
        pdb_files = glob.glob(os.path.join(os.path.join(self.out_folder, 'QM'), "*pdb"))

        # Computes the energetic difference between PELE and QM minimized structures
        d = {}
        energetic_diff = []
        for pdb_file in pdb_files:
            _, index = os.path.split(pdb_file)
            index = re.sub('.pdb', '', index)
            try:
                file_path = os.path.join(self.out_folder, 'PELE/{}/PELE_output.txt'.format(index))
                data = self._parse_output_file(file=file_path)
                e_ini = data[-1].get('ENERGY VACUUM + CONSTRAINTS')
                e_fin = data[-2].get('ENERGY VACUUM + CONSTRAINTS')
                e_diff = float(e_fin) - float(e_ini)
                energetic_diff.append(e_diff)
                d.update({index: tuple((e_ini, e_fin, e_diff))})
            except Exception:
                pass

        # Writes out a CSV file with the dictionary of the energies results.
        df = pd.DataFrame(d.items(), columns=['Ligand ID', 'Energies'])
        df.to_csv(os.path.join(self.out_folder, 'energies.csv'))

        # Plots an energetic histogram
        plt.figure(figsize=(10, 7))
        plt.hist(energetic_diff, bins=10, rwidth=0.8, align='mid', color='green')
        plt.xlabel('Energetic difference (kcal/mol)')
        plt.ylabel('Frequency')
        plt.title('Energetic Histogram')
        plt.savefig(os.path.join(self.out_folder, 'energies.png'))
