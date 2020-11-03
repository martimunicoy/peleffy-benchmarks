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

    def run(self, data, output_path='output'):
        """
        It runs a bunch of minimizations.

        Parameters
        ----------
        data : dict[dict]
            Parsed data containing all the records and their attributes
            from the retrieved dataset. All of them will be minimized
        output_path : str
            The output path where results will be saved
        """
        from multiprocessing import Pool
        import json
        from tqdm import tqdm

        self._output_path = output_path

        index_to_name = dict()

        with Pool(self.n_proc) as p:
            list(tqdm(p.imap(self._parallel_minimizer,
                             enumerate(data.items())),
                      total=len(data.items())))

        for index, name in enumerate(data.keys()):
            index_to_name[index] = name

        json_output = json.dumps(index_to_name)

        with open(os.path.join(self._output_path,
                               'index_to_name.json'), "w") as f:
            f.write(json_output)

    def _parallel_minimizer(self, iteration_data):
        """
        It runs a minimization in parallel.

        Parameters
        ----------
        iteration_data : tuple[int, tuple[str, list[str]]]
            It contains the data for a certain minimization iteration
        """
        minimizer = Minimizer(PELE_exec=self._PELE_exec,
                              PELE_src=self._PELE_src)

        index, (name, attributes) = iteration_data

        smiles = attributes['canonical_isomeric_explicit_hydrogen_smiles']
        try:
            minimizer.minimize(smiles, str(index),
                               output_path=self._output_path)
        except Exception as e:
            print('Exception found for molecule {} {}'.format(name, smiles)
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

    def __init__(self, dataset, out_folder):
        """
        It initializes a MinimizationBenchmark object.

        Parameters:
        ----------
        dataset: name of the collection you want to extract from QCPortal
        out_folder: str
                    name for the output folder

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

        self.dataset = dataset
        self.out_folder = out_folder

    def _get_molecule_minimized(self, index):
        """
        It minimized the molecule using PELE's minimization.
        """

        from peleffy.topology import Molecule
        from peleffybenchmarktools.minimization import PELEMinimization

        # Load and parameterize the molecule
        try:
            mol = Molecule(os.path.join(os.getcwd(),
                                        self.out_folder,
                                        'QM/' '{}.pdb'.format(index + 1)))
            mol.parameterize('openff_unconstrained-1.2.1.offxml',
                             charge_method='gasteiger')

            # Runs a PELE Minimization
            pele_minimization = PELEMinimization(
                PELE_exec='/home/municoy/builds/PELE/PELE-repo_serial/PELE-1.6',
                PELE_src='/home/municoy/repos/PELE-repo/',
                PELE_license='/home/municoy/builds/PELE/license')
            pele_minimization.run(mol)

        except Exception:
            print('Skipping minimization for molecule {}'.format(index))

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

    def _filter_structures(self, smiles_tag):
        """
        It filter structures to keep only those that employ one of the
        OpenFF dihedrals with a non null phase constant.

        Output
        ------
        - It the structure employs one of the OpenFF dihedrals with a
          non null phase constant: return True.
        - Else: return False
        """

        # Filter out all non interesting molecules
        from peleffy.topology import Molecule
        from simtk import unit

        var = False
        try:
            mol = Molecule(smiles=smiles_tag)
            mol.parameterize('openff_unconstrained-1.2.1.offxml',
                             charge_method='gasteiger')
            for p in mol.propers:
                if p.phase not in (unit.Quantity(0, unit.degree),
                                   unit.Quantity(180, unit.degree)):
                    var = True
        except Exception:
            var = False
        return var

    def _get_dataset_structures(self, filter_angles):
        """
            It gets the Dataset from the QCPortal and saves in a folder the optimized molecules as PDB.
        """
        from peleffybenchmarktools.utils import QCPortal
        import qcportal as ptl

        # Get the optimization dataset
        client = ptl.FractalClient()
        ds = client.get_collection('OptimizationDataset', self.dataset)
        nmols = len(list(ds.data.records.keys()))

        # Initializes the QCPortal class
        qc_portal = QCPortal(n_proc=4)

        # Handles output paths
        os.mkdir(os.path.join(os.getcwd(), self.out_folder))
        set_folder = os.path.join(os.getcwd(), self.out_folder, 'QM')
        os.mkdir(set_folder)

        # Build the optimized molecule and save it as a PDB
        for index in range(nmols):
            if (index % 50) == 0:
                print('{}/{} optimized molecules built'.format(index + 1, nmols))
            entry = ds.get_entry(ds.df.index[index])
            if filter_angles:
                if self._filter_structures(smiles_tag=entry.attributes.get('canonical_smiles')):
                    qc_portal._parallel_struct_getter(ds, set_folder, index)
            else:
                qc_portal._parallel_struct_getter(ds, set_folder, index)

    def run(self, filter_structures=False):
        """
            It generates the folders for the optimized with QM structures to PDB files and these PDBs minnimized with PELE.
        """
        import shutil
        import glob
        import re

        # Obtain the dataset
        self._get_dataset_structures(filter_angles=filter_structures)

        # Gets all the PDB files from the set
        pdb_files = glob.glob(os.path.join(os.path.join(self.out_folder, 'QM'), "*pdb"))

        # Loads the molecules from the Dataset and runs a PELEMinimization
        for pdb_file in pdb_files:
            _, index = os.path.split(pdb_file)
            index = re.sub('.pdb', '', index)
            self._get_molecule_minimized(index=int(index))

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
