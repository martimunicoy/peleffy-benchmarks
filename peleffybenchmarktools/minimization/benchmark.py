"""
This module contains classes and functions that handle minimization
benchmarks.
"""


class MinimizationBenchmark(object):
    """
    It defines a MinimizationBenchmark object.
    """

    def __init__(self, dataset_name, output_path, n_proc=1,
                 PELE_exec, PELE_src, PELE_license,
                 geometry_selection='optimized'):
        """
        It initializes a MinimizationBenchmark object.

        Parameters:
        ----------
        dataset_name: str
            The name of the collection you want to extract from QCPortal
        output_path: str
            The name for the output folder
        n_proc : int
            The number of processors to employ
        geometry_selection : str
            The geometry to feed to the molecule. One of
            ['initial', 'optimized']

        Examples:
        ----------

        >>> benchmark = MinimizationBenchmark(
                dataset = 'Kinase Inhibitors: WBO Distributions',
                output_path = 'SET1')
        >>> benchmark.run()

        """
        # Preload matplotlib to prevent Jupyter Notebook from hiding
        # the first plot
        import matplotlib.pyplot as plt
        plt.__name__

        # Supress INFO messages from peleffy
        from peleffy.utils import Logger
        log = Logger()
        log.set_level('WARNING')

        self.dataset_name = dataset_name
        self.output_path = output_path
        self.n_proc = n_proc
        self.geometry_selection = geometry_selection
        self.PELE_exec = PELE_exec
        self.PELE_src = PELE_src
        self.PELE_license = PELE_license

    def _get_molecule_minimized(self, pdb_path, output_path):
        """
        It minimized the molecule using PELE's minimization.

        Parameters
        ----------
        pdb_path : str
            The path to the PDB of the molecule to minimize
        output_path : str
            The path where to run the PELE simulation and save the output

        Returns
        -------
        output_file : str
            The path to the PELE minimized PDB file. It is None if PELE
            exits with a non-zero code
        """
        import os
        from peleffy.topology import Molecule
        from peleffybenchmarktools.utils.pele import PELEMinimization

        # Load and parameterize the molecule
        try:
            mol = Molecule(pdb_path)
            mol.parameterize('openff_unconstrained-1.2.1.offxml',
                             charge_method='gasteiger')

            # Runs a PELE Minimization
            pele_minimization = PELEMinimization(
                PELE_exec=self.PELE_exec,
                PELE_src=self.PELE_src,
                PELE_license=self.PELE_license,
                output_path=output_path)
            output_file = pele_minimization.run(mol,
                                                output_file='PELE_output.txt')

            # Return not the path to the PELE output file but the path to
            # the PELE minimized PDB file
            # Be careful, we are using a hardcoded name
            return os.path.join(os.path.dirname(output_file), 'minimized.pdb')

        except Exception as e:
            print('Skipping minimization for molecule {}: '.format(pdb_path)
                  + str(e))

            return None

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

    def _get_filtered_indexes(self, filter_dihedrals, ds):
        """
        It returns a list of booleans pointing out which structures can be
        filtered out. Those structures filtered will have only standard
        dihedrals.

        Parameters
        ----------
        filter_dihedrals : bool
            Whether to filter by dihedrals or not
        ds : a QCPortal.OptimizationDataset
            The dataset containing the entries to analyze

        Returns
        -------
        filtered_indexes : list[bool]
            The list of booleans. If False the corresponding structure
            must be filtered out
        """
        from multiprocessing import Pool
        from tqdm import tqdm

        filtered_indexes = list()
        nmols = len(list(ds.data.records.keys()))

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

        return filtered_indexes

    def retrieve_reference_pdbs(self, filter_dihedrals=False):
        """
        It gets the Dataset from the QCPortal and saves in a folder
        the optimized molecules as PDB.

        Parameters
        ----------
        filter_dihedrals : bool
            Whether to filter entries to keep only non-standard dihedrals
            or use them all

        Returns
        -------
        pdb_paths : list[str]
            The list containing the paths to generated PDB files
        """
        import os
        from peleffybenchmarktools.utils import QCPortal
        import qcportal as ptl
        from functools import partial
        from multiprocessing import Pool
        from tqdm import tqdm

        # Get the optimization dataset
        client = ptl.FractalClient()
        ds = client.get_collection('OptimizationDataset', self.dataset_name)

        # Initializes the QCPortal class
        qc_portal = QCPortal()

        # Handles output paths
        current_output = os.path.join(os.getcwd(), self.output_path, 'QM')
        os.makedirs(current_output, exist_ok=True)

        filtered_indexes = self._get_filtered_indexes(filter_dihedrals, ds)

        parallel_function = partial(qc_portal._parallel_struct_getter,
                                    ds, current_output,
                                    self.geometry_selection)

        # Build optimized molecules and save them as a PDB
        print(' - Generating optimized molecules')
        with Pool(self.n_proc) as p:
            pdb_paths = list(tqdm(p.imap(parallel_function,
                                         filtered_indexes),
                                  total=len(filtered_indexes)))

        # Filter out None pdb paths that belong to failing builds
        pdb_paths = [p for p in pdb_paths if p is not None]

        return pdb_paths

    def minimize_structures(self, pdb_paths):
        """
        It runs a PELE minimization to each of the PDB files from the
        the supplied list.

        Parameters
        ----------
        pdb_paths : list[str]
            The list containing the paths to generated PDB files

        Returns
        -------
        minimization_paths : list[str]
            The list containing the minimized PDB files.
        """
        import os
        from multiprocessing import Pool
        from tqdm import tqdm
        from functools import partial

        # Handles output paths
        current_output = os.path.join(os.getcwd(), self.output_path, 'PELE')
        os.makedirs(current_output, exist_ok=True)

        parallel_function = partial(self._get_molecule_minimized,
                                    current_output)

        # Loads the molecules from the Dataset and runs a PELEMinimization
        print(' - Minimizing molecules')
        with Pool(self.n_proc) as p:
            pdb_paths = list(tqdm(p.imap(parallel_function, pdb_paths),
                                  total=len(pdb_paths)))

        # Filter out None pdb paths that belong to failing builds
        pdb_paths = [p for p in pdb_paths if p is not None]

        return pdb_paths

    def _link_pdb_paths(self, paths_set1, paths_set2, labeling1, labeling2):
        """
        Given two sets of PDB paths, it links them based on the
        supplied labeling criteria.

        Parameters
        ----------
        paths_set1 : list[str]
            The first set of PDB paths
        paths_set2 : list[str]
            The second set of PDB paths
        labeling1 : str
            Labeling criteria for the first set. One of ['file', 'folder']
        labeling2 : str
            Labeling criteria for the first set. One of ['file', 'folder']

        Returns
        -------
        links : dict
            The dictionary keyed with PDB indexes and the corresponding
            pairing of PDB paths
        """
        import os

        if labeling1 not in ['file', 'folder']:
            raise NameError('Wrong selected labeling for set 1: '
                            + labeling1)

        if labeling2 not in ['file', 'folder']:
            raise NameError('Wrong selected labeling for set 2: '
                            + labeling2)

        set1 = set(paths_set1)
        set2 = set(paths_set2)
        links = dict()

        # Find link between a molecule from each set
        # Note that they are linked by either filename or folder
        while(len(set1) != 0):
            pdb_file1 = set1.pop()
            if labeling1 == 'file':
                name1 = os.path.basename(os.path.dirname(pdb_file1))
            else:
                name1 = os.path.basename(pdb_file1)
            for pdb_file2 in set2:
                if labeling2 == 'file':
                    name2 = os.path.basename(os.path.dirname(pdb_file2))
                else:
                    name2 = os.path.basename(pdb_file2)
                if name1 == name2:
                    links[name1] = (pdb_file1, pdb_file2)
                    set2.remove(pdb_file2)
                    break
            else:
                print('- Warning PDB from set 1: {} '.format(pdb_file1)
                      + 'not found in set 2')

        for pdb_file in set2:
            print('- Warning PDB from set 2: {} '.format(pdb_file)
                  + 'not found in set 1')

        return links

    def compute_RMSD(self, paths_set1, paths_set2,
                     labeling1='file', labeling2='file'):
        """
        For a collection of structures stored in two sets, it saves a
        CSV file with a dictionary of the RMSD comparison between
        PELE and QM minimized structures. It also generates an histogram
        of the computed RMSD values.

        Parameters
        ----------
        paths_set1 : list[str]
            The first set of PDB paths
        paths_set2 : list[str]
            The second set of PDB paths
        labeling1 : str
            Labeling criteria for the first set. One of ['file', 'folder']
        labeling2 : str
            Labeling criteria for the first set. One of ['file', 'folder']
        """
        import os
        import mdtraj as md
        import pandas as pd
        import matplotlib.pyplot as plt

        # Find links between the PDB paths from each set
        links = self._link_pdb_paths(paths_set1, paths_set2, labeling1,
                                     labeling2)

        # Computes the RMSD between PELE and QM minimized structures
        d = {}
        rmsd_results = []
        for idx, (pdb_file1, pdb_file2) in links.items():
            try:
                molQM = md.load_pdb(pdb_file1)
                molPELE = md.load_pdb(pdb_file2)
                rsmd = md.rmsd(molQM, molPELE)[0]
                d.update({idx: rsmd})
                rmsd_results.append(rsmd)
            except Exception as e:
                print('Skipping RMSD comparison between \''
                      + pdb_file1 + '\' \'' + pdb_file2 + '\': '
                      + str(e))

        # Writes out a CSV file with the dictionary of the RMSD results.
        df = pd.DataFrame(d.items(), columns=['Ligand ID', 'RMSD'])
        df.to_csv(os.path.join(self.out_folder, 'rmsd.csv'))

        # Plots an histogram of the computed RMSD values
        plt.figure(figsize=(7, 5))
        plt.hist(rmsd_results, bins=10, rwidth=0.8, align='mid',
                 range=(0, 0.5), color='gray')
        plt.xlabel('RMSD')
        plt.ylabel('Frequency')
        plt.title('Structural Histogram')
        plt.savefig(os.path.join(self.out_folder, 'rmsd.png'))

    def energetic_difference(self, minimized_pdb_paths):
        """
        For a collection of structures stored in two sets, it computes
        and displays the distribution of energetic differences in an
        histogram.

        Parameters
        ----------
        minimized_pdb_paths : list[str]
            The list of minimized PDB paths
        """
        import os
        import pandas as pd
        import matplotlib.pyplot as plt

        # Computes the energetic difference between PELE and QM minimized structures
        d = {}
        energetic_diff = []
        for pdb_file in minimized_pdb_paths:
            file_path = os.path.join(os.path.dirname(pdb_file),
                                     'PELE_output.txt')
            idx = os.path.basename(os.path.dirname(pdb_file))

            try:
                data = self._parse_output_file(file=file_path)
                e_ini = data[-1].get('ENERGY VACUUM + CONSTRAINTS')
                e_fin = data[-2].get('ENERGY VACUUM + CONSTRAINTS')
                e_diff = float(e_fin) - float(e_ini)
                energetic_diff.append(e_diff)
                d.update({idx: tuple((e_ini, e_fin, e_diff))})
            except Exception as e:
                print('Skipping energetic difference calculation for \''
                      + file_path + '\': '
                      + str(e))

        # Writes out a CSV file with the dictionary of the energies results.
        df = pd.DataFrame(d.items(), columns=['Ligand ID', 'Energies'])
        df.to_csv(os.path.join(self.out_folder, 'energies.csv'))

        # Plots an energetic histogram
        plt.figure(figsize=(7, 5))
        plt.hist(energetic_diff, bins=10, rwidth=0.8, align='mid', color='green')
        plt.xlabel('Energetic difference (kcal/mol)')
        plt.ylabel('Frequency')
        plt.title('Energetic Histogram')
        plt.savefig(os.path.join(self.out_folder, 'energies.png'))