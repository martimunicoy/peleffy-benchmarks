"""
This module contains classes and functions that handle minimization
benchmarks.
"""


class MinimizationBenchmark(object):
    """
    It defines a MinimizationBenchmark object.
    """

    def __init__(self, dataset_name, output_path,
                 PELE_exec, PELE_src, PELE_license,
                 geometry_selection='optimized',
                 n_proc=1,
                 forcefield='openff_unconstrained-1.3.0.offxml',
                 charge_method=None,
                 distort_bonds=False, distort_value_for_bonds=0.5,
                 distort_torsions=False, distort_value_for_torsions=20,
                 distort_dihedrals=False, distort_value_for_dihedrals=20,
                 random_distort=False,
                 bond_constants_factor=1.0,
                 torsion_constants_factor=1.0,
                 dihedral_constants_factor=1.0,
                 force_parameterization=False):
        """
        It initializes a MinimizationBenchmark object.

        Parameters:
        ----------
        dataset_name: str
            The name of the collection you want to extract from QCPortal
        output_path: str
            The name for the output folder
        geometry_selection : str
            The geometry to feed to the molecule. One of
            ['initial', 'optimized']
        n_proc : int
            The number of processors to employ. Default is 1
        PELE_exec : str
            Path to the PELE executable
        PELE_src : str
            Path to PELE source folder
        PELE_license : str
            Path to PELE license directory
        forcefield : str
            The force field to employ in the parameterization job. Defaut is
            'openff_unconstrained-1.3.0.offxml'
        charge_method : str
            The charge method to employ in the parameterization job. Default
            is None
        distort_bonds : bool
            Whether to distort structural bonds or not. Default is False
        distort_torsions : bool
            Whether to distort structural torsions or not. Default is False
        distort_dihedrals : bool
            Whether to distort structural dihedrals or not. Default is False
        distort_value_for_bonds : float
            The value to find distortion lenghts for bonds, in
            angstroms. Default is 0.5
        distort_value_for_torsions : float
            The value to find distortion angles for torsions, in
            degrees. Default is 20
        distort_value_for_dihedrals : float
            The value to find distortion angles for dihedrals, in
            degrees. Default is 20
        random_distort : bool
            Whether to apply a random distort or a fixed one. Default is
            False
        bond_constants_factor : float
            The factor to customize bond force constants. Default is 1.0
        torsion_constants_factor : float
            The factor to customize torsion force constants. Default is 1.0
        dihedral_constants_factor : float
            The factor to customize dihedral force constants. Default is 1.0
        force_parameterization : bool
            In case the template of a molecule already exists, its
            parameterization will be skipped unless force_parameterization
            is set to true. Default is False

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

        # Supress OpenForceField toolkit warnings
        import logging
        logging.getLogger().setLevel(logging.ERROR)

        self.dataset_name = dataset_name
        self.output_path = output_path
        self.n_proc = n_proc
        self.geometry_selection = geometry_selection
        self.PELE_exec = PELE_exec
        self.PELE_src = PELE_src
        self.PELE_license = PELE_license
        self.forcefield = forcefield
        self.charge_method = charge_method
        self.distort_bonds = distort_bonds
        self.distort_torsions = distort_torsions
        self.distort_dihedrals = distort_dihedrals
        self.distort_value_for_bonds = distort_value_for_bonds
        self.distort_value_for_torsions = distort_value_for_torsions
        self.distort_value_for_dihedrals = distort_value_for_dihedrals
        self.random_distort = random_distort
        self.bond_constants_factor = bond_constants_factor
        self.torsion_constants_factor = torsion_constants_factor
        self.dihedral_constants_factor = dihedral_constants_factor
        self.force_parameterization = force_parameterization

    def _distort_molecule(self, molecule_to_distort, output_path,
                          seed):
        """
        It distorts the input molecule, according to the supplied
        distortion settings.

        Parameters
        ----------
        molecule_to_distort : an peleffy.topology.Molecule object
            The molecule to distort
        output_path : str
            The path used as output in the benchmark
        seed : int
            Seed for the pseudo-random generator

        Returns
        -------
        output_file : str
            The path to the resulting PDB file with the distorted
            molecule
        """
        from copy import deepcopy

        # Make a copy of the input molecule
        mol = deepcopy(molecule_to_distort)

        if self.distort_bonds:
            from peleffybenchmarktools.structure import DistortBonds

            distorter = DistortBonds(mol, seed)
            if self.random_distort:
                distorted_mol = distorter.randomly(
                    self.distort_value_for_bonds)
            else:
                distorted_mol = distorter.fixed(
                    self.distort_value_for_bonds)
            mol._rdkit_molecule = distorted_mol

        if self.distort_torsions:
            from peleffybenchmarktools.structure import DistortAngles

            if seed is not None:
                seed += 1
            distorter = DistortAngles(mol, seed)
            if self.random_distort:
                distorted_mol = distorter.randomly(
                    self.distort_value_for_torsions)
            else:
                distorted_mol = distorter.fixed(
                    self.distort_value_for_torsions)
            mol._rdkit_molecule = distorted_mol

        if self.distort_dihedrals:
            from peleffybenchmarktools.structure import DistortDihedrals

            if seed is not None:
                seed += 1
            distorter = DistortDihedrals(mol, seed)
            if self.random_distort:
                distorted_mol = distorter.randomly(
                    self.distort_value_for_dihedrals)
            else:
                distorted_mol = distorter.fixed(
                    self.distort_value_for_dihedrals)

        from rdkit import Chem
        import os
        os.makedirs(os.path.join(output_path, mol.name), exist_ok=True)
        output_file = os.path.join(output_path, mol.name,
                                   'perturbed.pdb')
        Chem.rdmolfiles.MolToPDBFile(distorted_mol, output_file)

        return output_file

    def _get_molecule_minimized(self, output_path, pdb_path, seed=None):
        """
        It minimizes the molecule using PELE's minimization.

        Parameters
        ----------
        output_path : str
            The path where to run the PELE simulation and save the output
        pdb_path : str
            The path to the PDB of the molecule to minimize
        seed : int
            Seed for the pseudo-random generator. Default is None

        Returns
        -------
        output_file : str
            The path to the PELE minimized PDB file. It is None if PELE
            exits with a non-zero code
        """
        import os
        from peleffy.topology import Molecule
        from peleffybenchmarktools.utils.pele import PELEMinimization
        from peleffy.utils import OutputPathHandler
        from peleffy.forcefield import ForceFieldSelector
        from contextlib import suppress

        # Load and parameterize the molecule
        try:
            mol = Molecule(pdb_path)
            ff_selector = ForceFieldSelector()
            forcefield = ff_selector.get_by_name(self.forcefield)
            mol.set_forcefield(forcefield)

            output_handler = OutputPathHandler(
                mol,
                output_path=os.path.join(output_path, mol.name),
                as_datalocal=True)

            impact_output_path = output_handler.get_impact_template_path()
            solvent_output_path = output_handler.get_solvent_template_path()

            if (not os.path.exists(impact_output_path)
                    or not os.path.exists(solvent_output_path)
                    or self.force_parameterization):
                mol.parameterize(charge_method=self.charge_method)
                self._apply_parameter_factors(mol)
                with suppress(OSError):
                    os.remove(impact_output_path)
                    os.remove(solvent_output_path)

            # Distort molecule, if it is the case
            distorted_molecule_path = None
            if (self.distort_bonds or self.distort_torsions
                    or self.distort_dihedrals):
                distorted_molecule_path = self._distort_molecule(mol,
                                                                 output_path,
                                                                 seed)

            if 'openff' in self.forcefield.lower():
                forcefield_name = 'OpenForceField'
            else:
                forcefield_name = 'OPLS2005'

            # Runs a PELE Minimization
            pele_minimization = PELEMinimization(
                PELE_exec=self.PELE_exec,
                PELE_src=self.PELE_src,
                PELE_license=self.PELE_license,
                output_path=output_path,
                forcefield=forcefield_name)
            output_file = pele_minimization.run(
                mol, output_file='PELE_output.txt',
                pdb_path=distorted_molecule_path,
                force_parameterization=False,
                forcefield=self.forcefield,
                charge_method=self.charge_method)

            # Return not the path to the PELE output file but the path to
            # the PELE minimized PDB file
            # Be careful, we are using a hardcoded name
            return os.path.join(os.path.dirname(output_file), 'minimized.pdb')

        except Exception as e:
            print('Skipping minimization for molecule {}: '.format(pdb_path)
                  + str(e))

            return None

    def _get_molecule_minimized_with_openmm(self, output_path, pdb_paths,
                                            smiles_tags, index,
                                            seed=None):
        """
        It minimizes the molecule using OpenMM's minimization.

        Parameters
        ----------
        output_path : str
            The path where to run the PELE simulation and save the output
        pdb_paths : list[str]
            The list of paths to the PDBs of the molecules to minimize
        smiles_tags : list[str]
            The list of smiles tags of the molecules to minimize
        index : int
            The index that belongs to the molecule that is going to be
            minimized in the current iteration
        seed : int
            Seed for the pseudo-random generator. Default is None

        Returns
        -------
        output_file : str
            The path to the PELE minimized PDB file. It is None if PELE
            exits with a non-zero code
        """
        import os
        from openforcefield.topology import Molecule, Topology
        from openforcefield.typing.engines.smirnoff import ForceField
        from simtk.openmm.app import PDBFile
        from simtk import openmm, unit

        smiles_tag = smiles_tags[index]
        pdb_path = pdb_paths[index]

        mol = Molecule.from_pdb_and_smiles(pdb_path, smiles_tag)
        pdbfile = PDBFile(pdb_path)

        omm_topology = pdbfile.topology
        off_topology = Topology.from_openmm(omm_topology, unique_molecules=[mol])

        forcefield = ForceField('openff_unconstrained-1.2.1.offxml')

        system = forcefield.create_openmm_system(off_topology)

        time_step = 2 * unit.femtoseconds  # simulation timestep
        temperature = 300 * unit.kelvin  # simulation temperature
        friction = 1 / unit.picosecond  # collision rate
        integrator = openmm.LangevinIntegrator(temperature, friction, time_step)

        simulation = openmm.app.Simulation(omm_topology, system, integrator)

        positions = pdbfile.getPositions()
        simulation.context.setPositions(positions)

        simulation.minimizeEnergy()

        with open(os.path.join(output_path, str(index + 1) + '.pdb')) as f:
            f.write(
                PDBFile.writeModel(
                    simulation.topology,
                    simulation.context.getState(
                        getPositions=True).getPositions()))

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
        with Pool(self.n_proc) as pool:
            min_pdb_paths = list(tqdm(pool.imap(parallel_function, pdb_paths),
                                      total=len(pdb_paths)))

        # Filter out None pdb paths that belong to failing builds
        min_pdb_paths = [p for p in min_pdb_paths if p is not None]

        return min_pdb_paths

    def minimize_structures_with_openmm(self, pdb_paths):
        """
        It runs an OpenMM minimization for each of the PDB files from the
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
        from peleffybenchmarktools.utils import QCPortal

        # Handles output paths
        current_output = os.path.join(os.getcwd(), self.output_path, 'OpenMM')
        os.makedirs(current_output, exist_ok=True)

        qcportal = QCPortal()
        smiles_tags = [item['canonical_isomeric_explicit_hydrogen_smiles']
                       for item in
                       list(qcportal.get_data('OpenFF Optimization Set 1',
                                              'OptimizationDataset').values())]
        parallel_function = partial(self._get_molecule_minimized_with_openmm,
                                    current_output, pdb_paths, smiles_tags)

        # Loads the molecules from the Dataset and runs a PELEMinimization
        print(' - Minimizing molecules')
        with Pool(self.n_proc) as pool:
            min_pdb_paths = list(tqdm(pool.imap(parallel_function,
                                                range(0, len(pdb_paths))),
                                      total=len(pdb_paths)))

        # Filter out None pdb paths that belong to failing builds
        min_pdb_paths = [p for p in min_pdb_paths if p is not None]

        return min_pdb_paths

    def load_pdbs_from(self, search_path):
        """
        Given a path, it loads all the PDB files that it finds in it.

        Parameters
        ----------
        search_path : str
            The path where the PDB files will be searched at

        Returns
        -------
        pdb_paths : list[str]
            The list of paths pointing to the PDB files that are found
        """
        from glob import glob

        pdb_paths = list()
        for pdb_path in glob(search_path + '/*.pdb'):
            pdb_paths.append(pdb_path)

        if len(pdb_paths) == 0:
            print(' - Warning: no PDB file found at ' + search_path)

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
            if labeling1 == 'folder':
                name1 = os.path.basename(os.path.dirname(pdb_file1))
            else:
                name1 = os.path.splitext(os.path.basename(pdb_file1))[0]
            for pdb_file2 in set2:
                if labeling2 == 'folder':
                    name2 = os.path.basename(os.path.dirname(pdb_file2))
                else:
                    name2 = os.path.splitext(os.path.basename(pdb_file2))[0]
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
                     labeling1='file', labeling2='file',
                     output_name='rmsd'):
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
        output_name : str
            The name to use with the output files. Default is 'rmsd'
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
        df.to_csv(os.path.join(self.output_path,
                               '{}.csv'.format(output_name)))

        # Plots an histogram of the computed RMSD values
        plt.figure(figsize=(7, 5))
        plt.hist(rmsd_results, bins=10, rwidth=0.8, align='mid',
                 range=(0, 0.5), color='gray')
        plt.xlabel('RMSD')
        plt.ylabel('Frequency')
        plt.title('Structural Histogram')
        plt.savefig(os.path.join(self.output_path,
                                 '{}.png'.format(output_name)))

    def _apply_parameter_factors(self, mol):
        """
        It applies the parameter factors to the molecule.

        Parameters
        ----------
        mol : an peleffy.topology.Molecule object
            The molecule whose parameters will be customized
        """
        if self.bond_constants_factor != 1.0:
            for bond in mol.bonds:
                bond._spring_constant *= self.bond_constants_factor

        if self.torsion_constants_factor != 1.0:
            for torsion in mol.angles:
                torsion._spring_constant *= self.torsion_constants_factor

        if self.dihedral_constants_factor != 1.0:
            for dihedral in mol.propers:
                dihedral._constant *= self.dihedral_constants_factor

    def _get_bond_differences(self, data):
        """
        It calculates the mean bond differences of each linked pair of
        PDB files.

        Parameters
        ----------
        data : tuple[int, tuple[str, str]]
            The data to compute the bond differences of a linked pair

        Returns
        -------
        mean_difference : float
            The mean difference of bond lenghts for the current linked
            pair of PDB files
        """
        import numpy as np
        from peleffy.topology import Molecule
        from rdkit.Chem import rdMolTransforms
        from rdkit import Chem

        idx, (pdb_file1, pdb_file2) = data

        try:
            mol1 = Molecule(pdb_file1)
            mol1.parameterize('openff_unconstrained-1.2.1.offxml',
                              charge_method='gasteiger')
            rdkit_mol1 = mol1.rdkit_molecule
            conformer1 = rdkit_mol1.GetConformer()

            conformer2 = Chem.rdmolfiles.MolFromPDBFile(
                pdb_file2, proximityBonding=False,
                removeHs=False).GetConformer()

            bond_differences = []

            for bond in mol1.bonds:
                idx1 = bond.atom1_idx
                idx2 = bond.atom2_idx

                if rdkit_mol1.GetBondBetweenAtoms(idx1, idx2).IsInRing():
                    continue

                bond1 = rdMolTransforms.GetBondLength(conformer1,
                                                      idx1, idx2)

                bond2 = rdMolTransforms.GetBondLength(conformer2,
                                                      idx1, idx2)

                bond_differences.append(abs(bond1 - bond2))

            if len(bond_differences) > 0:
                mean_difference = np.mean(bond_differences)
            else:
                mean_difference = None

            return mean_difference

        except Exception as e:
            print('Skipping mean bond length computation between \''
                  + pdb_file1 + '\' \'' + pdb_file2 + '\': '
                  + str(e))

            return None

    def compute_bond_differences(self, paths_set1, paths_set2,
                                 labeling1='file', labeling2='file',
                                 output_name='bond_differences'):
        """
        For a collection of structures stored in two sets, it saves a
        CSV file with a dictionary of the bond lengths comparison between
        PELE and QM minimized structures. It also generates an histogram
        of the computed differences.

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
        output_name : str
            The name to use with the output files. Default is
            'bond_differences'
        """
        import os
        import pandas as pd
        import matplotlib.pyplot as plt
        from multiprocessing import Pool
        from tqdm import tqdm

        # Find links between the PDB paths from each set
        links = self._link_pdb_paths(paths_set1, paths_set2, labeling1,
                                     labeling2)

        # Computes the RMSD between PELE and QM minimized structures
        d = {}
        bond_difference_means = []
        with Pool(self.n_proc) as pool:
            results = list(tqdm(pool.imap(self._get_bond_differences,
                                          links.items()),
                                total=len(links)))

        for idx, mean_difference in zip(links, results):
            if mean_difference is None:
                continue
            d.update({idx: mean_difference})
            bond_difference_means.append(mean_difference)

        # Writes out a CSV file with the dictionary of the RMSD results.
        df = pd.DataFrame(d.items(),
                          columns=['Ligand ID',
                                   'Mean bond length differences'])
        df.to_csv(os.path.join(self.output_path,
                               '{}.csv'.format(output_name)))

        # Plots an histogram of the computed RMSD values
        plt.figure(figsize=(7, 5))
        plt.hist(bond_difference_means, bins=10, rwidth=0.8, align='mid',
                 color='gray')
        plt.xlabel('Mean bond length difference (Å)')
        plt.ylabel('Frequency')
        plt.title('Structural Histogram')
        plt.savefig(os.path.join(self.output_path,
                                 '{}.png'.format(output_name)))

    def _get_torsion_differences(self, data):
        """
        It calculates the mean torsion differences of each linked pair of
        PDB files.

        Parameters
        ----------
        data : tuple[int, tuple[str, str]]
            The data to compute the torsion differences of a linked pair

        Returns
        -------
        mean_difference : float
            The mean difference of torsion angles for the current linked
            pair of PDB files
        """
        import numpy as np
        from peleffy.topology import Molecule
        from rdkit.Chem import rdMolTransforms
        from rdkit import Chem

        idx, (pdb_file1, pdb_file2) = data

        try:
            mol1 = Molecule(pdb_file1)
            mol1.parameterize('openff_unconstrained-1.2.1.offxml',
                              charge_method='gasteiger')
            rdkit_mol1 = mol1.rdkit_molecule
            conformer1 = rdkit_mol1.GetConformer()

            conformer2 = Chem.rdmolfiles.MolFromPDBFile(
                pdb_file2, proximityBonding=False,
                removeHs=False).GetConformer()

            torsion_differences = []

            for angle in mol1.angles:
                idx1 = angle.atom1_idx
                idx2 = angle.atom2_idx
                idx3 = angle.atom3_idx

                if rdkit_mol1.GetBondBetweenAtoms(idx1, idx2).IsInRing():
                    continue

                if rdkit_mol1.GetBondBetweenAtoms(idx2, idx3).IsInRing():
                    continue

                angle1 = rdMolTransforms.GetAngleDeg(conformer1,
                                                     idx1, idx2, idx3)

                angle2 = rdMolTransforms.GetAngleDeg(conformer2,
                                                     idx1, idx2, idx3)

                torsion_differences.append(abs(angle1 - angle2))

            if len(torsion_differences) > 0:
                mean_difference = np.mean(torsion_differences)
            else:
                mean_difference = None

            return mean_difference

        except Exception as e:
            print('Skipping mean torsion computation between \''
                  + pdb_file1 + '\' \'' + pdb_file2 + '\': '
                  + str(e))

            return None

    def compute_torsion_differences(self, paths_set1, paths_set2,
                                    labeling1='file', labeling2='file',
                                    output_name='torsion_differences'):
        """
        For a collection of structures stored in two sets, it saves a
        CSV file with a dictionary of the torsion comparison between
        PELE and QM minimized structures. It also generates an histogram
        of the computed differences.

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
        output_name : str
            The name to use with the output files. Default is
            'torsion_differences'
        """
        import os
        import pandas as pd
        import matplotlib.pyplot as plt
        from multiprocessing import Pool
        from tqdm import tqdm

        # Find links between the PDB paths from each set
        links = self._link_pdb_paths(paths_set1, paths_set2, labeling1,
                                     labeling2)

        # Computes the RMSD between PELE and QM minimized structures
        d = {}
        torsion_difference_means = []
        with Pool(self.n_proc) as pool:
            results = list(tqdm(pool.imap(self._get_torsion_differences,
                                          links.items()),
                                total=len(links)))

        for idx, mean_difference in zip(links, results):
            if mean_difference is None:
                continue
            d.update({idx: mean_difference})
            torsion_difference_means.append(mean_difference)

        # Writes out a CSV file with the dictionary of the RMSD results.
        df = pd.DataFrame(d.items(),
                          columns=['Ligand ID',
                                   'Torsion angle differences'])
        df.to_csv(os.path.join(self.output_path,
                               '{}.csv'.format(output_name)))

        # Plots an histogram of the computed RMSD values
        plt.figure(figsize=(7, 5))
        plt.hist(torsion_difference_means, bins=10, rwidth=0.8, align='mid',
                 color='gray')
        plt.xlabel('Mean torsion angle difference (degrees)')
        plt.ylabel('Frequency')
        plt.title('Structural Histogram')
        plt.savefig(os.path.join(self.output_path,
                                 '{}.png'.format(output_name)))

    def _get_dihedral_differences(self, data):
        """
        It calculates the mean dihedral differences of each linked pair of
        PDB files.

        Parameters
        ----------
        data : tuple[int, tuple[str, str]]
            The data to compute the dihedral differences of a linked pair

        Returns
        -------
        mean_difference : float
            The mean difference of dihedral angles for the current linked
            pair of PDB files
        """
        import numpy as np
        from peleffy.topology import Molecule
        from rdkit.Chem import rdMolTransforms
        from rdkit import Chem

        idx, (pdb_file1, pdb_file2) = data

        try:
            mol1 = Molecule(pdb_file1)
            mol1.parameterize('openff_unconstrained-1.2.1.offxml',
                              charge_method='gasteiger')
            rdkit_mol1 = mol1.rdkit_molecule
            conformer1 = rdkit_mol1.GetConformer()

            conformer2 = Chem.rdmolfiles.MolFromPDBFile(
                pdb_file2, proximityBonding=False,
                removeHs=False).GetConformer()

            already_visited = set()
            dihedral_differences = []

            for dihedral in mol1.propers:
                idx1 = dihedral.atom1_idx
                idx2 = dihedral.atom2_idx
                idx3 = dihedral.atom3_idx
                idx4 = dihedral.atom4_idx
                if (idx2, idx3) in already_visited:
                    continue
                else:
                    already_visited.add((idx2, idx3))
                    already_visited.add((idx3, idx2))

                if rdkit_mol1.GetBondBetweenAtoms(idx2, idx3).IsInRing():
                    continue

                dihedral1 = rdMolTransforms.GetDihedralDeg(conformer1,
                                                           idx1, idx2,
                                                           idx3, idx4)

                dihedral2 = rdMolTransforms.GetDihedralDeg(conformer2,
                                                           idx1, idx2,
                                                           idx3, idx4)

                if dihedral1 < 0:
                    dihedral1 += 360
                if dihedral1 >= 360:
                    dihedral1 -= 360
                if dihedral2 < 0:
                    dihedral2 += 360
                if dihedral2 >= 360:
                    dihedral2 -= 360

                dihedral_differences.append(min([abs(dihedral1 - dihedral2),
                                                 abs(dihedral1 + 360 - dihedral2),
                                                 abs(dihedral1 - dihedral2 - 360)]))

            if len(dihedral_differences) > 0:
                mean_difference = np.mean(dihedral_differences)
            else:
                mean_difference = None

            return mean_difference

        except Exception as e:
            print('Skipping mean dihedral computation between \''
                  + pdb_file1 + '\' \'' + pdb_file2 + '\': '
                  + str(e))

            return None

    def compute_dihedral_differences(self, paths_set1, paths_set2,
                                     labeling1='file', labeling2='file',
                                     output_name='dihedral_differences'):
        """
        For a collection of structures stored in two sets, it saves a
        CSV file with a dictionary of the dihedral comparison between
        PELE and QM minimized structures. It also generates an histogram
        of the computed differences.

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
        output_name : str
            The name to use with the output files. Default is
            'dihedral_differences'
        """
        import os
        import pandas as pd
        import matplotlib.pyplot as plt
        from multiprocessing import Pool
        from tqdm import tqdm

        # Find links between the PDB paths from each set
        links = self._link_pdb_paths(paths_set1, paths_set2, labeling1,
                                     labeling2)

        # Computes the RMSD between PELE and QM minimized structures
        d = {}
        dihedral_difference_means = []
        with Pool(self.n_proc) as pool:
            results = list(tqdm(pool.imap(self._get_dihedral_differences,
                                          links.items()),
                                total=len(links)))

        for idx, mean_difference in zip(links, results):
            if mean_difference is None:
                continue
            d.update({idx: mean_difference})
            dihedral_difference_means.append(mean_difference)

        # Writes out a CSV file with the dictionary of the RMSD results.
        df = pd.DataFrame(d.items(),
                          columns=['Ligand ID',
                                   'Dihedral angle differences'])
        df.to_csv(os.path.join(self.output_path,
                               '{}.csv'.format(output_name)))

        # Plots an histogram of the computed RMSD values
        plt.figure(figsize=(7, 5))
        plt.hist(dihedral_difference_means, bins=10, rwidth=0.8, align='mid',
                 color='gray')
        plt.xlabel('Mean dihedral angle difference (degrees)')
        plt.ylabel('Frequency')
        plt.title('Structural Histogram')
        plt.savefig(os.path.join(self.output_path,
                                 '{}.png'.format(output_name)))

    def energetic_difference(self, minimized_pdb_paths,
                             output_name='energies'):
        """
        For a collection of structures stored in two sets, it computes
        and displays the distribution of energetic differences in an
        histogram.

        Parameters
        ----------
        minimized_pdb_paths : list[str]
            The list of minimized PDB paths
        output_name : str
            The name to use with the output files. Default is 'energies'
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
        df.to_csv(os.path.join(self.output_path, '{}.csv'.format(output_name)))

        # Plots an energetic histogram
        plt.figure(figsize=(7, 5))
        plt.hist(energetic_diff, bins=10, rwidth=0.8, align='mid', color='green')
        plt.xlabel('Energetic difference (kcal/mol)')
        plt.ylabel('Frequency')
        plt.title('Energetic Histogram')
        plt.savefig(os.path.join(self.output_path, '{}.png'.format(output_name)))

    def distort_structures(self, pdb_paths,
                           bond_distortion_range=0,
                           angle_distortion_range=0,
                           dihedral_distortion_range=0,
                           seed=None):
        """
        Given a set of PDB files and distortion ranges, it applies
        structural distortions and saves the obtained results.

        Parameters
        ----------
        pdb_paths : list[str]
            The list of paths to the PDB files to distort
        bond_distortion_range : float
            The upper limit that defines the range of values that
            can be used to distort bond lengths of the structure
        angle_distortion_range : float
            The upper limit that defines the range of values that
            can be used to distort torsional angles of the structure
        dihedral_distortion_range : float
            The upper limit that defines the range of values that
            can be used to distort dihedral angles of the structure
        seed : int
            Seed for the pseudo-random generator

        Returns
        -------
        output_files : list[str]
            The list of paths to the resulting PDB files with distorted
            structures
        """
        import os
        from tqdm import tqdm
        from multiprocessing import Pool
        from functools import partial

        output_path = os.path.join(self.output_path, 'distorted')
        os.makedirs(output_path, exist_ok=True)

        parallel_function = partial(self._distort_structure,
                                    output_path,
                                    bond_distortion_range,
                                    angle_distortion_range,
                                    dihedral_distortion_range,
                                    seed)

        with Pool(self.n_proc) as p:
            output_files = list(tqdm(p.imap(parallel_function, pdb_paths),
                                     total=len(pdb_paths)))

        return output_files

    def _distort_structure(self, output_path,
                           bond_distortion_range,
                           angle_distortion_range,
                           dihedral_distortion_range,
                           seed,
                           pdb_to_distort):
        """
        Given a set of distortion ranges, it distorts a PDB structure
        and saves the result to the supplied output path.

        Parameters
        ----------
        output_path : str
            The path where to save the resulting distorted structure
        bond_distortion_range : float
            The upper limit that defines the range of values that
            can be used to distort bond lengths of the structure
        angle_distortion_range : float
            The upper limit that defines the range of values that
            can be used to distort torsional angles of the structure
        dihedral_distortion_range : float
            The upper limit that defines the range of values that
            can be used to distort dihedral angles of the structure
        seed : int
            Seed for the pseudo-random generator
        pdb_to_distort : str
            The path to the PDB to distort

        Returns
        -------
        output_file : str
            The path to the resulting PDB file with the distorted
            structure
        """
        import os
        from peleffy.topology import Molecule
        from rdkit import Chem
        from peleffybenchmarktools.structure import (DistortBonds,
                                                     DistortAngles,
                                                     DistortDihedrals)

        idx = os.path.splitext(os.path.basename(pdb_to_distort))[0]
        molecule = Molecule(pdb_to_distort)
        molecule.parameterize('openff_unconstrained-1.2.0.offxml',
                              charge_method='gasteiger')

        if bond_distortion_range > 0:
            distort = DistortBonds(molecule, seed)
            distorted_mol = distort.randomly(
                range=bond_distortion_range)
            molecule._rdkit_molecule = distorted_mol

        if angle_distortion_range > 0:
            distort = DistortAngles(molecule, seed)
            distorted_mol = distort.randomly(
                range=angle_distortion_range)
            molecule._rdkit_molecule = distorted_mol

        if dihedral_distortion_range > 0:
            distort = DistortDihedrals(molecule, seed)
            distorted_mol = distort.randomly(
                range=dihedral_distortion_range)

        output_file = os.path.join(output_path,
                                   '{}.pdb'.format(idx))

        Chem.rdmolfiles.MolToPDBFile(distorted_mol, output_file)

        return output_file
