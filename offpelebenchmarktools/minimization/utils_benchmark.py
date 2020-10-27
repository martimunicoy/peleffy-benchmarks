# Local imports
from getter import QCPortal
from minimize import Minimizer
from offpelebenchmarktools.utils.pele import PELEBaseJob, PELEMinimization
from offpele.topology import Molecule


# External imports
import argparse
import os 
import pandas as pd 
import glob
import re
import matplotlib.pyplot as plt

class MinimizationBenchmark(object): 
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
        self.dataset = dataset
        self.out_folder = out_folder

    def _get_molecule_minimized(self, index):
        """
            It minimized the molecule using PELE's minimization. 
        """
        
        #Load and parameterize the molecule 
        try: 
            mol = Molecule(os.path.join(os.getcwd(), self.out_folder,'QM/' '{}.pdb'.format(index + 1)))
            mol.parameterize('openff_unconstrained-1.2.1.offxml',
                charges_method='gasteiger')
            
            # Runs a PELE Minimization
            pele_minimization = PELEMinimization(
                PELE_exec='/home/municoy/builds/PELE/PELE-repo_serial/PELE-1.6',
                PELE_src='/home/municoy/repos/PELE-repo/',
                PELE_license='/home/municoy/builds/PELE/license')
            pele_minimization.run(mol)
        except: 
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
            It filter structures to keep only those that employ one of the OpenFF dihedrals with a non null phase constant.
            Output:
            ------------
                - It the structure employs one of the OpenFF dihedrals with a non null phase constant: return True.
                - Else: return False
        """

        # Filter out all non interesting molecules
        from offpele.topology import Molecule
        from simtk import unit

        var = False
        try:
            mol = Molecule(smiles=smiles_tag)
            mol.parameterize('openff_unconstrained-1.2.1.offxml',
                            charges_method='gasteiger')
            for p in mol.propers:
                if p.phase not in (unit.Quantity(0, unit.degree),
                                unit.Quantity(180, unit.degree)):
                    var = True
        except:
            var = False 
        return var

    def _get_dataset_structures(self, filter_angles):
        """
            It gets the Dataset from the QCPortal and saves in a folder the optimized molecules as PDB.
        """
        import qcportal as ptl

        # Get the optimization dataset
        client = ptl.FractalClient()
        ds = client.get_collection('OptimizationDataset', self.dataset)
        nmols = len(list(ds.data.records.keys()))

        # Initializes the QCPortal class
        qc_portal = QCPortal(n_proc=4)
        
        #Handles output paths
        os.mkdir(os.path.join(os.getcwd(), self.out_folder))
        set_folder = os.path.join(os.getcwd(), self.out_folder, 'QM')
        os.mkdir(set_folder)
        
        # Build the optimized molecule and save it as a PDB
        for index in range(nmols):
            if (index % 50) == 0 : print('{}/{} optimized molecules built'.format(index + 1, nmols))
            entry =  ds.get_entry(ds.df.index[index])
            if filter_angles: 
                if self._filter_structures(smiles_tag = entry.attributes.get('canonical_smiles')): 
                    qc_portal._parallel_struct_getter(ds, set_folder, index)
            else: 
                qc_portal._parallel_struct_getter(ds, set_folder, index)



    def run(self, filter_structures = False): 
        """
            It generates the folders for the optimized with QM structures to PDB files and these PDBs minnimized with PELE.
        """
        import shutil 

        # Obtain the dataset  
        self._get_dataset_structures(filter_angles = filter_structures)

        # Gets all the PDB files from the set
        pdb_files = glob.glob(os.path.join(os.path.join(self.out_folder, 'QM'), "*pdb"))

        # Loads the molecules from the Dataset and runs a PELEMinimization
        for pdb_file in pdb_files:
            _ ,index = os.path.split(pdb_file)
            index = re.sub('.pdb','',index)
            self._get_molecule_minimized(index = int(index))

        # Moves the output folder(created by the PELEMinimization) to the desired output folder
        shutil.move(os.path.join(os.getcwd(),'output'), os.path.join(os.getcwd(), self.out_folder, 'PELE'))


    def compute_RMSD(self):
        """
            For a collection of structures, 
            it saves a CSV file with a dictionary of the RMSD comparison between PELE and QM minimized structures.
            It generats an histogram of the computed RMSD values. 
        """

        import mdtraj as md

        # Gets all the structures from the output folder
        pdb_files = glob.glob(os.path.join(os.path.join(self.out_folder, 'QM'), "*pdb"))

        # Computes the RMSD between PELE and QM minimized structures
        d = {}
        rmsd_results = []
        for pdb_file in pdb_files: 
            _ ,index = os.path.split(pdb_file)
            index = re.sub('.pdb','',index)
            try:
                molQM = md.load_pdb(os.path.join(self.out_folder, 'QM/{}.pdb'.format(index)))
                molPELE = md.load_pdb(os.path.join(self.out_folder, 'PELE/{}/minimized.pdb'.format(index)))
                rsmd = md.rmsd(molQM, molPELE)[0]
                d.update({index :rsmd})
                rmsd_results.append(rsmd)
            except: 
                pass
 
        # Writes out a CSV file with the dictionary of the RMSD results.
        df = pd.DataFrame(d.items(), columns=['Ligand ID', 'RMSD'])
        df.to_csv(os.path.join(self.out_folder,'rmsd.csv'))

        #Plots an histogram of the computed RMSD values
        plt.figure(figsize=(10, 7))
        plt.hist(rmsd_results, bins = 10, rwidth = 0.8, align ='mid', range = (0,1), color = 'gray')
        plt.xlabel('RMSD')
        plt.ylabel('Frequency')
        plt.title('Structural Histogram')
        plt.savefig(os.path.join(self.out_folder,'rmsd.png'))

    def energetic_difference(self): 
        """
            For a ollection of structures, gets the energy before and after the minimization and computes its difference. 
            It generates a file with the energies and an energetic histogram. 
        """


        # Gets all the structures from the output folder
        pdb_files = glob.glob(os.path.join(os.path.join(self.out_folder, 'QM'), "*pdb"))

        # Computes the energetic difference between PELE and QM minimized structures
        d = {}
        energetic_diff =  []
        for pdb_file in pdb_files: 
            _ ,index = os.path.split(pdb_file)
            index = re.sub('.pdb','',index)
            try:
                file_path = os.path.join(self.out_folder, 'PELE/{}/PELE_output.txt'.format(index))
                data = self._parse_output_file(file = file_path)
                e_ini = data[-1].get('ENERGY VACUUM + CONSTRAINTS')
                e_fin = data[-2].get('ENERGY VACUUM + CONSTRAINTS')
                e_diff = float(e_fin) - float(e_ini)
                energetic_diff.append(e_diff)
                d.update({index: tuple((e_ini, e_fin, e_diff))})
            except: 
                pass

        # Writes out a CSV file with the dictionary of the energies results.
        df = pd.DataFrame(d.items(), columns=['Ligand ID', 'Energies'])
        df.to_csv(os.path.join(self.out_folder,'energies.csv'))

       #Plots an energetic histogram
        plt.figure(figsize=(10, 7))
        plt.hist(energetic_diff, bins = 10, rwidth = 0.8, align ='mid', color = 'green')
        plt.xlabel('Energetic difference (kcal/mol)')
        plt.ylabel('Frequency')
        plt.title('Energetic Histogram')       
        plt.savefig(os.path.join(self.out_folder,'energies.png'))



