# Local imports
from getter import QCPortal
from minimize import Minimizer
from offpelebenchmarktools.utils.pele import PELEBaseJob, PELEMinimization
from offpele.topology import Molecule

# External imports
import argparse
import os 
import pandas as pd 

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
        
        #Load  and parameterize the molecule 
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
                print('PHASE', p.phase)
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
        import glob
        import re

        # Obtain the dataset  
        self._get_dataset_structures(filter_angles = filter_structures)

        # Gets the number of files to minimize
        pdb_files = glob.glob(os.path.join(os.path.join(self.out_folder, 'QM'), "*pdb"))

        # Loads the molecules from the Dataset and runs a PELEMinimization
        for pdb_file in pdb_files:
            _ ,index = os.path.split(pdb_file)
            index = re.sub('.pdb','',index)
            self._get_molecule_minimized(index = int(index))

        # Creates the output folder
        shutil.move(os.path.join(os.getcwd(),'output'), os.path.join(os.getcwd(), self.out_folder, 'PELE'))

    def compute_RMSD(self):
        """
            For a collection of structures, 
            it saves a CSV file with a dictionary of the RMSD comparison between PELE and QM minimized structures.
        """

        from rdkit.Chem import rdMolAlign
        from rdkit.Chem import rdmolfiles

        import matplotlib.pyplot as plt
        import re 
        import glob

        # Computes the RMSD between ∫PELE and QM minimized structures
        pdb_files = glob.glob(os.path.join(os.path.join(self.out_folder, 'QM'), "*pdb"))
        d = {}
        rmsd_results = []
        for pdb_file in pdb_files: 
            _ ,index = os.path.split(pdb_file)
            index = re.sub('.pdb','',index)
            try:
                molQM = rdmolfiles.MolFromPDBFile(os.path.join(self.out_folder, 'QM/{}.pdb'.format(index)))
                molPELE = rdmolfiles.MolFromPDBFile(os.path.join(self.out_folder, 'PELE/{}/minimized.pdb'.format(index)))
                rsmd = rdMolAlign.AlignMol(molQM, molPELE)
                print(index, rsmd)
                d.update({index :rsmd})
                rmsd_results.append(rsmd)
            except Exception as e: 
                print('Skipping ligand {}. {}'.format(index,e))

        
        # Writes out a CSV file with the dictionary of the RMSD results.
        df = pd.DataFrame(d.items(), columns=['Ligand ID', 'RMSD'])
        df.to_csv(os.path.join(self.out_folder,'rmsd.csv'))

        #Plots an histogram of the computed RMSD values
        plt.figure(figsize=(10, 7))
        plt.hist(rmsd_results, rwidth = 0.8, align ='mid', color = 'green')
        plt.xlabel('RMSD')
        plt.ylabel('Frequency')
        plt.savefig(os.path.join(self.out_folder,'rmsd.png'))




