"""
Main script for the hydration free energies benchmark.

 
Example:
-------------------
python main.py OFF -solvent OBC -c gasteiger -p -o OFF_out1

"""

import os
import argparse 
import pandas as pd
from offpele.topology import Molecule, RotamerLibrary
from offpele.template import Impact
from offpele.solvent import OBC2
from offpele.main import handle_output_paths
import DiffEnergies as DiffEnergies
from methods import method_OFF, method_OPLS, method_OFFOPLS



def parse_args():
    """
        Parse command line arguments
        :returns: object -- Object containing command line options
    """
    parser = argparse.ArgumentParser(description="Hydration Energies Benchmark")
    parser.add_argument("method", type=str, help="Method you want to use: OPLS, OFF, OPLS-OFF.")
    parser.add_argument("-solvent", "--solvent", type=str, default="", help="Implicid solvent.")
    parser.add_argument("-nb", "--non_bonds", action='store_true', help="Method for obtaining non bonds parameters.")
    parser.add_argument("-ba", "--bonds_angles", action='store_true', help="Method for obtaining  bonds and angles parameters.")
    parser.add_argument("-c", "--charges", type=str, default="", help="Charges method")
    parser.add_argument("-p", "--plots", action='store_true' ,help="For plotting the results." )
    parser.add_argument("-o", "--output_folder", type=str, default="./out", help="Name of the folder where the results will be saved.")
    args = parser.parse_args()
    return args.method, args.solvent, args.non_bonds, args.bonds_angles, args.charges, args.plots, args.output_folder


def check_args(method, solvent, non_bonds, bonds_angles, charges):
	"""
		It handels the different args needed for each method and returns an error if any are missing. 
	"""
	import sys

	if  (method != 'OPLS') and (method != 'OFF') and (method != 'OFFOPLS'):
		sys.exit('Error: Invalid method.') 
	if method == 'OPLS':
		if solvent == '': 
			sys.exit('Error: For OPLS you need to indicate the solvent.')
	if method == 'OFF':
		if solvent == '':
			sys.exit('Error: For OFF you need to indicate the solvent.')
		if charges == '': 
			sys.exit('Error: For OFF you need to indicate the charges method.')
	if method == 'OFFOPLS':
		if solvent == '': 
			sys.exit('Error: For OPLS-OFF you need to indicate the solvent.')
		if non_bonds == '': 
			sys.exit('Error: For OPLS-OFF you need to indicate which non bonding parameters will be used.')
		if bonds_angles == '': 
			sys.exit('Error: For OPLS-OFF you need to indicate which bonding and angles parameters will be used.')
		if charges == '': 
			sys.exit('Error: For OPLS-OFF you need to indicate the charges method.')

def read_dataset(PATH_TO_FREESOLV_DATABASE):
    """
        It reads the FreeSolv database and returns the lists of the needed values.
    """
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
    smiles_tags = freesolv_db['SMILES'].to_list()
    experimental_v = freesolv_db['experimental value'].to_list()
    return compound_ids, smiles_tags, experimental_v



def save_output(out_folder, energies): 
	"""
		It saves the output results.  
	"""
	df = pd.DataFrame(energies, columns = ['CID','Energetic Difference', 'Experimental value'])
	df.to_csv(os.path.join(out_folder,'results.txt'))

def plot_results(out_folder,differences, experimental_values):
	"""
		It generates and histogram and a regression for the comparision between the experimental values and the computed for the hydration free energy.
	"""
	import matplotlib.pyplot as plt 
	import numpy as np
	import pandas 
	from collections import Counter 
	import os 

	# Create a directory for the plots
	os.makedirs(os.path.join(out_folder,'Plots'), exist_ok=True)

	# Computes the fit for the regresion
	y = np.array(differences)
	x = np.array(experimental_values)
	coef = np.polyfit(x,y, deg = 1)

	#Histogram
	abs_val = abs(x - y)
	bins = np.arange(0,10,0.75)
	plt.figure()
	plt.grid(axis='y', alpha=0.75)
	plt.ylabel('Frequency');
	plt.xlabel('Absolute difference (kcal/mol)')
	plt.title('Absolute difference')
	_, bins, patches = plt.hist(np.clip(abs_val, bins[0], bins[-1]),                        
                                bins=bins, color=['#0504aa'], alpha=0.7, rwidth=0.8)
	plt.savefig(os.path.join(os.path.join(out_folder,'Plots'), 'histogram.png'))

	#Regression
	poly1d_fn = np.poly1d(coef)
	plt.figure()
	plt.title('Experimental value vs. Energetic Difference')
	plt.ylabel('Energetic difference (kcal/mol)')
	plt.xlabel('Experimental value (kcal/mol)')
	plt.plot(x,y, 'yo', x, poly1d_fn(x), '--k')
	plt.savefig(os.path.join(os.path.join(out_folder,'Plots'), 'regression.png'))


def main(method, solvent, non_bonds, bonds_angles, charges, plots, out_folder):
	
	CHARGES_METHOD = charges
	PATH_TO_FREESOLV_DATABASE = '/home/lauramalo/repos/offpele-benchmarks/benchmarks/data/FreeSolv/FreeSolv0.52.txt'
	OFF_FORCEFIELD = 'openff_unconstrained-1.2.0.offxml'
	
	PELE_EXEC = '/home/municoy/builds/PELE/PELE-repo_serial/PELE-1.6'
	PELE_SRC = '/home/municoy/repos/PELE-repo/'
	
	SCHRODINGER_SRC = '/data/general_software/schrodinger2017-2/'
	PLOPROTTEMP_SRC = '/home/municoy/repos/PlopRotTemp/main.py'
	IMPACT_TEMPLATE_PATH = 'DataLocal/Templates/OPLS2005/HeteroAtoms/'
	ROTAMER_LIBRARY_PATH = 'DataLocal/LigandRotamerLibs/'

	# It checks that the args are correct for each method
	check_args(method, solvent, non_bonds, bonds_angles, charges)

	# It reads the dataset FreeSolv
	compound_ids, smiles_tags, experimental_v = read_dataset(PATH_TO_FREESOLV_DATABASE)

	#It runs the selected method
	if method == 'OFF':

		VACUUM_CF = '/home/lauramalo/repos/offpele-benchmarks/benchmarks/solvent/Conf/VACUUM_minimization.conf'
		OBC_CF = '/home/lauramalo/repos/offpele-benchmarks/benchmarks/solvent/Conf/OBC_minimization.conf'

		energies, differences, experimental_values = \
			method_OFF(out_folder, compound_ids, smiles_tags, experimental_v, OFF_FORCEFIELD, CHARGES_METHOD, PELE_EXEC, VACUUM_CF, OBC_CF, PELE_SRC)

	if method == 'OPLS':

		VACUUM_CF = '/home/lauramalo/repos/offpele-benchmarks/benchmarks/solvent/Conf/OPLS_VACUUM_minimization.conf'
		OBC_CF = '/home/lauramalo/repos/offpele-benchmarks/benchmarks/solvent/Conf/OPLS_OBC_minimization.conf'
		VDGBNP_CF = '/home/lauramalo/repos/offpele-benchmarks/benchmarks/solvent/Conf/OPLS_VDGBNP_minimization.conf'

		energies, differences, experimental_values = \
			method_OPLS(out_folder,compound_ids, smiles_tags, experimental_v, solvent, PELE_EXEC, SCHRODINGER_SRC, PLOPROTTEMP_SRC, PELE_SRC, IMPACT_TEMPLATE_PATH, ROTAMER_LIBRARY_PATH, VACUUM_CF, OBC_CF,VDGBNP_CF)
	
	if method == 'OFFOPLS':

		VACUUM_CF = '/home/lauramalo/repos/offpele-benchmarks/benchmarks/solvent/Conf/OPLS_VACUUM_minimization.conf'
		OBC_CF = '/home/lauramalo/repos/offpele-benchmarks/benchmarks/solvent/Conf/OPLS_VDGBNP_minimization.conf'
		
		energies, differences, experimental_values = \
			method_OFFOPLS(out_folder, compound_ids, smiles_tags, experimental_v, solvent, non_bonds, bonds_angles, OFF_FORCEFIELD, PELE_EXEC, CHARGES_METHOD,  PELE_SRC, VACUUM_CF, OBC_CF, IMPACT_TEMPLATE_PATH, ROTAMER_LIBRARY_PATH)

	# Saves the results in a file
	save_output(out_folder,energies)

	# If selected, it create the plots for the analysis of the results. 
	if plots: plot_results(out_folder, differences, experimental_values)

if __name__ == "__main__":
    method, solvent, non_bonds, bonds_angles, charges, plots, out_folder = parse_args()
    main(method, solvent, non_bonds, bonds_angles, charges, plots, out_folder)

