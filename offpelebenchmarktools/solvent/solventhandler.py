

class SolventBenchmark(object):
    
    def __init__(self,method):
        """
        It initialized an SolventBenchmark object. 

        Parameters
        ----------

        method: an offpelebenchmarktools.solvent.method object 
            The Method object that will be used to run the benchmark. 
        """

        self._method = method



    def _read_dataset(self,PATH_TO_FREESOLV_DATABASE):
        """
            It reads the FreeSolv database and returns the lists of the needed values.
        """
        import pandas as pd 
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

    def _get_method(self):
        print('Hola',self._method)
        method_type = self._method[0]
        charges = self._method[1]
        solvent =  self._method[2]
        non_bonds = self._method[3]
        bonds_angles = self._method[4]
        return method_type, charges, solvent, non_bonds, bonds_angles

    def run(self):
        from solventhandler import SolventBenchmark
        from methods import method_OFF, method_OPLS, method_OFFOPLS
        
        PATH_TO_FREESOLV_DATABASE = '/home/lauramalo/repos/offpele-benchmarks/benchmarks/data/FreeSolv/FreeSolv0.52.txt'
        OFF_FORCEFIELD = 'openff_unconstrained-1.2.0.offxml'
        
        PELE_EXEC = '/home/municoy/builds/PELE/PELE-repo_serial/PELE-1.6'
        PELE_SRC = '/home/municoy/repos/PELE-repo/'
        
        SCHRODINGER_SRC = '/data/general_software/schrodinger2017-2/'
        PLOPROTTEMP_SRC = '/home/municoy/repos/PlopRotTemp/main.py'
        IMPACT_TEMPLATE_PATH = 'DataLocal/Templates/OPLS2005/HeteroAtoms/'
        ROTAMER_LIBRARY_PATH = 'DataLocal/LigandRotamerLibs/'

        method_type, CHARGES_METHOD, solvent, non_bonds, bonds_angles = self._get_method()
        compound_ids, smiles_tags, experimental_v = self._read_dataset(PATH_TO_FREESOLV_DATABASE)

        
        #It runs the selected method
        if method_type == 'OFF':

            VACUUM_CF = '/home/lauramalo/repos/offpele-benchmarks/benchmarks/solvent/Conf/VACUUM_minimization.conf'
            OBC_CF = '/home/lauramalo/repos/offpele-benchmarks/benchmarks/solvent/Conf/OBC_minimization.conf'

            energies, differences, experimental_values = \
                method_OFF(out_folder, compound_ids, smiles_tags, experimental_v, OFF_FORCEFIELD, CHARGES_METHOD, PELE_EXEC, VACUUM_CF, OBC_CF, PELE_SRC)

        if method_type == 'OPLS':

            VACUUM_CF = '/home/lauramalo/repos/offpele-benchmarks/benchmarks/solvent/Conf/OPLS_VACUUM_minimization.conf'
            OBC_CF = '/home/lauramalo/repos/offpele-benchmarks/benchmarks/solvent/Conf/OPLS_OBC_minimization.conf'
            VDGBNP_CF = '/home/lauramalo/repos/offpele-benchmarks/benchmarks/solvent/Conf/OPLS_VDGBNP_minimization.conf'

            energies, differences, experimental_values = \
                method_OPLS(out_folder,compound_ids, smiles_tags, experimental_v, solvent, PELE_EXEC, SCHRODINGER_SRC, PLOPROTTEMP_SRC, PELE_SRC, IMPACT_TEMPLATE_PATH, ROTAMER_LIBRARY_PATH, VACUUM_CF, OBC_CF,VDGBNP_CF)
        
        if method_type == 'OFFOPLS':

            VACUUM_CF = '/home/lauramalo/repos/offpele-benchmarks/benchmarks/solvent/Conf/OPLS_VACUUM_minimization.conf'
            OBC_CF = '/home/lauramalo/repos/offpele-benchmarks/benchmarks/solvent/Conf/OPLS_VDGBNP_minimization.conf'
            
            energies, differences, experimental_values = \
                method_OFFOPLS(out_folder, compound_ids, smiles_tags, experimental_v, solvent, non_bonds, bonds_angles, OFF_FORCEFIELD, PELE_EXEC, CHARGES_METHOD,  PELE_SRC, VACUUM_CF, OBC_CF, IMPACT_TEMPLATE_PATH, ROTAMER_LIBRARY_PATH)
        return energies, differences, experimental_values


    def save_output(self, out_folder, energies): 
        """
            It saves the output results.  
        """
        df = pd.DataFrame(energies, columns = ['CID','Energetic Difference', 'Experimental value'])
        df.to_csv(os.path.join(out_folder,'results.txt'))

    def plot_results(self, energies = None, differences = None , experimental_values = None, file = None):
        """
            It generates and histogram and a regression for the comparision between the experimental values and the computed for the hydration free energy.
        """
        import matplotlib.pyplot as plt 
        import numpy as np
        import pandas as pd 
        from collections import Counter 
        import os 

        if file != None: 
            data = pd.read_csv(file, usecols=['CID', 'Energetic Difference', 'Experimental value'])
            y = data['Energetic Difference'].to_numpy()
            x = data['Experimental value'].to_numpy()
        else:
            y = np.array(differences)
            x = np.array(experimental_values)
        
        # Computes the fit for the regresion
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

        #Regression
        poly1d_fn = np.poly1d(coef)
        plt.figure()
        plt.title('Experimental value vs. Energetic Difference')
        plt.ylabel('Energetic difference (kcal/mol)')
        plt.xlabel('Experimental value (kcal/mol)')
        plt.plot(x,y, 'yo', x, poly1d_fn(x), '--k')




