"""
This module contains classes and functions to validate the solvent
model using PELE with different force fields.
"""


from offpelebenchmarktools.utils import get_data_file_path
from offpelebenchmarktools.solvent.hydfreeenergycalculator import runner


FREESOLV_PATH = 'databases/FreeSolv0.52.txt'


class SolventBenchmark(object):
    def __init__(self, pele_exec, pele_src, pele_license,
                 ploprottemp_src, schrodinger_src,
                 charge_method='am1bcc', solvent='OBC',
                 opls_nonbonding=False, opls_bonds_angles=False,
                 n_proc=1, forcefield_name=None, forcefield=None):
        """
        It initialized an SolventBenchmark object.

        Parameters
        ----------
        PELE_exec : str
            Path to the PELE executable
        PELE_src : str
            Path to PELE source folder
        PELE_license : str
            Path to PELE license directory
        ploprottemp_src : str
            Path to PlopRotTemp source code
        schrodinger_src : str
            Path to Schrodinger source code
        charge_method : str
            The method to calculate partial charges
        solvent : str
            The solvent model to employ
        opls_nonbonding : bool
            Whether to use OPLS2005 to parameterize nonbonding terms or not
        opls_bonds_angles : bool
            Whether to use OPLS2005 to paramterize bonds and angles or not
        n_proc : int
            Number of parallel computing processors to employ. Default is 1
        forcefield_name : str
            The force field name to employ. Default is None
        forcefield : an offpele.forcefield._BaseForceField
            The forcefield representation to employ. Default is None
        """
        self.pele_exec = pele_exec
        self.pele_src = pele_src
        self.pele_license = pele_license
        self.charge_method = charge_method
        self.solvent = solvent
        self.opls_nonbonding = opls_nonbonding
        self.opls_bonds_angles = opls_bonds_angles
        self.ploprottemp_src = ploprottemp_src
        self.schrodinger_src = schrodinger_src
        self._n_proc = n_proc
        self.forcefield_name = forcefield_name
        self.forcefield = forcefield
        self._results = dict()

        # Deactivate offpele output
        from offpele.utils import Logger
        logger = Logger()
        logger.set_level('WARNING')

        # Supress OpenForceField toolkit warnings
        import logging
        logging.getLogger().setLevel(logging.ERROR)

    def _read_dataset(self):
        """
        It reads the FreeSolv database and returns the lists of the needed
        values.
        """
        import pandas as pd

        freesolv_path = get_data_file_path(FREESOLV_PATH)

        freesolv_db = pd.read_csv(freesolv_path, delimiter=';',
                                  skipinitialspace=True,
                                  skiprows=[0, 1, 2], header=0,
                                  names=['compound id', 'SMILES',
                                         'iupac name',
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

    def run(self, out_folder):
        compound_ids, smiles_tags, experimental_v = self._read_dataset()

        # It runs the selected method
        energies = runner(out_folder, compound_ids,
                          smiles_tags, experimental_v,
                          self.solvent, self.charge_method,
                          self.pele_exec, self.pele_src,
                          self.pele_license,
                          n_proc=self._n_proc,
                          forcefield_name=self.forcefield_name,
                          forcefield=self.forcefield)

        self.results['cids'] = list()
        self.results['differences'] = list()
        self.results['experimental_values'] = list()

        for cid, difference, experimental_value in energies:
            self.results['cids'].append(cid)
            self.results['differences'].append(difference)
            self.results['experimental_values'].append(experimental_value)

    def to_csv(self, out_folder):
        """It saves the output results as a csv file."""
        import pandas as pd
        import os

        df = pd.DataFrame(zip(self.results['cids'],
                              self.results['differences'],
                              self.results['experimental_values']),
                          columns=['CID', 'Energetic Difference',
                                   'Experimental value'])
        df.to_csv(os.path.join(out_folder))

    def from_csv(self, path_to_load):
        """It loads the results from a csv file."""
        import pandas as pd

        self._results = pd.read_csv(path_to_load).to_dict()

    def plot_results(self, energies=None, differences=None,
                     experimental_values=None):
        """
        It generates and histogram and a regression for the comparision
        between the experimental values and the computed for the hydration
        free energy.
        """
        import matplotlib.pyplot as plt
        import numpy as np

        y = np.array(self.results['differences'])
        x = np.array(self.results['experimental_values'])

        # Computes the fit for the regresion
        coef = np.polyfit(x, y, deg=1)

        # Histogram
        abs_val = abs(x - y)
        bins = np.arange(0, 10, 0.75)
        plt.figure()
        plt.grid(axis='y', alpha=0.75)
        plt.ylabel('Frequency')
        plt.xlabel('Absolute difference (kcal/mol)')
        plt.title('Absolute difference')
        _, bins, patches = plt.hist(np.clip(abs_val, bins[0], bins[-1]),
                                    bins=bins, color=['#0504aa'], alpha=0.7,
                                    rwidth=0.8)

        # Regression
        poly1d_fn = np.poly1d(coef)
        plt.figure()
        plt.title('Experimental value vs. Energetic Difference')
        plt.ylabel('Energetic difference (kcal/mol)')
        plt.xlabel('Experimental value (kcal/mol)')
        plt.plot(x, y, 'yo', x, poly1d_fn(x), '--k')

    @property
    def results(self):
        """The benchmark results."""
        return self._results

    def display_outlayers_above(self, threshold):
        """Displays the molecules that overcome the given threshold."""

        from offpele.topology import Molecule
        from IPython.display import display

        compound_ids, smiles_tags, _ = self._read_dataset()
        cid_to_smiles = dict()

        for cid, smiles_tag in zip(compound_ids, smiles_tags):
            cid_to_smiles[cid] = smiles_tag

        for cid, diff, expv in zip(self.results['cids'],
                                   self.results['differences'],
                                   self.results['experimental_values']):
            if abs(diff - expv) > 10:
                smiles = cid_to_smiles[cid]
                print('-' * (len(cid) + len(smiles) + 3))
                print(cid, '-', smiles)
                print('-' * (len(cid) + len(smiles) + 3))
                mol = Molecule(smiles=smiles)
                print(' - Experimental difference: '
                      + '{: 10.1f} kcal/mol'.format(expv))
                print(' - Predicted difference: '
                      + '{: 10.1f} kcal/mol'.format(diff))
                print(' - Absolute error: '
                      + '{: 10.1f} kcal/mol'.format(abs(diff - expv)))
                display(mol)
