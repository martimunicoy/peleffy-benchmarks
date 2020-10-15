from offpelebenchmarktools.utils import get_data_file_path


FREESOLV_PATH = 'databases/FreeSolv0.52.txt'


class SolventBenchmark(object):
    def __init__(self, method,
                 pele_exec, pele_src, pele_license,
                 ploprottemp_src, schrodinger_src,
                 off_forcefield='openff_unconstrained-1.2.0.offxml',
                 charges_method='am1bcc', solvent='OBC',
                 opls_nonbonding=False, opls_bonds_angles=False):
        """
        It initialized an SolventBenchmark object.

        Parameters
        ----------
        method : an offpelebenchmarktools.solvent.method object
            The Method object that will be used to run the benchmark
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
        off_forcefield : str
            The OpenFF force field
        charges_method : str
            The method to calculate partial charges
        solvent : str
            The solvent model to employ
        opls_nonbonding : bool
            Whether to use OPLS2005 to parameterize nonbonding terms or not
        opls_bonds_angles : bool
            Whether to use OPLS2005 to paramterize bonds and angles or not
        """
        self.pele_exec = pele_exec
        self.pele_src = pele_src
        self.pele_license = pele_license
        self.method = method
        self.off_forcefield = off_forcefield
        self.charges_method = charges_method
        self.solvent = solvent
        self.opls_nonbonding = opls_nonbonding
        self.opls_bonds_angles = opls_bonds_angles
        self.ploprottemp_src = ploprottemp_src
        self.schrodinger_src = schrodinger_src

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
        from methods import method_OFF, method_OPLS, method_OFFOPLS

        compound_ids, smiles_tags, experimental_v = self._read_dataset()

        # It runs the selected method
        if self.method == 'OFF':
            energies, differences, experimental_values = \
                method_OFF(out_folder, compound_ids, smiles_tags,
                           experimental_v, self.off_forcefield,
                           self.charges_method, self.pele_exec,
                           self.pele_src, self.pele_license)

        elif self.method == 'OPLS':

            VACUUM_CF = '/home/lauramalo/repos/offpele-benchmarks/benchmarks/solvent/Conf/OPLS_VACUUM_minimization.conf'
            OBC_CF = '/home/lauramalo/repos/offpele-benchmarks/benchmarks/solvent/Conf/OPLS_OBC_minimization.conf'
            VDGBNP_CF = '/home/lauramalo/repos/offpele-benchmarks/benchmarks/solvent/Conf/OPLS_VDGBNP_minimization.conf'

            energies, differences, experimental_values = \
                method_OPLS(out_folder, compound_ids, smiles_tags,
                            experimental_v, self.solvent, self.pele_exec,
                            self.schrodinger_src, self.ploprottemp_src,
                            self.pele_src)

        elif self.method == 'OFFOPLS':

            VACUUM_CF = '/home/lauramalo/repos/offpele-benchmarks/benchmarks/solvent/Conf/OPLS_VACUUM_minimization.conf'
            OBC_CF = '/home/lauramalo/repos/offpele-benchmarks/benchmarks/solvent/Conf/OPLS_VDGBNP_minimization.conf'

            energies, differences, experimental_values = \
                method_OFFOPLS(out_folder, compound_ids, smiles_tags,
                               experimental_v, self.solvent,
                               self.opls_nonbonding,
                               self.opls_bonds_angles,
                               self.off_openforcefield,
                               self.pele_src)
        return energies, differences, experimental_values

    def save_output(self, out_folder, energies):
        """It saves the output results."""
        import pandas as pd
        import os

        df = pd.DataFrame(energies, columns=['CID', 'Energetic Difference', 'Experimental value'])
        df.to_csv(os.path.join(out_folder, 'results.txt'))

    def plot_results(self, energies=None, differences=None, experimental_values=None, file=None):
        """
        It generates and histogram and a regression for the comparision
        between the experimental values and the computed for the hydration
        free energy.
        """
        import matplotlib.pyplot as plt
        import numpy as np
        import pandas as pd

        if file is not None:
            data = pd.read_csv(file, usecols=['CID',
                                              'Energetic Difference',
                                              'Experimental value'])
            y = data['Energetic Difference'].to_numpy()
            x = data['Experimental value'].to_numpy()
        else:
            y = np.array(differences)
            x = np.array(experimental_values)

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
