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
                 charges_method='am1bcc', output_path=None):
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
        charges_method : str
            The charges method to calculate the partial charges with
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
                                  charges_method=charges_method)

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
        import shutil
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
                             charges_method='am1bcc'):
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
        charges_method : str
            The charges method to calculate the partial charges with
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
        molecule.parameterize(forcefield, charges_method=charges_method)
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
