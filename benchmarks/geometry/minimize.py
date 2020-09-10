"""
It prepares the system files and runs a PELE minimization.
"""


class Minimizer(object):
    """
    It contains all the tools to minimize a molecule with PELE and the
    OpenForceField toolkit for PELE.
    """

    CONTROL_FILES = {'vacuum': 'data/VACUUM_minimization.conf',
                     'OBC': 'data/OBC_minimization.conf'}

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
        import offpele
        from offpele.topology import Molecule
        from offpele.template import Impact
        from offpele.solvent import OBC2
        from offpele.main import handle_output_paths

        # Create representation of a particular molecule
        molecule = Molecule(smiles=smiles, name=mol_id)

        # Saving paths
        rotamer_library_output_path, impact_output_path, \
            solvent_output_path = handle_output_paths(molecule=molecule,
                                                      output=output_path,
                                                      as_datalocal=True)

        # Generate its rotamer library
        rotamer_library = offpele.topology.RotamerLibrary(molecule)
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
