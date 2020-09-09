

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

    def minimize(self, input_PDB_file, solvent='vacuum',
                 forcefield='openff_unconstrained-1.2.0.offxml',
                 charges_method='am1bcc', output_path=None):
        """
        Given an input PDB file, it runs a minimization with PELE.

        Parameters
        ----------
        input_PDB_file : str
            The input PDB file to minimize with PELE
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

        self._set_output_path(input_PDB_file, output_path)

        self._create_directory(input_PDB_file)

        self._link_folders()

        self._generate_parameters(input_PDB_file, forcefield=forcefield,
                                  charges_method=charges_method)

        self._run_PELE_minimization(solvent)

    def _set_output_path(self, input_PDB_file, output_path=None):
        """
        It handles and sets output path

        Parameters
        ----------
        input_PDB_file : str
            The input PDB file to minimize with PELE
        output_path : str
            The output path where results will be saved
        """
        import os
        from pathlib import Path

        if output_path is None:
            self._output_path = os.path.join('output',
                                             Path(input_PDB_file).stem)

    def _create_directory(self, input_PDB_file):
        """
        It creates an output directory where all the results will be saved
        and copies the ligand's pdb inside.

        Parameters
        ----------
        input_PDB_file : str
            The input PDB file to minimize with PELE
        """
        import shutil
        import os

        # It makes the output directory
        os.makedirs(self._output_folder, exist_ok=True)
        shutil.copy(input_PDB_file, os.path.join(os.getcwd(),
                                                 self._output_folder,
                                                 'ligand.pdb'))

    def _link_folders(self):
        """
        It links the necessary folders to the output folder.
        """
        import os

        # Link to Data
        os.symlink(os.path.join(self._PELE_src, 'Data'),
                   os.path.join(os.getcwd(), self._output_path, 'Data'))

        # Link to Documents
        os.symlink(os.path.join(self._PELE_src, 'Documents'),
                   os.path.join(os.getcwd(), self._output_path, 'Documents'))

    def _generate_parameters(self,
                             input_PDB_file,
                             forcefield='openff_unconstrained-1.2.0.offxml',
                             charges_method='am1bcc',
                             path=None):
        """
        It generates the parameters of the molecule (from the input_file)
        as DataLocal in the output folder.

        Parameters
        ----------
        input_PDB_file : str
            The input PDB file to parameterize with the Open Force Field
            Toolkit for PELE
        forcefield : str
            The Open Force Field force field to generate the parameters
            with
        charges_method : str
            The charges method to calculate the partial charges with
        path : str
            Path where parameters will be saved
        """
        import offpele
        from offpele.topology import Molecule
        from offpele.template import Impact
        from offpele.solvent import OBC2
        from offpele.main import handle_output_paths
        from pathlib import Path

        if path is None:
            general_path = Path(input_PDB_file).parent

        # Create representation of a particular molecule
        molecule = Molecule(input_PDB_file)

        # Saving paths
        rotamer_library_output_path, impact_output_path, \
            solvent_output_path = handle_output_paths(molecule=molecule,
                                                      output=general_path,
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

    def _run_PELE_minimization(self, solvent):
        """
        It runs a PELE minimization.

        Parameters
        ----------
        solvent : str
            The solvent name. One of ['vacuum', 'OBC']. Default is 'vacuum'
        """
        import os

        # Minimization
        os.chdir(os.path.join(os.getcwd(), self._output_path))
        os.system("{} {} > {}_minimization.out".format(
            self._PELE_exec, self.CONTROL_FILES[solvent], solvent))
