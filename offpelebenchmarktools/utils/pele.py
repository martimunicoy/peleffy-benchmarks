"""
This module contains classes and functions to facilitate the execution
of PELE workflows.
"""


class PELEBaseJob(object):
    """
    It represents the base class of a PELE job.
    """

    _name = ''
    _CONTROL_FILES = dict()

    def __init__(self, PELE_exec, PELE_src, output_path=None):
        """
        It initializes a PELEBaseJob object.

        Parameters
        ----------
        PELE_exec : str
            Path to the PELE executable
        PELE_src : str
            Path to PELE source folder
        output_path : str
            The path to save the output coming from PELE
        """
        self._PELE_exec = PELE_exec
        self._PELE_src = PELE_src
        self._output_path = None

    def set_output_path(self, output_path):
        """
        It sets the path to save the output coming from PELE.

        Parameters
        ----------
        output_path : str
            The output path
        """
        self._output_path = output_path

    def _select_control_file(self):
        """
        It selects the PELE's control file to run.

        Returns
        -------
        control_file : str
            The path to the control file
        """
        raise NotImplementedError

    def _get_output_path(self, molecule):
        """
        It sets the output path.

        Parameters
        ----------
        molecule : an offpele.topology.Molecule object
            The molecule to run the PELE workflow on

        Returns
        -------
        output_path : str
            The output path that has been assigned
        """
        import os

        if self.output_path is None:
            output_path = os.path.join('output', molecule.name)
        else:
            output_path = os.path.join(output_path, molecule.name)

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

    def _generate_parameters(self, molecule, output_path,
                             forcefield='openff_unconstrained-1.2.0.offxml',
                             charges_method='am1bcc'):
        """
        It generates the parameters of the molecule (from the input_file)
        as DataLocal in the output folder.

        Parameters
        ----------
        molecule : an offpele.topology.Molecule object
            The molecule to run the PELE workflow on
        output_path : str
            The output path where results will be saved
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
        import os

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

    def run(self, molecule, pdb_path=None,
            forcefield='openff_unconstrained-1.2.0.offxml',
            charges_method='am1bcc'):
        """
        It runs a job with PELE.

        Parameters
        ----------
        molecule : an offpele.topology.Molecule object
            The molecule to run the PELE workflow on
        pdb_path : str
            The path to the PDB file to use as the input structure for
            the PELE job. Default is None
        forcefield : str
            The Open Force Field force field to generate the parameters
            with. Default is 'openff_unconstrained-1.2.0.offxml'
        charges_method : str
            The charges method to calculate the partial charges with.
            Default is 'am1bcc'
        """
        import os

        control_file = self._select_control_file()

        output_path = self._get_output_path(molecule)

        self._create_directory(output_path)

        self._link_folders(output_path)

        self._generate_parameters(molecule, output_path,
                                  forcefield=forcefield,
                                  charges_method=charges_method)

        if pdb_path is None:
            molecule.to_pdb_file(os.path.join(output_path, 'ligand.pdb'))
        else:
            from shutil import copyfile
            copyfile(pdb_path, os.path.join(output_path, 'ligand.pdb'))

        previous_dir = os.getcwd()
        os.chdir(os.path.join(os.getcwd(), output_path))
        os.system("{} {} > {}PELE_output.txt".format(
            self._PELE_exec, control_file))
        os.chdir(previous_dir)

    @property
    def name(self):
        """
        The name of this PELEBaseJob object.

        Returns
        -------
        name : str
            The corresponding name
        """
        return self._name

    @property
    def output_path(self):
        """
        The path to save the output coming from PELE.

        output_path : str
            The output path
        """
        return self._output_path


class PELESinglePoint(PELEBaseJob):
    """
    It represents the single point class of a PELE job which runs PELE
    to compute the energy of the current conformation of a molecule.
    """

    def __init__(self, PELE_exec, PELE_src, output_path=None):
        """
        It initializes a PELESinglePoint job.

        Parameters
        ----------
        PELE_exec : str
            Path to the PELE executable
        PELE_src : str
            Path to PELE source folder
        output_path : str
            The path to save the output coming from PELE
        """
        super().__init__(PELE_exec, PELE_src, output_path)

        from offpelebenchmarktools.utils import get_data_file_path

        single_point_cf = get_data_file_path('PELE_control/single_point.conf')

        self._add_control_file('single_point', single_point_cf)

    def _select_control_file(self):
        """
        It selects the PELE's control file to run.

        Returns
        -------
        control_file : str
            The path to the control file
        """
        return self._CONTROL_FILES['single_point']
