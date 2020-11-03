"""
This module contains classes and functions to validate the energetic
profiles of dihedrals.
"""


class EnergeticProfileBaseCalculator(object):
    """
    It represents the base class of any energetic profile calculator.
    """

    _name = ''

    def __init__(self, dihedral_benchmark):
        """
        It initializes an EnergeticProfileBaseCalculator object.

        Parameters
        ----------
        dihedral_benchmark : an peleffybenchmarktools.dihedrals.DihedralBenchmark object
            The DihedralBenchmark object that will be used to obtain the
            energetic profile
        """
        self._dihedral_benchmark = dihedral_benchmark

    def _get_plot_values(self, resolution):
        """
        It gets the required values for the plot.

        resolution : float
            The resolution, in degrees, that is applied when rotating
            the dihedral rotatable bond
        """
        raise NotImplementedError

    def get_energies(self):
        """
        It returns the energies of this energetic profile calculator.

        Returns
        -------
        energies : a simtk.unit.Quantity object
            The array of energies represented with a simtk's Quantity objec
        """
        raise NotImplementedError

    def plot_energies(self, resolution=30):
        """
        It plots the energies of this energetic profile calculator.

        Parameters
        ----------
        resolution : float
            The resolution, in degrees, that is applied when rotating
            the dihedral rotatable bond. Default is 30 degrees
        """
        from matplotlib import pyplot as plt

        x, y = self._get_plot_values(resolution)

        plt.title(self.name)
        plt.plot(x, y, 'r--', marker='x', markeredgecolor='k')
        plt.xlabel('Theta angle (degrees)')
        plt.ylabel('Dihedral energy (kcal/mol)')
        plt.xlim((0, 360))
        plt.show()

    @property
    def dihedral_benchmark(self):
        """
        A dihedral benchmark object that contains the system that will
        be validated.

        Returns
        -------
        dihedral_benchmark : an peleffybenchmarktools.dihedrals.DihedralBenchmark object
            The DihedralBenchmark object
        """
        return self._dihedral_benchmark

    @property
    def name(self):
        """
        The name of the current energetic profile calculator

        Returns
        -------
        name : str
            The energetic profile calculator name
        """
        return self._name


class OpenMMEnergeticProfile(EnergeticProfileBaseCalculator):
    """
    It represents an energetic profile calculator that employs OpenMM.
    """

    _name = 'OpenMM energetic profile'

    def __init__(self, dihedral_benchmark):
        """
        It initializes an OpenMMEnergeticProfile object.

        Parameters
        ----------
        dihedral_benchmark : an peleffybenchmarktools.dihedrals.DihedralBenchmark object
            The DihedralBenchmark object that will be used to obtain the
            energetic profile
        """
        super().__init__(dihedral_benchmark)

        from openforcefield.typing.engines.smirnoff import ForceField

        mol = self.dihedral_benchmark.molecule
        ff = ForceField(mol.forcefield + '.offxml')
        self._omm_top = mol.off_molecule.to_topology()
        self._omm_system = ff.create_openmm_system(self._omm_top)

    def _calc_energy(self, coords):
        """
        Given some coords, it calculates periodic torsion contribution
        to the energy, according to OpenMM.

        Parameters
        ----------
        coords : a simtk.unit.Quantity object
            The coordinates of a certain conformation of the molecule

        Returns
        -------
        omm_energy : a simtk.unit.Quantity
            The dihedral energy from OpenMM
        """
        from simtk.openmm import app, LangevinIntegrator, unit

        omm_idx_to_force = {}
        for idx, force in enumerate(self.omm_system.getForces()):
            force.setForceGroup(idx)
            omm_idx_to_force[str(type(force).__name__)] = idx

        omm_integrator = LangevinIntegrator(300 * unit.kelvin,
                                            1 / unit.picosecond,
                                            0.002 * unit.picoseconds)
        omm_simulation = app.Simulation(self.omm_top, self.omm_system,
                                        omm_integrator)
        omm_simulation.context.setPositions(coords)

        omm_energy = omm_simulation.context.getState(
            getEnergy=True,
            groups={omm_idx_to_force['PeriodicTorsionForce']}
        ).getPotentialEnergy()

        return omm_energy

    def get_energies(self, resolution=30, get_thetas=False):
        """
        It returns the energies of this energetic profile calculator.

        Parameters
        ----------
        resolution : float
            The resolution, in degrees, that is applied when rotating
            the dihedral rotatable bond. Default is 30 degrees
        get_thetas : bool
            Whether to return thetas (dihedral angles) or not. Default
            is False

        Returns
        -------
        dihedral_energies : a simtk.unit.Quantity object
            The array of energies represented with a simtk's Quantity
            object
        thetas : list[float]
            The array of thetas, only if requested in the corresponding
            function parameter
        """
        mol_with_conformers = \
            self.dihedral_benchmark.generate_dihedral_conformers(
                resolution)

        from openforcefield.topology import Molecule

        mol = Molecule.from_rdkit(mol_with_conformers)

        dihedral_energies = list()

        for i, conformer in enumerate(mol.conformers):
            energy = self._calc_energy(conformer)
            dihedral_energies.append(energy)

        if get_thetas:
            from rdkit.Chem import rdMolTransforms

            thetas = list()

            for conf in mol_with_conformers.GetConformers():
                thetas.append(rdMolTransforms.GetDihedralDeg(
                    conf, *self.dihedral_benchmark.atom_indexes))

            return dihedral_energies, thetas

        return dihedral_energies

    def _get_plot_values(self, resolution=30):
        """
        It plots the energies of this energetic profile calculator.

        Parameters
        ----------
        resolution : float
            The resolution, in degrees, that is applied when rotating
            the dihedral rotatable bond. Default is 30 degrees

        Returns
        -------
        thetas : list[float]
            The array of thetas
        dihedral_energies : list[float]
            The array of dihedral energies
        """
        from simtk import unit

        dihedral_energies, thetas = self.get_energies(resolution,
                                                      get_thetas=True)

        # Thetas reduction to first period
        reduced_thetas = list()

        for theta in thetas:
            if theta < 0:
                reduced_thetas.append(theta + 360)
            else:
                reduced_thetas.append(theta)

        # Pair thetas and energies
        theta_to_energy = dict(zip(reduced_thetas, dihedral_energies))

        # Get values in order
        ordered_thetas = sorted(theta_to_energy)
        ordered_dihedral_energies = [theta_to_energy[t].value_in_unit(
            unit.kilocalorie / unit.mole) for t in ordered_thetas]

        return ordered_thetas, ordered_dihedral_energies

    @property
    def omm_top(self):
        """
        The OpenMM topology corresponding to the molecule to benchmark.

        Returns
        -------
        omm_top : an openforcefield.topology.topology.Topology object
            The OpenMM topology
        """
        return self._omm_top

    @property
    def omm_system(self):
        """
        The OpenMM system to be used in the energy calculations.

        Returns
        -------
        omm_system : a simtk.openmm.openmm.System
            The OpenMM system
        """
        return self._omm_system


class OpenFFEnergeticProfile(EnergeticProfileBaseCalculator):
    """
    It represents an energetic profile calculator that employs the
    theoretical equations from the OpenFF Toolkit.
    """

    _name = 'OpenFF energetic profile'

    def __init__(self, dihedral_benchmark):
        """
        It initializes an OpenFFEnergeticProfile object.

        Parameters
        ----------
        dihedral_benchmark : an peleffybenchmarktools.dihedrals.DihedralBenchmark object
            The DihedralBenchmark object that will be used to obtain the
            energetic profile
        """
        super().__init__(dihedral_benchmark)

        from openforcefield.topology import Topology
        from openforcefield.typing.engines.smirnoff import ForceField

        mol = self.dihedral_benchmark.molecule
        topology = Topology.from_molecules([mol.off_molecule])
        ff = ForceField(mol.forcefield + '.offxml')
        parameters = ff.label_molecules(topology)[0]
        self._parameters = dict(parameters['ProperTorsions'])

    def get_energies(self, resolution=30, get_thetas=False):
        """
        It returns the energies of this energetic profile calculator.

        Parameters
        ----------
        resolution : float
            The resolution, in degrees, that is applied when rotating
            the dihedral rotatable bond. Default is 30 degrees
        get_thetas : bool
            Whether to return thetas (dihedral angles) or not. Default
            is False

        Returns
        -------
        dihedral_energies : a simtk.unit.Quantity object
            The array of energies represented with a simtk's Quantity
            object
        thetas : list[float]
            The array of thetas, only if requested in the corresponding
            function parameter
        """

        import numpy as np
        from copy import deepcopy
        from simtk import unit
        from rdkit.Chem import rdMolTransforms

        xs = unit.Quantity(np.arange(0, 360, resolution),
                           unit=unit.degrees)
        ys = unit.Quantity(np.zeros(len(xs)),
                           unit=(unit.kilocalorie / unit.mole))

        rdkit_mol = deepcopy(
            self.dihedral_benchmark.molecule.rdkit_molecule)
        conformer = rdkit_mol.GetConformer()

        for indexes in self.parameters:
            rot_bond = set(list(indexes)[1:3])
            a1, a2 = self.dihedral_benchmark.rotatable_bond

            if a1 in rot_bond and a2 in rot_bond:
                off_torsion = self.parameters[indexes]

                for i, k in enumerate(off_torsion.k):
                    if k == 0:
                        continue

                    phase = off_torsion.phase[i]
                    idivf = off_torsion.idivf[i]
                    periodicity = off_torsion.periodicity[i]

                    k /= idivf

                    for j, x in enumerate(xs):
                        rdMolTransforms.SetDihedralDeg(
                            conformer, *self.dihedral_benchmark.atom_indexes,
                            float(x.value_in_unit(unit.degree)))
                        theta = unit.Quantity(
                            value=rdMolTransforms.GetDihedralDeg(conformer,
                                                                 *indexes),
                            unit=unit.degree)

                        ys[j] += k * (1 + np.cos(periodicity
                                                 * theta.value_in_unit(
                                                     unit.radian)
                                                 - phase.value_in_unit(
                                                     unit.radian)))

        if get_thetas:
            return ys, xs

        return ys

    def _get_plot_values(self, resolution=30):
        """
        It plots the energies of this energetic profile calculator.

        Parameters
        ----------
        resolution : float
            The resolution, in degrees, that is applied when rotating
            the dihedral rotatable bond. Default is 30 degrees

        Returns
        -------
        thetas : list[float]
            The array of thetas
        dihedral_energies : list[float]
            The array of dihedral energies
        """
        dihedral_energies, thetas = self.get_energies(resolution,
                                                      get_thetas=True)

        return thetas, dihedral_energies

    @property
    def parameters(self):
        """
        The OpenFF force field to use to obtain the dihedral parameters.

        Returns
        -------
        parameters : dict
            The dictionary containing the dihedral parameters, keyed by
            atom indexes
        """
        return self._parameters


class PELEEnergeticProfile(EnergeticProfileBaseCalculator):
    """
    It represents an energetic profile calculator that employs PELE
    single point calculations.
    """

    _name = 'PELE energetic profile'

    def __init__(self, dihedral_benchmark, PELE_exec, PELE_src,
                 PELE_license):
        """
        It initializes an PELEEnergeticProfile object.

        Parameters
        ----------
        dihedral_benchmark : an peleffybenchmarktools.dihedrals.DihedralBenchmark object
            The DihedralBenchmark object that will be used to obtain the
            energetic profile
        PELE_exec : str
            Path to the PELE executable
        PELE_src : str
            Path to PELE source folder
        PELE_license : str
            Path to PELE license directory
        """
        super().__init__(dihedral_benchmark)

        self._PELE_exec = PELE_exec
        self._PELE_src = PELE_src
        self._PELE_license = PELE_license

    def get_energies(self, resolution=30, get_thetas=False):
        """
        It returns the energies of this energetic profile calculator.

        Parameters
        ----------
        resolution : float
            The resolution, in degrees, that is applied when rotating
            the dihedral rotatable bond. Default is 30 degrees
        get_thetas : bool
            Whether to return thetas (dihedral angles) or not. Default
            is False

        Returns
        -------
        dihedral_energies : a simtk.unit.Quantity object
            The array of energies represented with a simtk's Quantity
            object
        thetas : list[float]
            The array of thetas, only if requested in the corresponding
            function parameter
        """

        mol_with_conformers = \
            self.dihedral_benchmark.generate_dihedral_conformers(
                resolution)

        import tempfile
        from peleffybenchmarktools.utils import temporary_cd
        from peleffybenchmarktools.utils.pele import (PELESinglePoint,
                                                      PELEOutputParser)

        mol = self.dihedral_benchmark.molecule
        dihedral_energies = list()

        with tempfile.TemporaryDirectory() as tmpdir:
            with temporary_cd(tmpdir):
                for conf in mol_with_conformers.GetConformers():
                    # Write conformation to PDB file
                    mol.set_conformer(conf)
                    mol.to_pdb_file('ligand.pdb')

                    pele_job = PELESinglePoint(PELE_exec=self.PELE_exec,
                                               PELE_src=self.PELE_src,
                                               PELE_license=self.PELE_license)

                    output_file = pele_job.run(
                        self.dihedral_benchmark.molecule,
                        pdb_path='ligand.pdb',
                        charges_method='gasteiger')

                    parsed_output = PELEOutputParser(output_file)
                    dihedral_energies.append(parsed_output['torsion energy'])

        if get_thetas:
            from rdkit.Chem import rdMolTransforms

            thetas = list()

            for conf in mol_with_conformers.GetConformers():
                thetas.append(rdMolTransforms.GetDihedralDeg(
                    conf, *self.dihedral_benchmark.atom_indexes))

            return dihedral_energies, thetas

        return dihedral_energies

    def _get_plot_values(self, resolution=30):
        """
        It plots the energies of this energetic profile calculator.

        Parameters
        ----------
        resolution : float
            The resolution, in degrees, that is applied when rotating
            the dihedral rotatable bond. Default is 30 degrees

        Returns
        -------
        thetas : list[float]
            The array of thetas
        dihedral_energies : list[float]
            The array of dihedral energies
        """
        from simtk import unit

        dihedral_energies, thetas = self.get_energies(resolution,
                                                      get_thetas=True)

        # Thetas reduction to first period
        reduced_thetas = list()

        for theta in thetas:
            if theta < 0:
                reduced_thetas.append(theta + 360)
            else:
                reduced_thetas.append(theta)

        # Pair thetas and energies
        theta_to_energy = dict(zip(reduced_thetas, dihedral_energies))

        # Get values in order
        ordered_thetas = sorted(theta_to_energy)
        ordered_dihedral_energies = [theta_to_energy[t].value_in_unit(
            unit.kilocalorie / unit.mole) for t in ordered_thetas]

        return ordered_thetas, ordered_dihedral_energies

    @property
    def PELE_exec(self):
        """
        Path to the PELE executable

        Returns
        -------
        PELE_exec : str
            The path to the PELE executable
        """
        return self._PELE_exec

    @property
    def PELE_src(self):
        """
        Path to the PELE source code

        Returns
        -------
        PELE_src : str
            The path to the PELE source code
        """
        return self._PELE_src

    @property
    def PELE_license(self):
        """
        Path to PELE license directory

        Returns
        -------
        PELE_license : str
            The path to the PELE license directory
        """
        return self._PELE_license

    @property
    def forcefield(self):
        """
        The OpenFF force field to employ to parameterize the molecule.

        Returns
        -------
        forcefield : str
            The forcefield name
        """
        return self._forcefield


class OFFPELEEnergeticProfile(EnergeticProfileBaseCalculator):
    """
    It represents an energetic profile calculator that employs the
    theoretical equations from OFFPELE.
    """

    _name = 'OFF-PELE energetic profile'

    def __init__(self, dihedral_benchmark):
        """
        It initializes an OpenFFEnergeticProfile object.

        Parameters
        ----------
        dihedral_benchmark : an peleffybenchmarktools.dihedrals.DihedralBenchmark object
            The DihedralBenchmark object that will be used to obtain the
            energetic profile
        """
        super().__init__(dihedral_benchmark)

    def get_energies(self, resolution=30, get_thetas=False):
        """
        It returns the energies of this energetic profile calculator.

        Parameters
        ----------
        resolution : float
            The resolution, in degrees, that is applied when rotating
            the dihedral rotatable bond. Default is 30 degrees
        get_thetas : bool
            Whether to return thetas (dihedral angles) or not. Default
            is False

        Returns
        -------
        dihedral_energies : a simtk.unit.Quantity object
            The array of energies represented with a simtk's Quantity
            object
        thetas : list[float]
            The array of thetas, only if requested in the corresponding
            function parameter
        """

        import numpy as np
        from copy import deepcopy
        from simtk import unit
        from rdkit.Chem import rdMolTransforms

        xs = unit.Quantity(np.arange(0, 360, resolution),
                           unit=unit.degrees)
        ys = unit.Quantity(np.zeros(len(xs)),
                           unit=(unit.kilocalorie / unit.mole))

        rdkit_mol = deepcopy(
            self.dihedral_benchmark.molecule.rdkit_molecule)
        conformer = rdkit_mol.GetConformer()

        propers = self.dihedral_benchmark.molecule.propers

        a1, a2 = self.dihedral_benchmark.rotatable_bond

        for proper in propers:
            rot_bond = set((proper.atom2_idx, proper.atom3_idx))

            if a1 in rot_bond and a2 in rot_bond:
                if proper.constant == 0:
                    continue

                constant = proper.constant
                prefactor = proper.prefactor
                periodicity = proper.periodicity
                phase = proper.phase

                for j, x in enumerate(xs):
                    rdMolTransforms.SetDihedralDeg(
                        conformer, *self.dihedral_benchmark.atom_indexes,
                        float(x.value_in_unit(unit.degree)))
                    theta = unit.Quantity(
                        value=rdMolTransforms.GetDihedralDeg(conformer,
                                                             *indexes),
                        unit=unit.degree)

                    ys[j] += constant * (1 + prefactor * np.cos(
                        periodicity * theta.value_in_unit(unit.radian)
                        - phase.value_in_unit(unit.radian)))

        if get_thetas:
            return ys, xs

        return ys

    def _get_plot_values(self, resolution=30):
        """
        It plots the energies of this energetic profile calculator.

        Parameters
        ----------
        resolution : float
            The resolution, in degrees, that is applied when rotating
            the dihedral rotatable bond. Default is 30 degrees

        Returns
        -------
        thetas : list[float]
            The array of thetas
        dihedral_energies : list[float]
            The array of dihedral energies
        """
        dihedral_energies, thetas = self.get_energies(resolution,
                                                      get_thetas=True)

        return thetas, dihedral_energies
