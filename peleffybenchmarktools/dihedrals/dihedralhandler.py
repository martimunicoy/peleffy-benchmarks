"""
This module contains classes and functions to handle dihedral benchmarks.
"""


class DihedralBenchmark(object):
    """
    It defines a certain dihedral and allows to run several benchmarks
    on it.
    """

    def __init__(self, dihedral_atom_indexes, molecule):
        """
        It initializes a DihedralBenchmark object.

        Parameters
        ----------
        dihedral_atom_indexes : tuple[int]
            The indexes of atoms involved in the dihedral
        molecule : an peleffy.topology.Molecule
            The peleffy's Molecule object
        """

        # Hide peleffy output
        from peleffy.utils import Logger
        logger = Logger()
        logger.set_level('WARNING')

        molecule.assert_parameterized()

        self._atom_indexes = dihedral_atom_indexes
        self._molecule = molecule
        self._forcefield = molecule.forcefield

    def generate_dihedral_conformers(self, resolution):
        """
        Given a resolution, it generates all the possible conformers
        that come out when rotating the dihedral's rotatable bond.

        Parameters
        ----------
        resolution : float
            The resolution, in degrees, that is applied when rotating
            the dihedral rotatable bond. Default is 30 degrees

        Returns
        -------
        rdkit_mol : an rdkit.Chem.rdchem.Mol object
            The RDKit's Molecule object with the conformers that were
            generated
        """
        from copy import deepcopy
        from rdkit.Chem import rdMolTransforms

        rdkit_mol = deepcopy(self.molecule.rdkit_molecule)

        conformer = rdkit_mol.GetConformer()

        for angle in range(0, 360, resolution):
            rdMolTransforms.SetDihedralDeg(conformer, *self.atom_indexes,
                                           angle)
            rdkit_mol.AddConformer(conformer, assignId=True)

        # Remove initial conformer (which is repeated)
        rdkit_mol.RemoveConformer(conformer.GetId())

        return rdkit_mol

    def save_dihedral_conformers(self, resolution=30):
        """
        It saves the dihedral conformers that are obtained when
        rotating the dihedral's rotatable bond with a certain resolution.

        Parameters
        ----------
        resolution : float
            The resolution, in degrees, that is applied when rotating
            the dihedral rotatable bond. Default is 30 degrees
        """
        mol_with_conformers = self.generate_dihedral_conformers(resolution)

        from rdkit import Chem

        # Write structures
        w = Chem.SDWriter('conformers.sdf')
        for i, conf in enumerate(mol_with_conformers.GetConformers()):
            tm = Chem.Mol(mol_with_conformers, False, conf.GetId())
            w.write(tm)
        w.close()

    def display_dihedral(self):
        """
        It displays the dihedral that is being tracked
        """
        from rdkit import Chem
        from copy import deepcopy
        try:
            from IPython.display import display
            from IPython.display import SVG
        except ImportError:
            raise Exception('Jupyter notebook not found')

        rdk_mol = deepcopy(self.molecule.rdkit_molecule)

        Chem.rdDepictor.Compute2DCoords(rdk_mol)

        d = Chem.Draw.rdMolDraw2D.MolDraw2DSVG(250, 250)
        Chem.Draw.rdMolDraw2D.PrepareAndDrawMolecule(
            d, rdk_mol, highlightAtoms=self.atom_indexes)
        d.FinishDrawing()
        display(SVG(d.GetDrawingText()))

    def get_dihedral_parameters(self):
        """
        It returns the parameters of the dihedral that is being tracked,
        according to the OpenFF toolkit.

        Returns
        -------
        parameters : dict
            The dictionary with the parameters
        """
        from openforcefield.topology import Topology
        from openforcefield.typing.engines.smirnoff import ForceField

        topology = Topology.from_molecules([self.molecule.off_molecule])
        ff = ForceField(self.forcefield + '.offxml')
        all_parameters = ff.label_molecules(topology)[0]
        parameters = dict(all_parameters['ProperTorsions'])[self.atom_indexes]

        return parameters

    @property
    def atom_indexes(self):
        """
        The indexes of atoms involved in the dihedral.

        Returns
        -------
        atom_indexes : tuple[int]
            The tuple of atom indexes
        """
        return self._atom_indexes

    @property
    def molecule(self):
        """
        The molecule the dihedral belongs to.

        Returns
        -------
        molecule : an peleffy.topology.Molecule
            The peleffy's Molecule object
        """
        return self._molecule

    @property
    def forcefield(self):
        """
        The forcefield employed to parameterize the molecule.

        Returns
        -------
        forcefield : an openforcefield.typing.engines.smirnoff.ForceField
                     object
            The forcefield employed to parameterize this Molecule object
        """
        return self._forcefield

    @property
    def rotatable_bond(self):
        """
        It returns the atom indexes in the dihedral that belong to the
        rotatable bond.

        Returns
        -------
        rot_bond : tuple[int]
        The tuple of atom indexes
        """
        return tuple(list(self.atom_indexes)[1:3])
