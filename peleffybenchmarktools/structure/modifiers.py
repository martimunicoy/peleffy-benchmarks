"""
This module contains different methods to modify the structure of a
molecule.
"""


def _DistortionWrapper(object):
    """
    A wrapper for any structural distortion.
    """

    def __init__(self, molecule):
        """
        It initializes a distortion wrapper object.

        Parameters
        ----------
        molecule : an peleffy.topology.Molecule
            The peleffy's Molecule object. It needs to be previously
            parameterized
        """
        molecule.assert_parameterized()

        self._original_mol = molecule


def DistortBonds(_DistortionWrapper):
    """
    It defines a class that distorts a set of bonds. Note that bonds
    belonging to a ring will not be modified.
    """

    def randomly(self, range):
        """
        Given a distortion range, it applies a random distortion to bond
        lengths.

        Parameters
        ----------
        range : float
            The distortion range in angstroms to apply to the molecule

        Returns
        -------
        distorted_mol : an RDKit.molecule object
            The resulting distorted molecule
        """
        from copy import deepcopy
        from rdkit.Chem import rdMolTransforms
        from random import uniform

        rdkit_mol = deepcopy(self._original_mol.rdkit_molecule)
        conformer = rdkit_mol.GetConformer()

        for bond in self._original_mol.bonds:
            idx1 = bond.atom1_idx
            idx2 = bond.atom2_idx

            if rdkit_mol.GetBondBetweenAtoms(idx1, idx2).IsInRing():
                continue

            bond_length = rdMolTransforms.GetBondLength(conformer,
                                                        idx1, idx2)
            bond_delta = uniform(-range, range)

            rdMolTransforms.SetBondLength(conformer, idx1, idx2,
                                          bond_length + bond_delta)

        rdkit_mol.AddConformer(conformer, assignId=True)
        rdkit_mol.RemoveConformer(conformer.GetId())

        return rdkit_mol

    def fixed(self, value):
        """
        Given a distortion range, it applies a random distortion to bond
        lengths.

        Parameters
        ----------
        value : float
            The distortion value in angstroms to apply to the molecule

        Returns
        -------
        distorted_mol : an RDKit.molecule object
            The resulting distorted molecule
        """
        from copy import deepcopy
        from rdkit.Chem import rdMolTransforms

        rdkit_mol = deepcopy(self._original_mol.rdkit_molecule)
        conformer = rdkit_mol.GetConformer()

        for bond in self._original_mol.bonds:
            idx1 = bond.atom1_idx
            idx2 = bond.atom2_idx

            if rdkit_mol.GetBondBetweenAtoms(idx1, idx2).IsInRing():
                continue

            bond_length = rdMolTransforms.GetBondLength(conformer,
                                                        idx1, idx2)
            bond_delta = value

            rdMolTransforms.SetBondLength(conformer, idx1, idx2,
                                          bond_length + bond_delta)

        rdkit_mol.AddConformer(conformer, assignId=True)
        rdkit_mol.RemoveConformer(conformer.GetId())

        return rdkit_mol


def DistortAngles(_DistortionWrapper):
    """
    It defines a class that distorts a set of angles. Note that angles
    belonging to a ring will not be modified.
    """

    def randomly(self, range):
        """
        Given a distortion range, it applies a random distortion to
        torsion angles.

        Parameters
        ----------
        range : float
            The distortion range in degrees to apply to the molecule

        Returns
        -------
        distorted_mol : an RDKit.molecule object
            The resulting distorted molecule
        """
        from copy import deepcopy
        from rdkit.Chem import rdMolTransforms
        from random import uniform

        rdkit_mol = deepcopy(self._original_mol.rdkit_molecule)
        conformer = rdkit_mol.GetConformer()

        for angle in self._original_mol.angles:
            idx1 = angle.atom1_idx
            idx2 = angle.atom2_idx
            idx3 = angle.atom3_idx

            if rdkit_mol.GetBondBetweenAtoms(idx1, idx2).IsInRing():
                continue

            if rdkit_mol.GetBondBetweenAtoms(idx2, idx3).IsInRing():
                continue

            torsion_angle = rdMolTransforms.GetBondLength(conformer,
                                                          idx1, idx2, idx3)
            torsion_delta = uniform(-range, range)

            rdMolTransforms.SetAngleDeg(conformer, idx1, idx2, idx3,
                                        torsion_angle + torsion_delta)

        rdkit_mol.AddConformer(conformer, assignId=True)
        rdkit_mol.RemoveConformer(conformer.GetId())

        return rdkit_mol

    def fixed(self, value):
        """
        Given a distortion range, it applies a random distortion to
        torsion angles.

        Parameters
        ----------
        value : float
            The distortion value in degrees to apply to the molecule

        Returns
        -------
        distorted_mol : an RDKit.molecule object
            The resulting distorted molecule
        """
        from copy import deepcopy
        from rdkit.Chem import rdMolTransforms

        rdkit_mol = deepcopy(self._original_mol.rdkit_molecule)
        conformer = rdkit_mol.GetConformer()

        for angle in self._original_mol.angles:
            idx1 = angle.atom1_idx
            idx2 = angle.atom2_idx
            idx3 = angle.atom3_idx

            if rdkit_mol.GetBondBetweenAtoms(idx1, idx2).IsInRing():
                continue

            if rdkit_mol.GetBondBetweenAtoms(idx2, idx3).IsInRing():
                continue

            torsion_angle = rdMolTransforms.GetBondLength(conformer,
                                                          idx1, idx2, idx3)
            torsion_delta = value

            rdMolTransforms.SetAngleDeg(conformer, idx1, idx2, idx3,
                                        torsion_angle + torsion_delta)

        rdkit_mol.AddConformer(conformer, assignId=True)
        rdkit_mol.RemoveConformer(conformer.GetId())

        return rdkit_mol


def DistortDihedrals(_DistortionWrapper):
    """
    It defines a class that distorts a set of dihedrals. Note that
    dihedrals belonging to a ring will not be modified.
    """

    def randomly(self, range):
        """
        Given a distortion range, it applies a random distortion to
        dihedral angles.

        Parameters
        ----------
        range : float
            The distortion range in degrees to apply to the molecule

        Returns
        -------
        distorted_mol : an RDKit.molecule object
            The resulting distorted molecule
        """
        from copy import deepcopy
        from rdkit.Chem import rdMolTransforms
        from random import uniform

        rdkit_mol = deepcopy(self._original_mol.rdkit_molecule)
        conformer = rdkit_mol.GetConformer()

        for proper in self._original_mol.propers:
            idx1 = proper.atom1_idx
            idx2 = proper.atom2_idx
            idx3 = proper.atom3_idx
            idx4 = proper.atom4_idx

            if rdkit_mol.GetBondBetweenAtoms(idx2, idx3).IsInRing():
                continue

            proper_angle = rdMolTransforms.GetDihedralDeg(conformer,
                                                          idx1, idx2,
                                                          idx3, idx4)

            proper_delta = uniform(-range, range)

            rdMolTransforms.SetDihedralDeg(conformer,
                                           idx1, idx2, idx3, idx4,
                                           proper_angle + proper_delta)

        rdkit_mol.AddConformer(conformer, assignId=True)
        rdkit_mol.RemoveConformer(conformer.GetId())

        return rdkit_mol

    def fixed(self, value):
        """
        Given a distortion range, it applies a random distortion to
        dihedral angles.

        Parameters
        ----------
        value : float
            The distortion value in degrees to apply to the molecule

        Returns
        -------
        distorted_mol : an RDKit.molecule object
            The resulting distorted molecule
        """
        from copy import deepcopy
        from rdkit.Chem import rdMolTransforms

        rdkit_mol = deepcopy(self._original_mol.rdkit_molecule)
        conformer = rdkit_mol.GetConformer()

        for proper in self._original_mol.propers:
            idx1 = proper.atom1_idx
            idx2 = proper.atom2_idx
            idx3 = proper.atom3_idx
            idx4 = proper.atom4_idx

            if rdkit_mol.GetBondBetweenAtoms(idx2, idx3).IsInRing():
                continue

            proper_angle = rdMolTransforms.GetDihedralDeg(conformer,
                                                          idx1, idx2,
                                                          idx3, idx4)

            proper_delta = value

            rdMolTransforms.SetDihedralDeg(conformer,
                                           idx1, idx2, idx3, idx4,
                                           proper_angle + proper_delta)

        rdkit_mol.AddConformer(conformer, assignId=True)
        rdkit_mol.RemoveConformer(conformer.GetId())

        return rdkit_mol
