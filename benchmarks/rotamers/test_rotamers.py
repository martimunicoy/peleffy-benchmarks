"""
This module contains the test to check the rotamers library.
"""
import utils as u


class TestRotamers:

    def __init__(self, mae_file, pdb_file, file):
        self.mae_file = mae_file
        self.pdb_file = pdb_file
        self.file = file

    def test_bonds(self):
        bonds_offpele = u.get_bonds(self.pdb_file)
        bonds_PlopRotTemp = u.get_bonds(self.mae_file)
        diff1, diff2 = u.differentiate(bonds_offpele, bonds_PlopRotTemp)
        u.write_diff_bonds(self.file, diff1, diff2)

    def test_resolution(self):
        resolution_offpele, resolution_PlopRotTemp = \
            u.get_resolution_elements(self.pdb_file, self.mae_file)
        resolution_differences = u.compare_resolution(resolution_offpele,
                                                      resolution_PlopRotTemp)
        u.write_diff_resolution(self.file, resolution_differences)
