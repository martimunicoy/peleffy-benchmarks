"""
This module contains the test to check the rotamers library.
"""

class TestRotamers(mae_file, pdb_file):

	def test_bonds(self):
		bonds_peleffy = get_bonds(pdb_file)
		bonds_PlopRotTemp = get_bonds(mae_file)
		diff1, diff2 = differentiate(bonds_peleffy,bonds_PlopRotTemp)
		write_diff(diff1,diff2)

#	def test_resolution(self):
#		return

#	def test_groups(self):
#		return
