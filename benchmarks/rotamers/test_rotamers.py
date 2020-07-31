"""
This module contains the test to check the rotamers library.
"""

class TestRotamers(mae_file, pdb_file):
	def test_groups(self):
		return

	def test_bonds(self):
		bonds_offpele = get_bonds(pdb_file)
		bonds_PlopRotTemp = get_bonds(mae_file)
		diff1, diff2 = differentiate(bonds_offpele,bonds_PlopRotTemp)


	def test_resolution(self):
		return
