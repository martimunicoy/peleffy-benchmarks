

class MoleculeMinimized(object):

	def __init__(self, path = none):
		self.input_file = input_file
		self.PELE_version = PELE_version

	def _output_folder(self,input_file):
		os.system("mkdir output")
		os.system("cp %s ./output/" % (input_file))
		os.system("cd output")
		os.system("mv ./output/%s ./output/ligand.pdb" % (input_file))
		os.system("cd output")


	def _generate_parameters():
		os.system("python /home/lauramalo/repos/offpele/offpele/main.py ./output/ligand.pdb --with_solvent --as_DataLocal")
		os.system("mv ./DataLocal/ ./output/DataLocal/")

	def _link_folders():
		os.system("ln -s /home/municoy/repos/PELE-repo/Data ./output/")
		os.system("ln -s /home/municoy/repos/PELE-repo/Documents ./output/") 

	def minimize(input_file, PELE_version):
		_output_folder(input_file)
    	_generate_parameters()
    	_link_folders()
		os.system(" %s /home/lauramalo/tests/hydrationfreeenergies/VACUUM_minimization.conf > ./output/VACUUM_minimization.out" % (PELE_version))
		os.system(" %s /home/lauramalo/tests/hydrationfreeenergies/OBC_minimization.conf > ./output/OBC_minimization.out" % (PELE_version))

