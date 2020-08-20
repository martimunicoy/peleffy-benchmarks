import os 

class MoleculeMinimized:
	"""
	It contains all the tools to minimize a molecule with the OpenForceField toolkit for PELE.
	"""

	def __init__(self, input_file, PELE_version):
		"""
		It initializes a MocelueMinimized object. 

		Parameters: 
		----------
		input_file: PDB with the parameters of the molecule. 

		PELE_version: str
					  path of an executable version of PELE


		Examples: 
		----------

		Load a molecule from a PDB file and minimize it in the vacuum and OBC solvent. 
        >>> import MoleculeMinimized as MM

		>>> new_molecule = MM.MoleculeMinimized('ligand.pdb', '/home/municoy/builds/PELE/PELE-repo_serial/PELE-1.6')
		>>> new_molecule.minimize(input_file = 'ligand.pdb', PELE_version = '/home/municoy/builds/PELE/PELE-repo_serial/PELE-1.6')

		"""
		self.input_file = input_file
		self.PELE_version = PELE_version

	def _output_folder(self,input_file):
		"""
		It creates an output folder with a copy of the ligand's pdb where all the results will be saved. 
		"""
		os.system("mkdir output")
		os.system("cp %s ./output/" % (input_file))
		os.system("mv ./output/%s ./output/ligand.pdb" % (input_file))

	def _generate_parameters(self):
		"""
		It generates the parameters of the molecule (from the input_file) as DataLocal in the output folder.
		"""
		os.system("python /home/lauramalo/repos/offpele/offpele/main.py ./output/ligand.pdb --with_solvent --as_DataLocal")
		os.system("mv ./DataLocal/ ./output/DataLocal/")

	def _link_folders(self):
		"""
		It links the encessary folders to the output folder.
		"""
		os.system("ln -s /home/municoy/repos/PELE-repo/Data ./output/")
		os.system("ln -s /home/municoy/repos/PELE-repo/Documents ./output/") 

	def minimize(self,input_file, PELE_version):
		"""
		It minimized the molecule with the open forcefield (OFFPELE). 

		Parameters:
		----------
		input_file: PDB with the parameters of the molecule. 

		PELE_version: str
					  path of an executable version of PELE

		"""
		self._output_folder(input_file)
		self._link_folders()
		self._generate_parameters()
		os.chdir("./output/")
		os.system(" %s /home/lauramalo/tests/geometry/VACUUM_minimization.conf > VACUUM_minimization.out" % (PELE_version))
		os.system(" %s /home/lauramalo/tests/geometry/OBC_minimization.conf > OBC_minimization.out" % (PELE_version))