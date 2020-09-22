
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
		from pathlib import Path
		import shutil
		import os 


		# Handels path of the input file
		p = Path(input_file)
		file, folder = p.name, p.parents[0]

		# It makes the output directory
		os.makedirs('output',exist_ok = True)
		shutil.copy(p, os.path.join(os.getcwd(),'output', 'ligand.pdb'))


	def _generate_parameters(self):
		"""
		It generates the parameters of the molecule (from the input_file) as DataLocal in the output folder.
		"""
		import offpele
		from offpele.topology import Molecule
		from offpele.template import Impact
		from offpele.solvent import OBC2
		from offpele.main import handle_output_paths
		import os 

	
		# Forcefield and charges method
		forcefield = 'openff_unconstrained-1.2.0.offxml'
		charges_method = 'am1bcc'

		# Create representation of a particular molecule
		PATH_molecule = os.path.join(os.getcwd(),'output', 'ligand.pdb')
		molecule = Molecule(PATH_molecule)

		# Saving paths
		rotamer_library_output_path, impact_output_path, solvent_output_path = \
			handle_output_paths(molecule = molecule, output =os.path.join(os.getcwd(),'output'), as_datalocal = True )

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


	def _link_folders(self):
		"""
		It links the encessary folders to the output folder.
		"""
		import os 

		PELE_SRC = '/home/municoy/repos/PELE-repo/'
		
		os.symlink('{}Data'.format(PELE_SRC), os.path.join(os.getcwd(),'output','Data'))
		os.symlink('{}Documents'.format(PELE_SRC), os.path.join(os.getcwd(),'output', 'Documents'))


	def minimize(self,input_file, PELE_version):
		"""
		It minimized the molecule with the open forcefield (OFFPELE). 

		Parameters:
		----------
		input_file: PDB with the parameters of the molecule. 

		PELE_version: str
					  path of an executable version of PELE

		"""
		import os
		import requests
		from pathlib import Path

		VACUUM_CF = '/home/lauramalo/repos/offpele-benchmarks/benchmarks/solvent/Conf/VACUUM_minimization.conf'
		OBC_CF = '/home/lauramalo/repos/offpele-benchmarks/benchmarks/solvent/Conf/OBC_minimization.conf'

		self._output_folder(input_file)
		self._link_folders()
		self._generate_parameters()
		
		# Minimization
		os.chdir("./output/")
		os.system(" %s %s > VACUUM_minimization.out" % (PELE_version, VACUUM_CF))
		os.system(" %s %s > OBC_minimization.out" % (PELE_version, OBC_CF))

		# Rename the output folder to the molecule name
		os.chdir("..")
		p = Path(input_file)
		file, folder = p.name, p.parents[0]
		new_folder = os.path.splitext(file)[0]
		os.rename('output', new_folder)



