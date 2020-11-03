"""
It takes to peleffy's molecule representations, aligns them and calculates
the RMSD using RDKit.
"""


from multiprocessing import Pool
from copy import deepcopy

from rdkit.Chem import AllChem
from rdkit import Chem


class Aligner(object):
    """
    It aligns two molecule representations as much as possible and
    calculates the resulting RMSD difference.
    """

    def __init__(self, prb_mol, ref_mol, n_proc=1, max_iter=1000,
                 threshold=1.0):
        """
        It initializes an Aligner object.

        Parameters
        ----------
        prb_mol : an rdkit.Chem.rdchem.Mol object
            The RDKit's Molecule object to use as probe
        ref_mol : an rdkit.Chem.rdchem.Mol object
            The RDKit's Molecule object to use as reference
        n_proc : int
            The number of processors to employ for such alignment
        max_iter : int
            The maximum number of conformations to generate during the
            alignment process
        threshold : float
            The minimum threshold that is used as exit condition
        """
        # First checkings
        if not isinstance(prb_mol, Chem.rdchem.Mol):
            raise TypeError('Wrong type for probe molecule')
        if not isinstance(ref_mol, Chem.rdchem.Mol):
            raise TypeError('Wrong type for reference molecule')

        self._prb_mol = prb_mol
        self._ref_mol = ref_mol
        self._n_proc = n_proc
        self._max_iter = max_iter
        self._threshold = threshold
        self._results = None

    def _get_lowest_RMSD(self, seed):
        prb_mol = deepcopy(self._prb_mol)
        _ = AllChem.EmbedMultipleConfs(prb_mol, 15, randomSeed=seed)

        best_RMSD = None
        best_conf = None
        for i in range(0, prb_mol.GetNumConformers()):
            rmsd = Chem.rdMolAlign.AlignMol(prbMol=prb_mol,
                                            refMol=self._ref_mol,
                                            prbCid=i)
            if best_RMSD is None or rmsd < best_RMSD:
                best_RMSD = rmsd
                best_conf = i

        return {'best_RMSD': best_RMSD, 'best_conformation': best_conf,
                'molecule': prb_mol}

    def align(self):

        match_found = False
        self._results = list()

        # TODO save all results

        for iteration_id in range(0, self._max_iter, self._n_proc):
            with Pool(self._n_proc) as pool:
                results = pool.map(self._get_lowest_RMSD,
                                   [iteration_id + i + 1 for i
                                    in range(0, self._n_proc)])

                for result in results:
                    self._results.append(result)
                    if result['best_RMSD'] < self._threshold:
                        match_found = True

                if match_found:
                    break

    def get_results(self):
        if self._results is not None:
            return self._results
        else:
            return None

    def get_best_results(self):
        if self._results is not None:
            return sorted(self._results, key=lambda d: d['best_RMSD'])[0]
        else:
            return None

    def to_pdb(self, prb_pdb_name, ref_pdb_name):
        best_results = self.get_best_results()
        mol = best_results['molecule']
        cid = best_results['best_conformation']

        Chem.rdmolfiles.MolToPDBFile(mol, prb_pdb_name, confId=cid)
        Chem.rdmolfiles.MolToPDBFile(self._ref_mol, ref_pdb_name)
