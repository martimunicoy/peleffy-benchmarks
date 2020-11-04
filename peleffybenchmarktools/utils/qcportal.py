"""
It handles queries to the QCPortal.
"""


class QCPortal(object):
    """
    It handles a QCPortal query to download and parse a dataset.
    """

    def __init__(self, n_proc=1):
        """
        It initializes a QCPortal object.

        Parameters
        ----------
        n_proc : int
            The number of processors to employ to gather and parse data
        """
        self.n_proc = n_proc

        # Supress OpenForceField toolkit warnings
        import logging
        logging.getLogger().setLevel(logging.ERROR)

    def _parallel_getter(self, attributes, record_item):
        """
        It gets the attributes from a record and builds the data dict.

        Parameters
        ----------
        attributes : list[str]
            List of attributes to keep from each record
        record_item : tuple[str, qcportal.collections.optimization_dataset.OptEntry]
            The entry to parse

        Returns
        -------
        data : dict
            The dictionary that contains the parsed data
        """
        name, record = record_item

        data = dict()
        record_attributes = dict()

        for attribute in attributes:
            record_attributes[attribute] = record.attributes[attribute]

        data[name] = record_attributes

        return data

    def get_data(self, collection_name,
                 collection_type='OptimizationDataset',
                 attributes=['canonical_isomeric_explicit_hydrogen_smiles']):
        """
        It downloads and parses a certain dataset from the QCPortal.

        Parameters
        ----------
        collection_name : str
            Name of the collection to retrieve
        collection_type : str
            Type of collection to retrieve
        attributes : list[str]
            List of attributes to keep from each record

        Returns
        -------
        data : dict[dict]
            Parsed data containing all the records and their attributes
            from the retrieved dataset
        """

        import qcportal as ptl
        from functools import partial
        from tqdm import tqdm
        from multiprocessing import Pool

        client = ptl.FractalClient()
        ds = client.get_collection(collection_type, collection_name)

        parallel_function = partial(self._parallel_getter, attributes)

        with Pool(self.n_proc) as p:
            listed_data = list(tqdm(p.imap(parallel_function,
                                           ds.data.records.items()),
                                    total=len(ds.data.records.items())))

        data = dict()
        for dictionary_elements in listed_data:
            for key, value in dictionary_elements.items():
                data[key] = value

        return data

    def _parallel_struct_getter(self, ds, output_path, geometry_selection,
                                index):
        """
        It builds the optimized molecule from an Optimization Dataset
        record and saves it to a PDB file.

        Parameters
        ----------
        ds : a qcportal.collections.optimization_dataset.OptimizationDataset
            The optimization dataset the record belongs to
        output_path : str
            The output path where the PDB file will be saved
        geometry_selection : str
            The geometry to feed to the molecule. One of
            ['initial', 'optimized']
        index : int
            The index of the record name to extract. It belongs to the
            index of the list: list(ds.data.records.keys())
        """
        from openforcefield.topology import Molecule
        from simtk import unit
        from rdkit import Chem
        import os

        # Check supplied geometry selection string
        if geometry_selection.lower() not in ['initial', 'optimized']:
            raise NameError('The supplied geometry selection is not valid')

        # Get records and record names
        records = ds.data.records
        record_names = list(ds.data.records.keys())

        # Get record
        record = records[record_names[index]]

        # Get the corresponding mapped smiles and the optimized geometry
        mapped_smiles = self._get_mapped_smiles(record)
        if geometry_selection.lower() == 'initial':
            geometry = self._get_initial_geometry(ds, record)
        else:
            geometry = self._get_optimized_geometry(ds, record)

        if geometry is None:
            print('Warning: the {} geometry of '.format(geometry_selection)
                  + '{}-{} '.format(record_names[index], mapped_smiles)
                  + 'could not be retrieved')
            return

        # Obtain a molecule from mapped smiles
        off_mol = Molecule.from_mapped_smiles(mapped_smiles)

        # Add a conformer with the optimized geometry
        off_mol.add_conformer(geometry.in_units_of(unit.angstrom))

        # Convert it to an RDKit representation and save it to a PDB file
        rdkit_mol = off_mol.to_rdkit()
        Chem.rdmolfiles.MolToPDBFile(
            rdkit_mol,
            os.path.join(output_path, '{}.pdb'.format(index + 1)))

    def get_structures(self, collection_name, output_path='out',
                       collection_type='OptimizationDataset',
                       filter_out_standard_dihedrals=False):
        """
        It writes down the optimized structures from the QCPortal.

        Parameters
        ----------
        collection_name : str
            Name of the collection to retrieve
        output_path : str
            The output path where structures will be saved
        collection_type : str
            Type of collection to retrieve
        filter_out_standard_dihedrals : bool
            Whether to filter out standard dihedrals to only keep
            non-standard ones, or not
        """
        import os
        import qcportal as ptl
        import json
        from functools import partial
        from tqdm import tqdm
        from multiprocessing import Pool

        os.makedirs(output_path, exist_ok=True)

        client = ptl.FractalClient()
        ds = client.get_collection(collection_type, collection_name)

        record_names = list(ds.data.records.keys())

        parallel_function = partial(self._parallel_struct_getter,
                                    ds, output_path)

        with Pool(self.n_proc) as p:
            list(tqdm(p.imap(parallel_function, range(0, len(record_names))),
                      total=len(ds.data.records.items())))

        # Store the index pointing to the corresponding record name
        index_to_name = dict()
        for index, name in enumerate(ds.data.records.keys()):
            index_to_name[index + 1] = '-'.join(name.split('-')[:-1])

        json_output = json.dumps(index_to_name)

        with open(os.path.join(output_path, 'index_to_name.json'), "w") as f:
            f.write(json_output)

    def _get_mapped_smiles(self, record):
        """
        It returns the mapped smiles of an Optimization Dataset record.

        Parameters
        ----------
        record : a qcportal.collections.optimization_dataset.OptEntry
            The optimization entry whose mapped smiles is requested

        Returns
        -------
        mapped_smiles : str
            The mapped smiles defining the molecule of this record
            following a specific atom ordering
        """
        mapped_smiles = \
            record.dict()['attributes'][
                'canonical_isomeric_explicit_hydrogen_mapped_smiles']

        return mapped_smiles

    def _get_optimized_geometry(self, ds, record):
        """
        It obtains the optimized geometry of an Optimization Dataset record.

        Parameters
        ----------
        ds : a qcportal.collections.optimization_dataset.OptimizationDataset
            The optimization dataset the record belongs to
        record : a qcportal.collections.optimization_dataset.OptEntry
            The optimization entry whose optimized geometry is requested

        Returns
        -------
        geometry : a simtk.unit.Quantity object
            The array of 3D coordinates that define the optimitzed geometry
            of the record's molecule. Note that the order of the atoms
            matches with the mapped smiles
        """
        from simtk import unit
        import numpy as np

        optimization_record = ds.get_record(record.dict()['name'],
                                            specification='default')
        try:
            geometry = optimization_record.get_final_molecule().geometry
        except TypeError:
            return None

        geometry = unit.Quantity(np.array(geometry, np.float), unit.bohr)

        return geometry

    def _get_initial_geometry(self, ds, record):
        """
        It obtains the initial geometry of an Optimization Dataset record.

        Parameters
        ----------
        ds : a qcportal.collections.optimization_dataset.OptimizationDataset
            The optimization dataset the record belongs to
        record : a qcportal.collections.optimization_dataset.OptEntry
            The optimization entry whose optimized geometry is requested

        Returns
        -------
        geometry : a simtk.unit.Quantity object
            The array of 3D coordinates that define the initial geometry
            of the record's molecule. Note that the order of the atoms
            matches with the mapped smiles
        """
        from simtk import unit
        import numpy as np

        optimization_record = ds.get_record(record.dict()['name'],
                                            specification='default')
        try:
            geometry = optimization_record.get_initial_molecule().geometry
        except TypeError:
            return None

        geometry = unit.Quantity(np.array(geometry, np.float), unit.bohr)

        return geometry
