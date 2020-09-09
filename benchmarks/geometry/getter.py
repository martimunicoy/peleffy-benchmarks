"""
It prepares the system files and runs a PELE minimization.
"""


class QCPortal(object):
    """
    It handles a QCPortal query to download and parse a dataset.
    """

    def __init__(self):
        """
        It initializes a QCPortal object.
        """
        pass

    def get_data(self, collection_type, collection_name, attributes):
        """
        It downloads and parses a certain dataset from the QCPortal.

        Parameters
        ----------
        collection_type : str
            Type of collection to retrieve
        collection_name : str
            Name of the collection to retrieve
        attributes : list[str]
            List of attributes to keep from each record

        Returns
        -------
        data : dict[dict]
            Parsed data containing all the records and their attributes
            from the retrieved dataset
        """

        import qcportal as ptl

        client = ptl.FractalClient()
        ds = client.get_collection(collection_type, collection_name)

        data = dict()

        for name, record in ds.data.records.items():
            record_attributes = dict()
            for attribute in attributes:
                record_attributes[attribute] = record.attributes[attribute]

            data[name] = record_attributes

        return data

    def get_structures(collection_type, collection_name, output_path):
        """
        It writes down the optimized structures from the QCPortal.

        Parameters
        ----------
        collection_type : str
            Type of collection to retrieve
        collection_name : str
            Name of the collection to retrieve
        output_path : str
            The output path where structures will be saved
        """
        import os
        import qcportal as ptl
        import json

        os.mkdirs(output_path, exist_ok=True)

        client = ptl.FractalClient()
        ds = client.get_collection(collection_type, collection_name)

        index = int(1)
        index_to_name = dict()
        for name, record in ds.data.records.items():
            optimization_record = ds.get_record(name, specification='default')
            opt_molecule = optimization_record.get_final_molecule()
            opt_molecule.to_file('{}/{}.xyz'.format(output_path, index))
            index_to_name[index] = name
            index += int(1)

        json_output = json.dumps(index_to_name)

        with open(os.path.join(output_path, 'index_to_name.json'), "w") as f:
            f.write(json_output)
