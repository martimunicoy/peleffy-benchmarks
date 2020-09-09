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

        for name, record in ds.df.records.items():
            record_attributes = dict()
            for attribute in attributes:
                record_attributes[attribute] = record.attributes[attribute]

            data[name] = record_attributes

        return data
