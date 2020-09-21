"""
This module contains miscellaneous set of handy classes and functions.
"""


__all__ = ["get_data_file_path"]


from pkg_resources import resource_filename
import os


def get_data_file_path(relative_path):
    """
    It returns the path in the package's data location.

    Parameters
    ----------
    relative_path : str
        The relative path to the file that is required

    Returns
    -------
    output_path : str
        The path in the package's data location, if found
    """
    output_path = resource_filename('offpelebenchmarktools', os.path.join(
        'data', relative_path))

    if not os.path.exists(output_path):
        raise ValueError(
            "Sorry! {output_path} does not exist. If you just added it, "
            + "you'll have to re-install")

    return output_path
