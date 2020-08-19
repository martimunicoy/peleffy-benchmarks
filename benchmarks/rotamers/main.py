"""
This module is designed in order to run the different test for the rotamers libraries between using
PlopRotTemp and offpele.
"""

import argparse
from test_rotamers import TestRotamers
import utils
import os.path as os


def parse_args():
    """
        Parse command line arguments
        :returns: object -- Object containing command line options
    """
    parser = argparse.ArgumentParser(description="Compute Statistics for "
                                     + "the merged files")
    parser.add_argument("rotamer_schrodinguer", type=str,
                        help="File obtained by PELE with Schodinguer")
    parser.add_argument("rotamer_offpele", type=str,
                        help="File obtained by PELE with offpele")
    parser.add_argument("-o", "--output_path", type=str, default="",
                        help="Path where to write the results")
    parser.add_argument("-of", "--output_file", type=str, default="diff.csv",
                        help="File with differences in the two files, if "
                        + "there are")
    parser.add_argument("-t", "--test_type", metavar="NAME", type=str,
                        default="resolution", help="Type of test to "
                        + "perform with the rotamers libreries")

    args = parser.parse_args()
    return args.rotamer_schrodinguer, args.rotamer_offpele, \
        args.output_path, args.output_file, args.test_type


def perform_test(file, test_type, mae_input_file, pdb_input_file):
    """
    Creates the object rotamers with the rotameres library obtained with
    offpele and with PlopRotTemp and performs the selected test on the data.
    """
    rotamers = TestRotamers(mae_input_file, pdb_input_file, file)
    if test_type == 'bonds':
        rotamers.test_bonds()

    if test_type == 'groups':
        rotamers.test_groups()

    if test_type == 'resolution':
        rotamers.test_resolution()

    if test_type == 'all':
        rotamers.test_bonds()
        rotamers.test_resolution()
        rotamers.test_groups()


def main(mae_input_file, pdb_input_file, output_path, output_file, test_type):
    """
    Main function that performs the rotamer libraries tests.
    """
    if output_path != "":
        output_file = os.join(output_path, output_file)

    with open(output_file, 'w') as f:
        with open(mae_input_file, 'r') as mae_file:
            with open(pdb_input_file, 'r') as pdb_file:
                perform_test(f, test_type, mae_file, pdb_file)


if __name__ == "__main__":
    mae_input_file, pdb_input_file, out_path, out_file, test_file = \
        parse_args()
    main(mae_input_file, pdb_input_file, out_path, out_file, test_file)
