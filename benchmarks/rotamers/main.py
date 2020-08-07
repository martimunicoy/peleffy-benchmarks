
"""
This module is designed in order to run the different test for the rotamers
libraries between using PlopRotTemp and offpele.
"""

import os
import argparse
from test_rotamers import TestRotamers


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
    parser.add_argument("-test" "--test_type", type=str, default="bonds",
                        help="Type of test to perform with the "
                        + "rotamers libreries")
    parser.add_argument("-o", "--output_path", type=str, default="",
                        help="Path where to write the results")
    parser.add_argument("-of", "--output_file", type=str, default="diff.csv",
                        help="File with differences in the two files, "
                        + "if there are")
    parser.add_argument("-t", "--test_type", metavar="NAME", type=str,
                        default="groups",
                        help="Type of test to perform with the rotamers "
                        + "libraries")
    args = parser.parse_args()

    return args


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


def main(mae_input_file, pdb_input_file, test_type, output_path,
         output_file):
    """
    Main function that performs the rotamer libraries tests.
    """
    if output_path != "":
        output_file = os.join(output_path, output_file)
    with open(output_file, 'w') as f:
        perform_test(f, test_type, mae_input_file, pdb_input_file)


if __name__ == "__main__":
    args = parse_args()
    main(args.mae_input_file, args.pdb_input_file, args.test_type,
         args.out_path, args.out_file)
