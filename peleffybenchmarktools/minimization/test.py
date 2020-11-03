# Local imports
from getter import QCPortal
from minimize import Minimizer
from peleffybenchmarktools.utils.pele import PELEBaseJob, PELEMinimization
from peleffy.topology import Molecule
from utils_benchmark import MinimizationBenchmark
# External imports
import argparse
import os

benchmark = MinimizationBenchmark(
        dataset = 'SMIRNOFF Coverage Set 1',
        out_folder = 'SET3')
#benchmark.run(filter_structures = True)
benchmark.compute_RMSD()
