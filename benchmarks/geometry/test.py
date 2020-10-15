# Local imports
from getter import QCPortal
from minimize import Minimizer
from offpelebenchmarktools.utils.pele import PELEBaseJob, PELEMinimization
from offpele.topology import Molecule
from utils_benchmark import MinimizationBenchmark
# External imports
import argparse
import os 

benchmark = MinimizationBenchmark( 
        dataset = 'SMIRNOFF Coverage Set 1',
        out_folder = 'SET2')
benchmark.run(filter_structures = True)
