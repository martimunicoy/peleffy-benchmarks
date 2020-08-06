"""
It extracts the xyz coordinate files for all the ligands found in a
certain QCArchive dataset. It also constructs a json file with
pairs of indexes and SMARTS tags.
"""

import os
import json
import qcportal as ptl
from functools import partial
from multiprocessing import Pool
from tqdm import tqdm


NUMBER_OF_PROCESSORS = 8
DATASET_NAME = "SMIRNOFF Coverage Set 1"
OUTPUT_FOLDER = "QM"
JSON_NAME = "ids_to_smarts.json"


def extract_xyz_file(record_tags, index):
    record_tag = record_tags[index]
    try:
        opt_record = ds.get_record(name=record_tag, specification="default")
        opt_molecule = opt_record.get_final_molecule()
    except TypeError:
        print(' - Skipping ligand', record_tag)
        return None

    opt_molecule.to_file('../{}/{}.xyz'.format(OUTPUT_FOLDER, index))

    return record_tag


client = ptl.FractalClient()
ds = client.get_collection("OptimizationDataset", DATASET_NAME)

record_tags = ds.df.index.to_list()

os.makedirs("../{}".format(OUTPUT_FOLDER), exist_ok=True)

ids_to_smarts = dict()
results = list()

with Pool(NUMBER_OF_PROCESSORS) as pool:
    parallel_function = partial(extract_xyz_file, record_tags)
    for result in tqdm(pool.imap(func=parallel_function, iterable=range(0, len(record_tags))), total=len(record_tags)):
        results.append(result)

for index, tag in enumerate(results):
    if tag is not None:
        assert '-' in tag, 'Unexpected record tag {}'.format(tag)
        ids_to_smarts[index] = ''.join(tag.split('-')[:-1])

json_output = json.dumps(ids_to_smarts)

with open("../{}/{}".format(OUTPUT_FOLDER, JSON_NAME), "w") as f:
    f.write(json_output)

