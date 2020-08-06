"""
It obtains the unique molecules in a QCArchive dataset according to
the SMARTS tags.
"""

import json


OUTPUT_FOLDER = "QM"
JSON_NAME = "ids_to_smarts.json"


with open('../{}/{}'.format(OUTPUT_FOLDER, JSON_NAME), 'r') as f:
    tags = json.load(f)

unique_tags = set(tags.values())

print("Number of unique molecules:", len(unique_tags))
print("Tags:")
for tag in sorted(unique_tags):
    print(" -", tag)
