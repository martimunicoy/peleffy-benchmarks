import json

with open('ids_to_smarts.json') as f:
    results = f.read()

results_dict = json.loads(results)

for index, smiles in results_dict.items():
    if smiles == "C[C@H](c1[nH][nH]c(=O)n1)NC#N":
        print(index)
