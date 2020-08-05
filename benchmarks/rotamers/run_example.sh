


python /home/lauramalo/repos/offpele/offpele/main.py /home/lauramalo/repos/offpele-benchmarks/benchmarks/Plop_ligands/pdb/3.pdb -o /home/lauramalo/tests/results/pdb/

PlopRotTemp /home/lauramalo/repos/offpele-benchmarks/benchmarks/Plop_ligands/mae/3.mae  --rotamerdir /home/lauramalo/tests/results/mae/

python /home/lauramalo/repos/offpele-benchmarks/benchmarks/rotamers/main.py  /home/lauramalo/tests/results/mae/LIG.rot.assign /home/lauramalo/tests/results/pdb/LIG.rot.assign -of /home/lauramalo/tests/results/diff/

rm /home/lauramalo/tests/results/pdb/LIG.rot.assign

rm /home/lauramalo/tests/results/mae/LIG.rot.assign
