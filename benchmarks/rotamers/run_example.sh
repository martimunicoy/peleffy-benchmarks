


python /home/lauramalo/repos/offpele/offpele/main.py /home/lauramalo/repos/offpele-benchmarks/benchmarks/Plop_ligands/pdb/[file_name.pdb] -o RIG.rot.temp /home/lauramalo/tests/results/pdb/

PlopRotTemp /home/lauramalo/repos/offpele-benchmarks/benchmarks/Plop_ligands/mae/[file_name.mae]  --rotamerdir /home/lauramalo/tests/results/mae/

python /home/lauramalo/repos/offpele-benchmarks/benchmarks/rotamers/main.py /home/lauramalo/tests/results/pdb/LIG.rot.assign /home/lauramalo/tests/results/mae/LIG.rot.assign -of /home/lauramalo/tests/results/diff/

rm /home/lauramalo/tests/results/pdb/LIG.rot.assign

rm /home/lauramalo/tests/results/mae/LIG.rot.assign
