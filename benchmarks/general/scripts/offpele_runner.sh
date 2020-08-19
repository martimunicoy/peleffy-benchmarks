#!/bin/bash

# Author: Martí Municoy
# Date: August 7th, 2020
# Description: It calls offpele to generate the parameters for all the supplied ligand structures.
# Usage: ./PlopRotTemp_runner.sh path_to_maes/*.mae

if [[ $# -eq 0 ]] ; then
    echo 'The path to one or multiple mae files must be supplied'
    exit 0
fi

for pdb_file in "$@"
do
	if [ ! -f $pdb_file ] | [[ $pdb_file != *.pdb ]]
	then
		echo "Warning! Skipping $pdb_file because it is not a valid file"
	else
		file_name=$(basename -s .pdb $pdb_file)
		mkdir -p offpele_out/$file_name
		python -m offpele.main $pdb_file -o offpele_out/$file_name -c gasteiger
	fi
	echo " - Template and rotamer library for $pdb_file were generated successfully"
done

