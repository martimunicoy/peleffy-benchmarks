#!/bin/bash

# Author: Martí Municoy
# Date: July 28th, 2020
# Description: It fixes the connects of a set of structures written in PDB files.
# Usage: ./fix_connects.sh path_to_pdbs/*.pdb

# Stop when any fail occurs
set -e

if [ -z ${SCHRODINGER+x} ]
then
	echo "SCHRODINGER variable is unset. Please, specify SCHRODINGER's root directory using this variable."i
	exit 0
else
	echo "Using Schrodinger installation from $SCHRODINGER"
fi

if [[ $# -eq 0 ]] ; then
    echo 'The path to one or multiple PDB files must be supplied'
    exit 0
fi

PDB_OUT="pdb/"
MAE_OUT="mae/"

rm -rf $PDB_OUT $MAE_OUT
mkdir -p $PDB_OUT $MAE_OUT

for pdb_file in "$@"
do
	if [ ! -f $pdb_file ] | [[ $pdb_file != *.pdb ]]
	then
		echo "Warning! Skipping $pdb_file because it is not a valid file"
	else
		absolute_path=$(realpath $pdb_file)
		file_name=$(basename -s .pdb $pdb_file)
		$SCHRODINGER/utilities/prepwizard $absolute_path $file_name.pdb -noepik -noprotassign -noccd -noimpref > /dev/null
		while [ ! -f $file_name.pdb ]; do sleep 1; done
		sleep 1
		mv $file_name.pdb $PDB_OUT
		rm -f $file_name.log
		$SCHRODINGER/utilities/prepwizard $absolute_path $file_name.mae -noepik -noprotassign -noccd -noimpref > /dev/null
		while [ ! -f $file_name.mae ]; do sleep 1; done
                sleep 1
                mv $file_name.mae $MAE_OUT
                rm -f $file_name.log
	fi
	echo " - Connects to $pdb_file were added successfully"
done
