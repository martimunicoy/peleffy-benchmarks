#!/bin/bash

# Author: Martí Municoy
# Date: August 7th, 2020
# Description: It calls PlopRotTemp to generate the parameters for all the supplied ligand structures.
# Usage: ./PlopRotTemp_runner.sh path_to_maes/*.mae

if [ -z ${SCHRODINGER+x} ]
then
	echo "SCHRODINGER variable is unset. Please, specify SCHRODINGER's root directory using this variable."i
	exit 0
else
	echo "Using Schrodinger installation from $SCHRODINGER"
fi

if [ -z ${PLOPROTTEMP+x} ]
then
        echo "PLOPROTTEMP variable is unset. Please, specify PLOPROTTEMP's root directory using this variable."i
        exit 0
else
        echo "Using PlopRotTemp installation from $PLOPROTTEMP"
fi

if [[ $# -eq 0 ]] ; then
    echo 'The path to one or multiple mae files must be supplied'
    exit 0
fi

for mae_file in "$@"
do
	if [ ! -f $mae_file ] | [[ $mae_file != *.mae ]]
	then
		echo "Warning! Skipping $mae_file because it is not a valid file"
	else
		file_name=$(basename -s .mae $mae_file)
		mkdir -p PlopRotTemp_out/$file_name
		$SCHRODINGER/utilities/python $PLOPROTTEMP/PlopRotTemp/main.py $mae_file --rotamerdir PlopRotTemp_out/$file_name --templatedir PlopRotTemp_out/$file_name
	fi
	echo " - Template and rotamer library for $mae_file were generated successfully"
done

