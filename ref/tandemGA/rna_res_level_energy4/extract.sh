#!/bin/bash

MODELDIR=models
DIRS=$1
OUTFILE=$2

for DIR in ${DIRS}*; do
    cd ${DIR}
    extract_lowscore_decoys.py ${OUTFILE} $3
    for PDB in ${OUTFILE}.*.pdb; do
	cp ${PDB} ../${MODELDIR}/${DIR}_${PDB}
    done
    cd ..
done
