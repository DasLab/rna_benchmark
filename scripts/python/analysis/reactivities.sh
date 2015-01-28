#!/bin/bash


PDB_PREFIX=$1

for PDB in ${PDB_PREFIX}*; do
    echo ${PDB}
    ./reactivity.py ${PDB}
done