#!/bin/bash


PDB_PREFIX=$1

for PDB in ${PDB_PREFIX}*; do
    echo ${PDB}
    ~/src/stepwise_benchmark/scripts/python/analysis/reactivity.py ${PDB}
done