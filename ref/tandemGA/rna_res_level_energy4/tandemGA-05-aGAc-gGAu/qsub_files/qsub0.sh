#!/bin/bash
#PBS -N _biox3_home_geniesse_src_stepwise_benchmark_new_tandemGA_rna_res_level_energy4_tandemGA-05-aGAc-gGAu_SWM_0
#PBS -o /dev/null
#PBS -e /dev/null
#PBS -m n
#PBS -M nobody@stanford.edu
#PBS -l walltime=48:00:00

cd /biox3/home/geniesse/src/stepwise_benchmark/new/tandemGA/rna_res_level_energy4/tandemGA-05-aGAc-gGAu

/home/geniesse/src/rosetta//main/source/bin/stepwise -s tandemGA-05-aGAc-gGAu_HELIX1.pdb tandemGA-05-aGAc-gGAu_HELIX2.pdb -native tandemGA-05-aGAc-gGAu_1mis_RNA_05-aGAc-gGAu.pdb -terminal_res A:1 A:16 -extra_min_res A:3 A:6 A:11 A:14 -superimpose_over_all -fasta tandemGA-05-aGAc-gGAu.fasta -cycles 200 -nstruct 20 -intermolecular_frequency 0.0 -save_times -score:weights stepwise/rna/rna_res_level_energy4.wts -cycles 200  -out:file:silent SWM/0/swm_rebuild.out > 0.out 2> 0.err 
