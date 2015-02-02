#!/bin/bash
#PBS -N _biox3_home_geniesse_src_stepwise_benchmark_new_tandemGA_rna_res_level_energy4_tandemGA-19-gGAc-gGAc_SWM_7
#PBS -o /dev/null
#PBS -e /dev/null
#PBS -m n
#PBS -M nobody@stanford.edu
#PBS -l walltime=48:00:00

cd /biox3/home/geniesse/src/stepwise_benchmark/new/tandemGA/rna_res_level_energy4/tandemGA-19-gGAc-gGAc

/home/geniesse/src/rosetta//main/source/bin/stepwise -s tandemGA-19-gGAc-gGAc_HELIX1.pdb tandemGA-19-gGAc-gGAc_HELIX2.pdb -native tandemGA-19-gGAc-gGAc_1mis_RNA_19-gGAc-gGAc.pdb -terminal_res A:1 A:16 -extra_min_res A:3 A:6 A:11 A:14 -superimpose_over_all -fasta tandemGA-19-gGAc-gGAc.fasta -cycles 200 -nstruct 20 -intermolecular_frequency 0.0 -save_times -score:weights stepwise/rna/rna_res_level_energy4.wts -cycles 200  -out:file:silent SWM/7/swm_rebuild.out > 7.out 2> 7.err 
