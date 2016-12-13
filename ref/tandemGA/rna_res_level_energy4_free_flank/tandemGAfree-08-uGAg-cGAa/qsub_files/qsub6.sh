#!/bin/bash
#PBS -N _biox3_home_geniesse_src_stepwise_benchmark_new_tandemGA_rna_res_level_energy4_free_flank_tandemGAfree-08-uGAg-cGAa_SWM_6
#PBS -o /dev/null
#PBS -e /dev/null
#PBS -m n
#PBS -M nobody@stanford.edu
#PBS -l walltime=48:00:00

cd /biox3/home/geniesse/src/stepwise_benchmark/new/tandemGA/rna_res_level_energy4_free_flank/tandemGAfree-08-uGAg-cGAa

/home/geniesse/src/rosetta//main/source/bin/stepwise -s tandemGAfree-08-uGAg-cGAa_HELIX1.pdb tandemGAfree-08-uGAg-cGAa_HELIX2.pdb -native tandemGAfree-08-uGAg-cGAa_1mis_RNA_08-uGAg-cGAa.pdb -terminal_res A:1 A:16 -extra_min_res A:2 A:7 A:10 A:15 -superimpose_over_all -fasta tandemGAfree-08-uGAg-cGAa.fasta -cycles 200 -nstruct 20 -intermolecular_frequency 0.0 -save_times -score:weights stepwise/rna/rna_res_level_energy4.wts -cycles 200  -out:file:silent SWM/6/swm_rebuild.out > 6.out 2> 6.err 
