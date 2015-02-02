#!/bin/bash
#PBS -N _biox3_home_geniesse_src_stepwise_benchmark_new_tandemGA_rna_res_level_energy4_vary_polar_hydrogen_geometry_2k_cycles_tandemGA-11-uGAg-uGAa_SWM_8
#PBS -o /dev/null
#PBS -e /dev/null
#PBS -m n
#PBS -M nobody@stanford.edu
#PBS -l walltime=72:00:00

cd /biox3/home/geniesse/src/stepwise_benchmark/new/tandemGA/rna_res_level_energy4_vary_polar_hydrogen_geometry_2k_cycles/tandemGA-11-uGAg-uGAa

/home/geniesse/src/rosetta//main/source/bin/stepwise -s tandemGA-11-uGAg-uGAa_HELIX1.pdb tandemGA-11-uGAg-uGAa_HELIX2.pdb -native tandemGA-11-uGAg-uGAa_1mis_RNA_11-uGAg-uGAa.pdb -terminal_res A:1 A:16 -extra_min_res A:3 A:6 A:11 A:14 -superimpose_over_all -fasta tandemGA-11-uGAg-uGAa.fasta -nstruct 20 -intermolecular_frequency 0.0 -save_times -score:weights stepwise/rna/rna_res_level_energy4.wts -cycles 2000 -vary_polar_hydrogen_geometry true -out:file:silent SWM/8/swm_rebuild.out > 8.out 2> 8.err 
