#!/bin/bash
#PBS -N _biox3_home_geniesse_src_stepwise_benchmark_new_tandemGA_rna_res_level_energy4_vary_polar_hydrogen_geometry_2k_cycles_tandemGA-04-uGAu-aGAa_SWM_9
#PBS -o /dev/null
#PBS -e /dev/null
#PBS -m n
#PBS -M nobody@stanford.edu
#PBS -l walltime=72:00:00

cd /biox3/home/geniesse/src/stepwise_benchmark/new/tandemGA/rna_res_level_energy4_vary_polar_hydrogen_geometry_2k_cycles/tandemGA-04-uGAu-aGAa

/home/geniesse/src/rosetta//main/source/bin/stepwise -s tandemGA-04-uGAu-aGAa_HELIX1.pdb tandemGA-04-uGAu-aGAa_HELIX2.pdb -native tandemGA-04-uGAu-aGAa_1mis_RNA_04-uGAu-aGAa.pdb -terminal_res A:1 A:16 -extra_min_res A:3 A:6 A:11 A:14 -superimpose_over_all -fasta tandemGA-04-uGAu-aGAa.fasta -nstruct 20 -intermolecular_frequency 0.0 -save_times -score:weights stepwise/rna/rna_res_level_energy4.wts -cycles 2000 -vary_polar_hydrogen_geometry true -out:file:silent SWM/9/swm_rebuild.out > 9.out 2> 9.err 
