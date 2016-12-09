#!/bin/bash
#PBS -N _biox3_home_geniesse_src_stepwise_benchmark_new_tandemGA_rna_res_level_energy4_free_flank_vary_polar_hydrogen_geometry_2k_cycles_tandemGAfree-36-uGAu-gGAg_SWM_0
#PBS -o /dev/null
#PBS -e /dev/null
#PBS -m n
#PBS -M nobody@stanford.edu
#PBS -l walltime=72:00:00

cd /biox3/home/geniesse/src/stepwise_benchmark/new/tandemGA/rna_res_level_energy4_free_flank_vary_polar_hydrogen_geometry_2k_cycles/tandemGAfree-36-uGAu-gGAg

/home/geniesse/src/rosetta//main/source/bin/stepwise -s tandemGAfree-36-uGAu-gGAg_HELIX1.pdb tandemGAfree-36-uGAu-gGAg_HELIX2.pdb -native tandemGAfree-36-uGAu-gGAg_1mis_RNA_36-uGAu-gGAg.pdb -terminal_res A:1 A:16 -extra_min_res A:2 A:7 A:10 A:15 -superimpose_over_all -fasta tandemGAfree-36-uGAu-gGAg.fasta -nstruct 20 -intermolecular_frequency 0.0 -save_times -score:weights stepwise/rna/rna_res_level_energy4.wts -cycles 2000 -vary_polar_hydrogen_geometry true -out:file:silent SWM/0/swm_rebuild.out > 0.out 2> 0.err 
