#!/bin/bash
#PBS -N _biox3_scratch_users_rhiju_projects_stepwise_benchmark_new_rna_res_level_energy4_RPT3_gagua_pentaloop_SWM_2
#PBS -o /dev/null
#PBS -e /dev/null
#PBS -m n
#PBS -M nobody@stanford.edu
#PBS -l walltime=48:00:00

cd /biox3/scratch/users/rhiju/projects/stepwise/benchmark/new/rna_res_level_energy4_RPT3/gagua_pentaloop

/home/rhiju/src/rosetta//main/source/bin/stepwise -s gagua_pentaloop_HELIX1.pdb -fasta gagua_pentaloop.fasta -terminal_res A:20 A:28 -extra_min_res A:21 A:27 -superimpose_over_all -cycles 200 -nstruct 20 -intermolecular_frequency 0.0 -save_times -native gagua_pentaloop_1xjr_RNA.pdb -score:weights stepwise/rna/rna_res_level_energy4.wts -cycles 200 -out:file:silent SWM/2/swm_rebuild.out > /dev/null 2> /dev/null 
