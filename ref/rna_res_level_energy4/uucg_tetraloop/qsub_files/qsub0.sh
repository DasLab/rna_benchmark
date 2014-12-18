#!/bin/bash
#PBS -N _biox3_scratch_users_rhiju_projects_stepwise_benchmark_new_rna_res_level_energy4_RPT3_uucg_tetraloop_SWM_0
#PBS -o /dev/null
#PBS -e /dev/null
#PBS -m n
#PBS -M nobody@stanford.edu
#PBS -l walltime=48:00:00

cd /biox3/scratch/users/rhiju/projects/stepwise/benchmark/new/rna_res_level_energy4_RPT3/uucg_tetraloop

/home/rhiju/src/rosetta//main/source/bin/stepwise -s uucg_tetraloop_HELIX1.pdb -fasta uucg_tetraloop.fasta -terminal_res A:4 A:11 -extra_min_res A:5 A:10 -superimpose_over_all -cycles 200 -nstruct 20 -intermolecular_frequency 0.0 -save_times -native uucg_tetraloop_2koc_RNA.pdb -score:weights stepwise/rna/rna_res_level_energy4.wts -cycles 200 -out:file:silent SWM/0/swm_rebuild.out > /dev/null 2> /dev/null 
