#!/bin/bash
#PBS -N _biox3_scratch_users_rhiju_projects_stepwise_benchmark_new_rna_res_level_energy4_RPT3_srl_fixed_SWM_5
#PBS -o /dev/null
#PBS -e /dev/null
#PBS -m n
#PBS -M nobody@stanford.edu
#PBS -l walltime=48:00:00

cd /biox3/scratch/users/rhiju/projects/stepwise/benchmark/new/rna_res_level_energy4_RPT3/srl_fixed

/home/rhiju/src/rosetta//main/source/bin/stepwise -s srl_fixed_START1_1q9a_RNA.pdb -fasta srl_fixed.fasta -terminal_res A:651 A:669 -extra_min_res A:652 A:658 A:663 A:668 -cycles 200 -nstruct 20 -intermolecular_frequency 0.0 -save_times -native srl_fixed_1q9a_RNA.pdb -score:weights stepwise/rna/rna_res_level_energy4.wts -cycles 200 -out:file:silent SWM/5/swm_rebuild.out > /dev/null 2> /dev/null 
