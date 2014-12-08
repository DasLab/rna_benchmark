#!/bin/bash
#PBS -N _biox3_scratch_users_rhiju_projects_stepwise_benchmark_new_rna_res_level_energy4_RPT3_j44a_p4p6_SWM_8
#PBS -o /dev/null
#PBS -e /dev/null
#PBS -m n
#PBS -M nobody@stanford.edu
#PBS -l walltime=48:00:00

cd /biox3/scratch/users/rhiju/projects/stepwise/benchmark/new/rna_res_level_energy4_RPT3/j44a_p4p6

/home/rhiju/src/rosetta//main/source/bin/stepwise -s j44a_p4p6_HELIX1.pdb j44a_p4p6_HELIX2.pdb -fasta j44a_p4p6.fasta -terminal_res A:111 A:118 A:203 A:209 -extra_min_res A:112 A:116 A:205 A:208 -superimpose_over_all -cycles 200 -nstruct 20 -intermolecular_frequency 0.0 -save_times -native j44a_p4p6_1gid_RNAA.pdb -score:weights stepwise/rna/rna_res_level_energy4.wts -cycles 200 -out:file:silent SWM/8/swm_rebuild.out > /dev/null 2> /dev/null 
