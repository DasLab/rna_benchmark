#!/bin/bash
#PBS -N _biox3_scratch_users_rhiju_projects_stepwise_benchmark_new_rna_res_level_energy4_RPT3_srp_domainIV_SWM_2
#PBS -o /dev/null
#PBS -e /dev/null
#PBS -m n
#PBS -M nobody@stanford.edu
#PBS -l walltime=48:00:00

cd /biox3/scratch/users/rhiju/projects/stepwise/benchmark/new/rna_res_level_energy4_RPT3/srp_domainIV

/home/rhiju/src/rosetta//main/source/bin/stepwise -s srp_domainIV_HELIX1.pdb srp_domainIV_HELIX2.pdb -fasta srp_domainIV.fasta -terminal_res A:3 A:10 B:15 B:22 -extra_min_res A:4 A:9 B:16 B:21 -superimpose_over_all -cycles 200 -nstruct 20 -intermolecular_frequency 0.0 -save_times -native srp_domainIV_native_1lnt_RNA.pdb -score:weights stepwise/rna/rna_res_level_energy4.wts -cycles 200 -out:file:silent SWM/2/swm_rebuild.out > /dev/null 2> /dev/null 
