#!/bin/bash
#PBS -N _biox3_scratch_users_rhiju_projects_stepwise_benchmark_new_rna_res_level_energy4_RPT3_tandem_ga_sheared_SWM_9
#PBS -o /dev/null
#PBS -e /dev/null
#PBS -m n
#PBS -M nobody@stanford.edu
#PBS -l walltime=48:00:00

cd /biox3/scratch/users/rhiju/projects/stepwise/benchmark/new/rna_res_level_energy4_RPT3/tandem_ga_sheared

/home/rhiju/src/rosetta//main/source/bin/stepwise -s tandem_ga_sheared_HELIX1.pdb tandem_ga_sheared_HELIX2.pdb -fasta tandem_ga_sheared.fasta -terminal_res A:2 A:7 B:10 B:15 -extra_min_res A:3 A:6 B:11 B:14 -superimpose_over_all -cycles 200 -nstruct 20 -intermolecular_frequency 0.0 -save_times -native tandem_ga_sheared_1yfv_RNA.pdb -score:weights stepwise/rna/rna_res_level_energy4.wts -cycles 200 -out:file:silent SWM/9/swm_rebuild.out > /dev/null 2> /dev/null 
