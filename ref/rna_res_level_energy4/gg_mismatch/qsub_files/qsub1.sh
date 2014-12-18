#!/bin/bash
#PBS -N _biox3_scratch_users_rhiju_projects_stepwise_benchmark_new_rna_res_level_energy4_RPT3_gg_mismatch_SWM_1
#PBS -o /dev/null
#PBS -e /dev/null
#PBS -m n
#PBS -M nobody@stanford.edu
#PBS -l walltime=48:00:00

cd /biox3/scratch/users/rhiju/projects/stepwise/benchmark/new/rna_res_level_energy4_RPT3/gg_mismatch

/home/rhiju/src/rosetta//main/source/bin/stepwise -s gg_mismatch_HELIX1.pdb gg_mismatch_HELIX2.pdb -fasta gg_mismatch.fasta -terminal_res A:2 A:6 B:5 B:9 -extra_min_res A:3 A:5 B:6 B:8 -superimpose_over_all -cycles 200 -nstruct 20 -intermolecular_frequency 0.0 -save_times -native gg_mismatch_1f5g.pdb -score:weights stepwise/rna/rna_res_level_energy4.wts -cycles 200 -out:file:silent SWM/1/swm_rebuild.out > /dev/null 2> /dev/null 
