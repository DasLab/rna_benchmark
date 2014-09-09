#!/bin/bash
#PBS -N _biox3_scratch_users_rhiju_projects_stepwise_benchmark_rna_res_level_energy_rnatorsion1_gagua_pentaloop_SWM_7
#PBS -o /dev/null
#PBS -e /dev/null
#PBS -m n
#PBS -M nobody@stanford.edu
#PBS -l walltime=16:00:00

cd /biox3/scratch/users/rhiju/projects/stepwise/benchmark/rna_res_level_energy_rnatorsion1/gagua_pentaloop

/home/rhiju/src/rosetta//main/source/bin/stepwise -s gagua_pentaloop_HELIX1.pdb -fasta gagua_pentaloop.fasta -terminal_res A:20 A:28 -extra_min_res A:21 A:27 -superimpose_over_all -cycles 200 -nstruct 20 -intermolecular_frequency 0.0 -save_times -native gagua_pentaloop_1xjr_RNA.pdb -score:rna_torsion_potential RNA11_based_new -chemical::enlarge_H_lj -set_weights rna_torsion 1.0 -cycles 200 -out:file:silent SWM/7/swm_rebuild.out > /dev/null 2> /dev/null 
