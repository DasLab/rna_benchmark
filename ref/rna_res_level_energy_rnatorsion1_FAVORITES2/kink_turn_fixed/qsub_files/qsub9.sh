#!/bin/bash
#PBS -N _biox3_scratch_users_rhiju_projects_stepwise_benchmark_rna_res_level_energy_rnatorsion1_favorites2_kink_turn_fixed_SWM_9
#PBS -o /dev/null
#PBS -e /dev/null
#PBS -m n
#PBS -M nobody@stanford.edu
#PBS -l walltime=16:00:00

cd /biox3/scratch/users/rhiju/projects/stepwise/benchmark/rna_res_level_energy_rnatorsion1_favorites2/kink_turn_fixed

/home/rhiju/src/rosetta//main/source/bin/stepwise -s kink_turn_fixed_START1_2gis_RNA.pdb -fasta kink_turn_fixed.fasta -terminal_res A:16 A:22 30 39 -extra_min_res A:17 A:21 31 38 -cycles 200 -nstruct 20 -intermolecular_frequency 0.0 -save_times -native kink_turn_fixed_2gis_RNA.pdb -score:rna_torsion_potential RNA11_based_new -chemical::enlarge_H_lj -set_weights rna_torsion 1.0 -cycles 200 -out:file:silent SWM/9/swm_rebuild.out > /dev/null 2> /dev/null 
