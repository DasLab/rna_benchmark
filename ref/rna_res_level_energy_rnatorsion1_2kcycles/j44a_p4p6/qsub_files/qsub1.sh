#!/bin/bash
#PBS -N _biox3_scratch_users_rhiju_projects_stepwise_benchmark_rna_res_level_energy_rnatorsion1_2kcycles_j44a_p4p6_SWM_1
#PBS -o /dev/null
#PBS -e /dev/null
#PBS -m n
#PBS -M nobody@stanford.edu
#PBS -l walltime=48:00:00

cd /biox3/scratch/users/rhiju/projects/stepwise/benchmark/rna_res_level_energy_rnatorsion1_2kcycles/j44a_p4p6

/home/rhiju/src/rosetta//main/source/bin/stepwise -s j44a_p4p6_HELIX1.pdb j44a_p4p6_HELIX2.pdb -fasta j44a_p4p6.fasta -terminal_res A:111 A:118 A:203 A:209 -extra_min_res A:112 A:116 A:205 A:208 -superimpose_over_all -cycles 200 -nstruct 20 -intermolecular_frequency 0.0 -native j44a_p4p6_1gid_RNAA.pdb -score:rna_torsion_potential RNA11_based_new -chemical::enlarge_H_lj -set_weights rna_torsion 1.0 -cycles 2000 -out:file:silent SWM/1/swm_rebuild.out > /dev/null 2> /dev/null 
