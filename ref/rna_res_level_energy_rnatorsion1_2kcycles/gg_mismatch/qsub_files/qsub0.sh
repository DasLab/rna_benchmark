#!/bin/bash
#PBS -N _biox3_scratch_users_rhiju_projects_stepwise_benchmark_rna_res_level_energy_rnatorsion1_2kcycles_gg_mismatch_SWM_0
#PBS -o /dev/null
#PBS -e /dev/null
#PBS -m n
#PBS -M nobody@stanford.edu
#PBS -l walltime=48:00:00

cd /biox3/scratch/users/rhiju/projects/stepwise/benchmark/rna_res_level_energy_rnatorsion1_2kcycles/gg_mismatch

/home/rhiju/src/rosetta//main/source/bin/stepwise -s gg_mismatch_HELIX1.pdb gg_mismatch_HELIX2.pdb -fasta gg_mismatch.fasta -terminal_res A:2 A:6 B:5 B:9 -extra_min_res A:3 A:5 B:6 B:8 -superimpose_over_all -cycles 200 -nstruct 20 -intermolecular_frequency 0.0 -native gg_mismatch_1f5g.pdb -score:rna_torsion_potential RNA11_based_new -chemical::enlarge_H_lj -set_weights rna_torsion 1.0 -cycles 2000 -out:file:silent SWM/0/swm_rebuild.out > /dev/null 2> /dev/null 
