#!/bin/bash
#PBS -N _biox3_scratch_users_rhiju_projects_stepwise_benchmark_ref_rna_res_level_energy_rnatorsion1_2kcycles_FAVORITES2_j55a_P4P6_fixed_SWM_9
#PBS -o /dev/null
#PBS -e /dev/null
#PBS -m n
#PBS -M nobody@stanford.edu
#PBS -l walltime=48:00:00

cd /biox3/scratch/users/rhiju/projects/stepwise/benchmark/ref/rna_res_level_energy_rnatorsion1_2kcycles_FAVORITES2/j55a_P4P6_fixed

/home/rhiju/src/rosetta//main/source/bin/stepwise -s j55a_P4P6_fixed_START1_2r8s_RNA.pdb -fasta j55a_P4P6_fixed.fasta -terminal_res R:118 R:128 R:194 R:203 -extra_min_res R:121 R:127 R:195 R:200 -cycles 200 -nstruct 20 -intermolecular_frequency 0.0 -save_times -native j55a_P4P6_fixed_2r8s_RNA.pdb -score:rna_torsion_potential RNA11_based_new -chemical::enlarge_H_lj -set_weights rna_torsion 1.0 -allow_skip_bulge -cycles 2000 -out:file:silent SWM/9/swm_rebuild.out > /dev/null 2> /dev/null 
