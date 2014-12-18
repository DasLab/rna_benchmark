#!/bin/bash
#PBS -N _biox3_scratch_users_rhiju_projects_stepwise_benchmark_ref_rna_res_level_energy_rnatorsion1_2kcycles_FAVORITES2_puzzle1_alt_fixed_SWM_3
#PBS -o /dev/null
#PBS -e /dev/null
#PBS -m n
#PBS -M nobody@stanford.edu
#PBS -l walltime=48:00:00

cd /biox3/scratch/users/rhiju/projects/stepwise/benchmark/ref/rna_res_level_energy_rnatorsion1_2kcycles_FAVORITES2/puzzle1_alt_fixed

/home/rhiju/src/rosetta//main/source/bin/stepwise -s puzzle1_alt_fixed_START1_3mei_with_symm_RNA.pdb -fasta puzzle1_alt_fixed.fasta -terminal_res A:12 A:19 B:5 B:11 C:8 C:10 D:13 D:15 -extra_min_res A:13 A:18 B:6 B:10 -cycles 200 -nstruct 20 -intermolecular_frequency 0.0 -save_times -native puzzle1_alt_fixed_3mei_with_symm_RNA.pdb -score:rna_torsion_potential RNA11_based_new -chemical::enlarge_H_lj -set_weights rna_torsion 1.0 -allow_skip_bulge -cycles 2000 -out:file:silent SWM/3/swm_rebuild.out > /dev/null 2> /dev/null 
