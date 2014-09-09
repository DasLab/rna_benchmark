#!/bin/bash
#PBS -N _biox3_scratch_users_rhiju_projects_stepwise_benchmark_ref_rna_res_level_energy_rnatorsion1_2kcycles_FAVORITES2_kiss_add_L3_fixed_SWM_6
#PBS -o /dev/null
#PBS -e /dev/null
#PBS -m n
#PBS -M nobody@stanford.edu
#PBS -l mem=500Mb
#PBS -l walltime=48:00:00

cd /biox3/scratch/users/rhiju/projects/stepwise/benchmark/ref/rna_res_level_energy_rnatorsion1_2kcycles_FAVORITES2/kiss_add_L3_fixed

/home/rhiju/src/rosetta//main/source/bin/stepwise -s kiss_add_L3_fixed_START1_1y26_RNA.pdb -fasta kiss_add_L3_fixed.fasta -terminal_res X:29 X:41 X:58 X:68 -extra_min_res X:59 X:67 -cycles 200 -nstruct 20 -intermolecular_frequency 0.0 -save_times -native kiss_add_L3_fixed_1y26_RNA.pdb -score:rna_torsion_potential RNA11_based_new -chemical::enlarge_H_lj -set_weights rna_torsion 1.0 -cycles 2000 -out:file:silent SWM/6/swm_rebuild.out > /dev/null 2> /dev/null 
