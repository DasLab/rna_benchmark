#!/bin/bash
#PBS -N _biox3_scratch_users_rhiju_projects_stepwise_benchmark_rna_res_level_energy_rnatorsion1_2kcycles_srp_domainIV_fixed_SWM_1
#PBS -o /dev/null
#PBS -e /dev/null
#PBS -m n
#PBS -M nobody@stanford.edu
#PBS -l walltime=48:00:00

cd /biox3/scratch/users/rhiju/projects/stepwise/benchmark/rna_res_level_energy_rnatorsion1_2kcycles/srp_domainIV_fixed

/home/rhiju/src/rosetta//main/source/bin/stepwise -s srp_domainIV_fixed_START1_native_1lnt_RNA.pdb -fasta srp_domainIV_fixed.fasta -terminal_res A:3 A:10 B:15 B:22 -extra_min_res A:4 A:9 B:16 B:21 -cycles 200 -nstruct 20 -intermolecular_frequency 0.0 -native srp_domainIV_fixed_native_1lnt_RNA.pdb -score:rna_torsion_potential RNA11_based_new -chemical::enlarge_H_lj -set_weights rna_torsion 1.0 -cycles 2000 -out:file:silent SWM/1/swm_rebuild.out > /dev/null 2> /dev/null 
