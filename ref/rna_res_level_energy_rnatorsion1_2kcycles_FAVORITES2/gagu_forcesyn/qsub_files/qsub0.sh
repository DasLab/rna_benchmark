#!/bin/bash
#PBS -N _biox3_scratch_users_rhiju_projects_stepwise_benchmark_ref_rna_res_level_energy_rnatorsion1_2kcycles_FAVORITES2_gagu_forcesyn_SWM_0
#PBS -o /dev/null
#PBS -e /dev/null
#PBS -m n
#PBS -M nobody@stanford.edu
#PBS -l walltime=48:00:00

cd /biox3/scratch/users/rhiju/projects/stepwise/benchmark/ref/rna_res_level_energy_rnatorsion1_2kcycles_FAVORITES2/gagu_forcesyn

/home/rhiju/src/rosetta//main/source/bin/stepwise -s gagu_forcesyn_HELIX1.pdb gagu_forcesyn_HELIX2.pdb -fasta gagu_forcesyn.fasta -terminal_res A:2 A:9 B:13 B:20 -extra_min_res A:3 A:8 B:14 B:19 -superimpose_over_all -cycles 200 -nstruct 20 -intermolecular_frequency 0.0 -save_times -native gagu_forcesyn_2lx1_RNA.pdb -force_syn_chi_res_list A:4 A:6 B:15 B:17 -score:rna_torsion_potential RNA11_based_new -chemical::enlarge_H_lj -set_weights rna_torsion 1.0 -allow_skip_bulge -cycles 2000 -out:file:silent SWM/0/swm_rebuild.out > /dev/null 2> /dev/null 
