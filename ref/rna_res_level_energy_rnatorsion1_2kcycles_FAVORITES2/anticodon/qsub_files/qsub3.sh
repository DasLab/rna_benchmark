#!/bin/bash
#PBS -N _biox3_scratch_users_rhiju_projects_stepwise_benchmark_ref_rna_res_level_energy_rnatorsion1_2kcycles_FAVORITES2_anticodon_SWM_3
#PBS -o /dev/null
#PBS -e /dev/null
#PBS -m n
#PBS -M nobody@stanford.edu
#PBS -l mem=500Mb
#PBS -l walltime=48:00:00

cd /biox3/scratch/users/rhiju/projects/stepwise/benchmark/ref/rna_res_level_energy_rnatorsion1_2kcycles_FAVORITES2/anticodon

/home/rhiju/src/rosetta//main/source/bin/stepwise -s anticodon_HELIX1.pdb -fasta anticodon.fasta -terminal_res A:29 A:41 -extra_min_res A:31 A:39 -superimpose_over_all -cycles 200 -nstruct 20 -intermolecular_frequency 0.0 -save_times -native anticodon_3l0u_RNA.pdb -score:rna_torsion_potential RNA11_based_new -chemical::enlarge_H_lj -set_weights rna_torsion 1.0 -cycles 2000 -out:file:silent SWM/3/swm_rebuild.out > /dev/null 2> /dev/null 
