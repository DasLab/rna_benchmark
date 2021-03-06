#!/bin/bash
#PBS -N _biox3_scratch_users_rhiju_projects_stepwise_benchmark_rna_res_level_energy_rnatorsion1_favorites2_loopE_SWM_7
#PBS -o /dev/null
#PBS -e /dev/null
#PBS -m n
#PBS -M nobody@stanford.edu
#PBS -l walltime=16:00:00

cd /biox3/scratch/users/rhiju/projects/stepwise/benchmark/rna_res_level_energy_rnatorsion1_favorites2/loopE

/home/rhiju/src/rosetta//main/source/bin/stepwise -s loopE_HELIX1.pdb loopE_HELIX2.pdb -fasta loopE.fasta -terminal_res A:70 A:80 B:96 B:106 -extra_min_res A:71 A:79 B:97 B:105 -superimpose_over_all -cycles 200 -nstruct 20 -intermolecular_frequency 0.0 -save_times -native loopE_354d_RNA.pdb -score:rna_torsion_potential RNA11_based_new -chemical::enlarge_H_lj -set_weights rna_torsion 1.0 -cycles 200 -out:file:silent SWM/7/swm_rebuild.out > /dev/null 2> /dev/null 
