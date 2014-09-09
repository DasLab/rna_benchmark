#!/bin/bash
#PBS -N _biox3_scratch_users_rhiju_projects_stepwise_benchmark_rna_res_level_energy_rnatorsion1_srl_free_SWM_5
#PBS -o /dev/null
#PBS -e /dev/null
#PBS -m n
#PBS -M nobody@stanford.edu
#PBS -l walltime=16:00:00

cd /biox3/scratch/users/rhiju/projects/stepwise/benchmark/rna_res_level_energy_rnatorsion1/srl_free

/home/rhiju/src/rosetta//main/source/bin/stepwise -s srl_free_HELIX1.pdb srl_free_HELIX2.pdb -fasta srl_free.fasta -terminal_res A:651 A:669 -extra_min_res A:652 A:668 -superimpose_over_all -cycles 200 -nstruct 20 -intermolecular_frequency 0.0 -save_times -native srl_free_1q9a_RNA.pdb -score:rna_torsion_potential RNA11_based_new -chemical::enlarge_H_lj -set_weights rna_torsion 1.0 -cycles 200 -out:file:silent SWM/5/swm_rebuild.out > /dev/null 2> /dev/null 
