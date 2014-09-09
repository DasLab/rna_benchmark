#!/bin/bash
#PBS -N _biox3_scratch_users_rhiju_projects_stepwise_benchmark_rna_res_level_energy_rnatorsion1_2kcycles_gcaa_tetraloop_SWM_2
#PBS -o /dev/null
#PBS -e /dev/null
#PBS -m n
#PBS -M nobody@stanford.edu
#PBS -l walltime=48:00:00

cd /biox3/scratch/users/rhiju/projects/stepwise/benchmark/rna_res_level_energy_rnatorsion1_2kcycles/gcaa_tetraloop

/home/rhiju/src/rosetta//main/source/bin/stepwise -s gcaa_tetraloop_HELIX1.pdb -fasta gcaa_tetraloop.fasta -terminal_res A:3 A:10 -extra_min_res A:4 A:9 -superimpose_over_all -cycles 200 -nstruct 20 -intermolecular_frequency 0.0 -native gcaa_tetraloop_1zih_RNA.pdb -score:rna_torsion_potential RNA11_based_new -chemical::enlarge_H_lj -set_weights rna_torsion 1.0 -cycles 2000 -out:file:silent SWM/2/swm_rebuild.out > /dev/null 2> /dev/null 
