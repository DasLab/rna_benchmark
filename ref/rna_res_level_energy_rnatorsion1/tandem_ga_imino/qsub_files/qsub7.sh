#!/bin/bash
#PBS -N _biox3_scratch_users_rhiju_projects_stepwise_benchmark_rna_res_level_energy_rnatorsion1_tandem_ga_imino_SWM_7
#PBS -o /dev/null
#PBS -e /dev/null
#PBS -m n
#PBS -M nobody@stanford.edu
#PBS -l walltime=16:00:00

cd /biox3/scratch/users/rhiju/projects/stepwise/benchmark/rna_res_level_energy_rnatorsion1/tandem_ga_imino

/home/rhiju/src/rosetta//main/source/bin/stepwise -s tandem_ga_imino_HELIX1.pdb tandem_ga_imino_HELIX2.pdb -fasta tandem_ga_imino.fasta -terminal_res A:2 A:7 B:10 B:15 -extra_min_res A:3 A:6 B:11 B:14 -superimpose_over_all -cycles 200 -nstruct 20 -intermolecular_frequency 0.0 -save_times -native tandem_ga_imino_1mis_RNA.pdb -score:rna_torsion_potential RNA11_based_new -chemical::enlarge_H_lj -set_weights rna_torsion 1.0 -cycles 200 -out:file:silent SWM/7/swm_rebuild.out > /dev/null 2> /dev/null 
