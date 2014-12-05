#!/bin/bash
#PBS -N _biox3_scratch_users_rhiju_projects_stepwise_benchmark_ref_rna_res_level_energy_rnatorsion1_2kcycles_FAVORITES2_tl_tr_P4P6_SWM_8
#PBS -o /dev/null
#PBS -e /dev/null
#PBS -m n
#PBS -M nobody@stanford.edu
#PBS -l walltime=48:00:00

cd /biox3/scratch/users/rhiju/projects/stepwise/benchmark/ref/rna_res_level_energy_rnatorsion1_2kcycles_FAVORITES2/tl_tr_P4P6

/home/rhiju/src/rosetta//main/source/bin/stepwise -s tl_tr_P4P6_HELIX3.pdb tl_tr_P4P6_START1_2r8s_RNA.pdb -fasta tl_tr_P4P6.fasta -terminal_res R:222 R:229 R:245 R:251 -extra_min_res R:223 R:227 R:247 R:250 -cycles 200 -nstruct 20 -intermolecular_frequency 0.0 -save_times -native tl_tr_P4P6_2r8s_RNA.pdb -score:rna_torsion_potential RNA11_based_new -chemical::enlarge_H_lj -set_weights rna_torsion 1.0 -allow_skip_bulge -cycles 2000 -out:file:silent SWM/8/swm_rebuild.out > /dev/null 2> /dev/null 
