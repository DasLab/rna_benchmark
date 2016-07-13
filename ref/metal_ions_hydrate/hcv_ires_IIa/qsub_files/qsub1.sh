#!/bin/bash
#PBS -N _biox3_scratch_users_rhiju_projects_stepwise_benchmark_new_metal_ions_hydrate_hcv_ires_IIa_SWM_1
#PBS -o /dev/null
#PBS -e /dev/null
#PBS -m n
#PBS -M nobody@stanford.edu
#PBS -l walltime=48:00:00

cd /biox3/scratch/users/rhiju/projects/stepwise/benchmark/new/metal_ions_hydrate/hcv_ires_IIa

/home/rhiju/src/rosetta//main/source/bin/stepwise -s hcv_ires_IIa_START1_2PN3.pdb -native hcv_ires_IIa_NATIVE_2PN3.pdb -terminal_res A:51 A:59 B:109 B:112 -block_stack_above_res A:59 B:112 -block_stack_below_res A:51 B:109 -extra_min_res A:52 A:58 -fasta hcv_ires_IIa.fasta -nstruct 20 -save_times -score:weights stepwise/rna/rna_res_level_energy4.wts -set_weights mg_lig 1.0 mg_sol 0.2 hoh_ref 1.0 -cycles 2000 -magnesium:hydrate -motif_mode -out:file:silent SWM/1/swm_rebuild.out > 1.out 2> 1.err 
