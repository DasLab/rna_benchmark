#!/bin/bash
#PBS -N _Users_rhiju_Dropbox_projects_stepwise_benchmark_new_metal_ions_test_just_Arich_bulge_fixedMg_SWM_0
#PBS -o /dev/null
#PBS -e /dev/null
#PBS -m n
#PBS -M nobody@stanford.edu
#PBS -l walltime=48:00:00

cd /Users/rhiju/Dropbox/projects/stepwise/benchmark/new/metal_ions_test/just_Arich_bulge_fixedMg

/Users/rhiju/src/rosetta//main/source/bin/stepwise -s just_Arich_bulge_fixedMg_START1_2R8S.pdb -native just_Arich_bulge_fixedMg_NATIVE_2R8S.pdb -terminal_res R:133 R:137 R:181 R:189-190 -block_stack_above_res R:137 R:190 -block_stack_below_res R:133 R:181 R:189 -extra_min_res R:182 -fasta just_Arich_bulge_fixedMg.fasta -nstruct 20 -save_times -score:weights stepwise/rna/rna_res_level_energy4.wts -set_weights mg_lig 1.0 mg_sol 0.2 hoh_ref 1.0 -cycles 2000 -motif_mode -out:file:silent SWM/0/swm_rebuild.out > 0.out 2> 0.err 
