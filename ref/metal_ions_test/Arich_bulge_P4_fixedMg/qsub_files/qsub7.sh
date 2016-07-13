#!/bin/bash
#PBS -N _Users_rhiju_Dropbox_projects_stepwise_benchmark_new_metal_ions_test_Arich_bulge_P4_fixedMg_SWM_7
#PBS -o /dev/null
#PBS -e /dev/null
#PBS -m n
#PBS -M nobody@stanford.edu
#PBS -l walltime=48:00:00

cd /Users/rhiju/Dropbox/projects/stepwise/benchmark/new/metal_ions_test/Arich_bulge_P4_fixedMg

/Users/rhiju/src/rosetta//main/source/bin/stepwise -s Arich_bulge_P4_fixedMg_START1_2R8S.pdb -native Arich_bulge_P4_fixedMg_NATIVE_2R8S.pdb -terminal_res R:108 R:111 R:133 R:137 R:181 R:189-190 R:210 R:213 -block_stack_above_res R:111 R:137 R:190 R:213 -block_stack_below_res R:108 R:133 R:181 R:189 R:210 -extra_min_res R:182 -fasta Arich_bulge_P4_fixedMg.fasta -nstruct 20 -save_times -score:weights stepwise/rna/rna_res_level_energy4.wts -set_weights mg_lig 1.0 mg_sol 0.2 hoh_ref 1.0 -cycles 2000 -motif_mode -out:file:silent SWM/7/swm_rebuild.out > 7.out 2> 7.err 
