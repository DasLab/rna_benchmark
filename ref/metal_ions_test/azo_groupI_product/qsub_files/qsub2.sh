#!/bin/bash
#PBS -N _Users_rhiju_Dropbox_projects_stepwise_benchmark_new_metal_ions_test_azo_groupI_product_SWM_2
#PBS -o /dev/null
#PBS -e /dev/null
#PBS -m n
#PBS -M nobody@stanford.edu
#PBS -l walltime=48:00:00

cd /Users/rhiju/Dropbox/projects/stepwise/benchmark/new/metal_ions_test/azo_groupI_product

/Users/rhiju/src/rosetta//main/source/bin/stepwise -s azo_groupI_product_START1_3BO3_altA_renumber.pdb -native azo_groupI_product_NATIVE_3BO3_altA_renumber.pdb -terminal_res B:8 B:12 B:127 B:137 B:169 B:179 D:1 -block_stack_above_res B:12 B:137 B:179 -block_stack_below_res B:8 B:127 B:169 D:1 -extra_min_res D:2 -fasta azo_groupI_product.fasta -nstruct 20 -save_times -score:weights stepwise/rna/rna_res_level_energy4.wts -set_weights mg_lig 1.0 mg_sol 0.2 hoh_ref 1.0 -cycles 2000 -motif_mode -out:file:silent SWM/2/swm_rebuild.out > 2.out 2> 2.err 
