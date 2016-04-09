#######################################################
# Nearest neighbor energy-scale checks on Rosetta RNA #
#######################################################
# R. Das, April 2016, Updating README from 2013.

# Create a bunch of models using rna_helix.py but 'finishing' (final minimization)
#  with your new weight set.
# Models are mainly helices, but can also generate models with dangles.
# A simple MATLAB scripts compares energies of these models with standard
# Turner parameters (37C, 1M NaCl), and also suggests how to refit.

#######################################################
# Make the models

#go into helix/
#
#run:
 
README_HELIX

# This is a wrapper around generate_helix.py, and
# creates a bunch of duplexes, duplexes with dangling ends, etc.
#
# Note that this can take a while. For just canonical base pair step nearest neighbor rules
#  (no dangles), can get rid of the code block in generate_helix.py that
#  is under 'singlets with bulge', and all the seq_sets.append(...) code for special cases.
#
# Can also save time in first runs by getting rid of gu and ug pairs (edit line beginning with "bps = ['au',..  ")
#
# Last, currently the routine is hard coded to assume you have 4 cores at your disposal... probably
#  should make this a flag. 
#
# Models show up in pdb/ file. Old runs are in pdb_TRY18, etc. ...
# To see what results they gave in matlab analysis (see next), can copy their results into pdb.

#######################################################
# launch MATLAB, and change directory into matlab/.
#  
# To extracts NN parameters by subtracting energies of different models in pdb/ directory

fit_weights

# Computes energies of each model from NN rules and compares to rosetta energies:
fit_weights_to_full_model_energies

# Both the fits above are least-squares fits, and include some regularization to make sure new weights 
#  are not radically different from old weights. Otherwise, you get crazy new weights (negative numbers) and
#  garbage upon reminimizing models.

#
# TODO: this analysis is probably better encoded in python+numpy. Would take a half hour to convert the code --
#   Kalli, perhaps we can do this quickly?

#######################################################
# Example output running fit_weights_to_full_model_energies
# on my most recent version of Rosetta (Jan 2016):

Multiplying input scores by 0.616000!

WEIGHTS
    1.02  fa_atr
    0.85  fa_rep
    0.88  fa_sol
    1.00  fa_intra_rep
    0.97  fa_elec_rna_phos_phos
    0.99  rna_torsion
    1.00  suiteness_bonus
    0.93  rna_sugar_close
    1.13  fa_stack
    0.93  stack_elec
    0.65  geom_sol_fast
    0.95  hbond
    0.63  free_base
    0.10  free_dof
    1.00  intermol
    1.00  other_pose
    1.00  loop_close
    1.00  linear_chainbreak
    1.05  unfolded_g
    1.14  unfolded_a
    1.16  unfolded_c
    1.27  unfolded_u

fit   rmsd:    0.564
raw   rmsd:    1.273

# Example output image is in 

# Note that in above, the scores are multiplied by 0.616 to convert Rosetta units (assumed to be close to kB T with T = 37 C) 
# to kcal/mol


# Weights are actually scalefactors to apply to individual weight components. If you apply the numbers above, then need to rerun
# modeling afterwards to reminimizine, and then iterate a bit.

######################################################################
# TODO: Fang's RECCES framework is way better than this, but is appropriate only if 
#  you use the energy function in the context of modeling ensembles. This code
#  allows checks on the effective energy function that is minimized in stepwise assembly or fragment assembly.
  
