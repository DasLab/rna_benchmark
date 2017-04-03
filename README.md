Management and analysis of benchmarks of RNA motif structure prediction
==================

The ab initio prediction of the structure of RNA loops is an excellent stress test for RNA structure prediction algorithms. This benchmark was originally developed for the stepwise Monte Carlo (SWM) and stepwise assembly (SWA) algorithms, and since then it has developed to permit the generation of equivalent "fragment assembly of RNA" (with optional "full-atom refinement"; FARNA and FARFAR, respectively) jobs. 

As a result, this benchmark may be used to gauge the performance, on common problems of broad interest, of new algorithms or new configurations thereof, or to help optimize new scoring terms or to re-weight existing scoring functions. 

This repository includes python scripts for constructing benchmark sets of loops from text specifications, setting up runs (either locally or on a cluster) for a particular benchmark set, and analyzing the results (via plots or tables, both in the styles employed in a publication currently under review).

How to contribute to this repository
------------------------------------
For the most part, the structure and contents of this repository is likely to remain fixed for some time: it is intended to serve as a reproducible starting point for the results obtained in the aforementioned paper. Afterwards, however, the SHA1 of that starting point will be memorialized here and further contributions will be encouraged.

What might be a valuable contribution? Maybe you've got a set of structures that you think stepwise ought to be able to recover well (or poorly); maybe you've got a fun (preferably: informative) visualization or ensemble analysis method you'd like to add to our analysis and plotting scripts.

Any new benchmarks that use the `stepwise` applications -- motif rebuilding, all tetraloops, miniproteins, helix nearest neighbor rules, etc. -- should go in here. Some or all of this repo may eventually be used for Rosetta scientific tests.

What's in here
--------------
`ref/`  holds reference benchmark data

`new/` will hold new benchmark runs -- will not be checked in, however

`scripts/` holds some scripts for visualizing runs

`input_files/` holds input PDBs, and text files like `favorites.txt`, which define motifs for benchmark, including sequence, secondary structure, reference PDBs, etc.

Setting up benchmark tools  
--------------------------
- Make sure you have Rosetta compiled. You also need the Rosetta tools repository, and to have set up the rna_tools subdirectory therein. Follow directions <a href="https://www.rosettacommons.org/docs/latest/RNA-tools.html">here</a>.
- Run `source ./INSTALL` from inside the `stepwise_benchmark/` directory.

Setting up a new benchmark run
------------------------------
(Look  inside `ref/` for checked in examples)
- Go inside `new/`
- Create a new folder with a descriptive name like `rna_res_level_energy_rnatorsion1_synGbonus`
- Create `extra_flags_benchmark.txt` that describes the interesting new flags in your run.
- Create a `README_SETUP` with a command-line to setup the benchmark like:
```
setup_stepwise_benchmark.py ../../input_files/<favorites.txt>
```
 and run it with `source README_SETUP`. You should see subdirectories with all the target names and input files.
- If you are on Stanford's `biox3`, or another cluster running PBS, you can type `source qsubMINI` to queue up all the jobs. Ditto for if you are on Stanford's `sherlock`, or another cluster running slurm; you can then type `source sbatchMINI`.
- While it is running or after it is done, run from the command-line:
```
easy_cat.py SWM
```
to concatenate models from various subdirectories for each target.
- Copy or rsync files to your local computer for visualization.

Regression tests
----------------

You should run regression tests before making any merge to the `master` branch. (Once we begin accepting contributions, we will be the ones making merges to `master` -- but you should explicitly list any regression test output in your PR.) Since these regression tests are VERY FAST (~15s), ideally you should run these regression tests more often than that.
- Go to the `test/` directory.
- Run `test.sh` on master to generate reference results (in `ref``).
- Check out your new branch, commit of interest, or whatever.
- Re-run `test.sh` to generate new results (in `new/`). The script will report on the diff automatically.

Visualizing benchmark runs
--------------------------
Plotting with Python
- In the `new/` directory, you can compare two runs by running a command like:
```
make_plots.py rna_res_level_energy4 rna_res_level_energy_rnatorsion1 
```
Plotting with MATLAB
- In the `scripts/matlab/` directory, you can compare two runs by running a command like:
```
make_plots( {'../new/rna_res_level_energy_rnatorsion1_novirtualo2prime_synGbonus_suitenessbonus','../new/rna_res_level_energy_rnatorsion1_novirtualo2prime_synGbonus_suitenessbonus_varypolarHgeom'} );
```
