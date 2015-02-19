stepwise_benchmark
==================

Python scripts &amp; reference data for stepwise modeling benchmarks

What to use this repo for
------------------------
Any new benchmarks that use the `stepwise` applications -- motif rebuilding, all tetraloops, miniproteins, helix nearest neighbor rules, etc. -- should go in here. Some or all of this repo will eventually move into the RosettaCommons scientific tests repo for regular testing.

What's in here
--------------
`ref/`  holds reference benchmark data

`new/` will hold new benchmark runs -- will not be checked in, however

`scripts/` holds some scripts for visualizing runs

`input_files/` holds input PDBs, and text files like `favorites.txt`, which define motifs for benchmark, including sequence, secondary structure, reference PDBs, etc.

Setting up benchmark tools  
--------------------------
- Make sure you have rosetta compiled and rna_tools setup. Follow directions <a href="https://www.rosettacommons.org/docs/latest/RNA-tools.html">here</a>.
- Edit the path to your copy of `stepwise_benchmark` in `INSTALL` and run `source ./INSTALL`

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
- If you are on `biox3`, you can type `source qsubMINI` to queue up all the jobs.
- While it is running or after it is done, run from command-line:
```
easy_cat.py SWM
```
to concatenate models from various subdirectories for each target.
- Copy or rsync files to your local computer to run MATLAB.

Visualizing benchmark runs
--------------------------
Plotting with MATLAB
- In the `scripts/matlab/` directory, you can compare two runs by running a command like:
```
make_plots( {'../new/rna_res_level_energy_rnatorsion1_novirtualo2prime_synGbonus_suitenessbonus','../new/rna_res_level_energy_rnatorsion1_novirtualo2prime_synGbonus_suitenessbonus_varypolarHgeom'} );
```
Plotting with Python
- In the `new/` directory, you can compare two runs by running a command like:
```
make_plots.py rna_res_level_energy4 rna_res_level_energy_rnatorsion1 
```
- In the `scripts/python/` directory, you can compare two runs by running a command like:
```
make_plots.py ../../new/rna_res_level_energy4 ../../new/rna_res_level_energy_rnatorsion1 
```
