if [[ -d ref ]] ; then
	if [[ -d new ]]; then
		rm -rf new
	fi
	
	mkdir new
	cd new
	for benchmark in RNA_loop_motifs_PS2011 favorites favorites2 challenges ; do
	    echo "Running setup for: $benchmark.txt"
	    ../../scripts/python/benchmark_util/setup_stepwise_benchmark.py ../../input_files/$benchmark.txt >/dev/null
	done
	cd ..
	diff -r -u -I '#PBS.*' -I '/Users/amw579/nucleic/stepwise_benchmark/test' -x '*.sh' -x 'qsubMINI' -x 'bsubMINI' -x 'sbatchMINI' -x 'sbatchMPI' -x 'condorMINI' -x 'MPI_ONEBATCH.job' -x 'sbatch_files' -x 'sbatch_files_MPI' -x 'sbatchMPI_ONEBATCH' -x 'README_SWM' --exclude="*/qsub_files" ref new > tempdiff
	
	if [[ -s tempdiff ]]; then
	    if type "colordiff" >& /dev/null; then
		cat tempdiff | colordiff
	    else
		cat tempdiff
		echo "(If you install colordiff, the diffs above will look pretty.)"
	    fi
	    echo "Directories were meaningfully different; do not commit these changes until you are certain the diff is correct."
	else
	    echo "Clean diff; go ahead and commit."
	fi
else
	# generate ref
	mkdir ref
	cd ref
	for benchmark in RNA_loop_motifs_PS2011 favorites favorites2 challenges ; do
		../../scripts/python/benchmark_util/setup_stepwise_benchmark.py ../../input_files/$benchmark.txt
	done
	cd ..
fi


