#!/bin/bash

for FILE in `ls *.txt`; do ../scripts/python/benchmark_util/rewrite_stepwise_benchmark_info_file.py $FILE --tsv; done
