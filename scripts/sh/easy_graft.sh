#!/bin/bash

if [ $# -lt 2 ]; then
    echo "Must supply the name of a native pdb and target for easy_graft.sh to look in"
    exit 1
fi

NATIVE=$(readlink -ne $1)
TARGET=${2%*/}
MODELS=3

if [ $# -eq 3 ]
then
    MODELS=$3   
fi

WORKDIR=$(pwd)

cd $TARGET


# Extract models from silent file into a model directory

if [ -d "models" ]; then
  rm -r models
fi

mkdir models 
cp swm_rebuild.out models/${TARGET}.out
cd models
out=$(extract_lowscore_decoys.py ${TARGET}.out $MODELS)

#echo ">${NATIVE}"
#echo $(get_sequence.py ${NATIVE})


# Graft into Native
for pdb in ${TARGET}.out.*.pdb ; do

	#compare length of native sequence with working native pdb sequence
	working_native="../${TARGET}_$(basename $NATIVE)"
	native_sequence=$(get_sequence.py $working_native)
	design_sequence=$(get_sequence.py $pdb)

	native_size=${#native_sequence}
	size=${#design_sequence}
	#echo "Native: $native_sequence Size: $native_size"
	#echo "Design: $design_sequence Size: $size"


	grafted_pdb=${pdb/.out./.out.grafted.}
	
	out=$(
		rna_graft.linuxgccrelease \
		-s $NATIVE $pdb \
		-superimpose_res $(get_res_num.py ../*START1*.pdb) \
		-o $grafted_pdb
	) 

	echo ">${grafted_pdb}  ss_score=$(get_ss_score.py $grafted_pdb) missing=$(($native_size - $size))"
	echo $(get_sequence.py $grafted_pdb)

done

# cd back into original directory
cd $WORKDIR