#!/bin/bash

if [ $# -lt 2 ]; then
    echo "Must supply the name of a native pdb and target for easy_graft.sh to look in"
    exit 1
fi

NATIVE=$(readlink -ne $1)
STARTDIR=$(pwd)

WORKDIR=$2
if [ "$WORKDIR" == "." ]
then
	WORKDIR=$(pwd)
fi

MODELS=3
if [ $# -eq 3 ]
then
    MODELS=$3   
fi

# cd back into original directory
cd $WORKDIR

echo ">${NATIVE} ss_score=$(get_ss_score.py $NATIVE)"
echo $(get_sequence.py ${NATIVE})

for TARGET in */ ; do
	TARGET=${TARGET%*/}

	if [ "$TARGET" == "." ]
	then
		continue
	fi 

	if [ "$TARGET" == "models" ]
	then
		continue
	fi 
	easy_graft.sh $NATIVE $TARGET $MODELS
	#echo $TARGET
done

cd $STARTDIR

