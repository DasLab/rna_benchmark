#!/bin/bash

# This script is useful for running easy_cat.py in specified directories only.

if [ $# -eq 0 ]; then
    echo "Must supply the name of a directory for easy_cat.py to look in"
    exit 1
fi

if [ $# -eq 2 ]
then
    DIRS=$2*   
else
    DIRS=*
fi

echo "Running 'easy_cat.py $1' in the following directories: "

for DIR in $DIRS; do
    echo $DIR
done
 
for DIR in $DIRS; do
    cd $DIR
    pwd
    easy_cat.py $1
    cd ..
done