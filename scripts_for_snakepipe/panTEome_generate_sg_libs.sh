#!/bin/bash

source activate EDTA2_fromyml
export LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libjemalloc.so.2
export BLASTDB_LMDB_MAP_SIZE=100000000


genome=$1
threads=$2
directory=$3

mkdir $directory
cd $directory
/ebio/abt6/ssushko/EDTA/EDTA.pl --genome $genome --overwrite 0 --anno 0 --threads $threads
