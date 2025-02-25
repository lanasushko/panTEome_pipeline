#!/bin/bash

## This script reannotates the genomes with their own updated EDTA library

folder=$1
genome=$2
threads=$3

## 1 move around and change names to adapt for EDTA annotation run
cd $folder

name=$(basename $folder)
# rename original lib to save it as such
# mv $name'.Chr_scaffolds.fa.mod.EDTA.TElib.fa' $name'.Chr_scaffolds.fa.mod.EDTA.TElib.originalsgEDTA.fa'
# # rename updated lib for EDTA run
# mv $name'.Chr_scaffolds.fa.mod.EDTA.TElib.LTRupdated.fa.wgn.fa' $name'.Chr_scaffolds.fa.mod.EDTA.TElib.fa'
# # copy updated lib to final/ folder for EDTA run
# cp $name'.Chr_scaffolds.fa.mod.EDTA.TElib.fa' $name'.Chr_scaffolds.fa.mod.EDTA.final'

## 2 annotate genomes from their own EDTAlib
source activate EDTA2_fromyml
export LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libjemalloc.so.2
export BLASTDB_LMDB_MAP_SIZE=100000000

/ebio/abt6/ssushko/EDTA/EDTA.pl --genome $genome --overwrite 0 --step anno --threads $threads --anno 1

mv $name'.Chr_scaffolds.fa.mod.EDTA.TEanno.gff3' $name'.Chr_scaffolds.fa.mod.EDTA.TEanno.sgannot.gff3'

# Remove previous versions of file in case they exist

rm -f $name.Chr_scaffolds.fa.mod_divergence_plot.pdf_* 
rm -f $name.Chr_scaffolds.fa.mod.EDTA.intact.fa_*
rm -f $name.Chr_scaffolds.fa.mod.EDTA.intact.gff3_*
rm -f $name.Chr_scaffolds.fa.mod.EDTA.TEanno.gff3_*
rm -f $name.Chr_scaffolds.fa.mod.EDTA.TEanno.sum_*
rm -f $name.Chr_scaffolds.fa.mod.EDTA.TElib.fa_* 