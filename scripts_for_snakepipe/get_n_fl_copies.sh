#!/bin/bash

## This script extracts the TE families with at least $fl_copy number of copies per genome

RMout=$1
fl_copy=$2
find_flTE=$3


workdir=$(echo $RMout | rev | cut -d'/' -f3- | rev)
name=$(basename $workdir)

## get fl-TE with ≥ $fl_copy copies in each genome

perl $find_flTE $RMout | awk '{print $10}' | sort | uniq -c | perl -snale 'my ($count, $id) = (split); next if $count < $fl_copy; print $_' -- -fl_copy=$fl_copy | awk '{print $2"#"}' > $workdir/$name.Chr_scaffolds.fa.mod.EDTA.TElib.fa.keep.list

## extract pan-TE library candidate sequences [CHANGE THIS TO SEQKIT]
sed -i 's/#/#\.\*/g' $workdir/$name.Chr_scaffolds.fa.mod.EDTA.TElib.fa.keep.list # adapt to regex pattern

seqkit grep -n -r -f $workdir/$name.Chr_scaffolds.fa.mod.EDTA.TElib.fa.keep.list $workdir/$name.Chr_scaffolds.fa.mod.EDTA.TElib.fa -w 0 -o $workdir/$name.Chr_scaffolds.fa.mod.EDTA.TElib.mincopyfiltered.fa # extract sequences with seqkit

