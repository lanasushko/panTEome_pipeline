#!/bin/bash

# Input
edtabed=$1
edtagff=$2
trashgff=$3

# Internal variables
name_bed=$(basename $edtabed)
name_gff=$(basename $edtagff)

# Subtract repeats and remove FPs (Penelope)
bedtools subtract -A -a $edtabed -b $trashgff | grep -v 'Penelope' > $name_bed.cleaned.bed
bedtools subtract -A -a $edtagff -b $trashgff | grep -v 'Penelope' > $name_gff.cleaned.gff

# Reformat .bed and add standarized labels to .gff


# Merge 