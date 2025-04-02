#!/bin/bash

# This script takes panTElib and genome with previous single-genome EDTA annotation as input and annotates the genome using families from panTElib

# Input
genome=$1
panTElib=$2
workdir=$3
threads=$4

# Internal variables
name=$(basename $genome | cut -f1 -d'.')
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

# load EDTA prerequisites

export LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libjemalloc.so.2
export BLASTDB_LMDB_MAP_SIZE=100000000

# run EDTA reannotation
cd $workdir
perl $SCRIPT_DIR/EDTA/EDTA.pl --genome $genome --threads $threads --step final --anno 1 --curatedlib $panTElib

# rename for pipeline streamlinness
mv $name.Chr_scaffolds.fa.mod_divergence_plot.pdf $name.Chr_scaffolds.fa.mod_divergence_plot.panTEannot.pdf
mv $name.Chr_scaffolds.fa.mod.EDTA.TEanno.density_plots.pdf $name.Chr_scaffolds.fa.TEanno.density_plots.panTEannot.pdf
mv $name.Chr_scaffolds.fa.mod.EDTA.intact.fa $name.Chr_scaffolds.fa.mod.EDTA.intact.panTEannot.fa
mv $name.Chr_scaffolds.fa.mod.EDTA.intact.gff3 $name.Chr_scaffolds.fa.mod.EDTA.intact.panTEannot.gff3
mv $name.Chr_scaffolds.fa.mod.EDTA.TEanno.gff3 $name.Chr_scaffolds.fa.mod.EDTA.TEanno.panTEannot.gff3
mv $name.Chr_scaffolds.fa.mod.EDTA.TEanno.sum $name.Chr_scaffolds.fa.mod.EDTA.TEanno.panTEannot.sum
mv $name.Chr_scaffolds.fa.mod.EDTA.TElib.fa $name.Chr_scaffolds.fa.mod.EDTA.TElib.panTEannot.fa

mv $name.Chr_scaffolds.fa.mod_divergence_plot.pdf_* $name.Chr_scaffolds.fa.mod_divergence_plot.pdf
mv $name.Chr_scaffolds.fa.mod.EDTA.TEanno.density_plots.pdf_* $name.Chr_scaffolds.fa.TEanno.density_plots.pdf
mv $name.Chr_scaffolds.fa.mod.EDTA.intact.fa_* $name.Chr_scaffolds.fa.mod.EDTA.intact.fa
mv $name.Chr_scaffolds.fa.mod.EDTA.intact.gff3_* $name.Chr_scaffolds.fa.mod.EDTA.intact.gff3
mv $name.Chr_scaffolds.fa.mod.EDTA.TEanno.gff3_* $name.Chr_scaffolds.fa.mod.EDTA.TEanno.gff3
mv $name.Chr_scaffolds.fa.mod.EDTA.TEanno.sum_* $name.Chr_scaffolds.fa.mod.EDTA.TEanno.sum
mv $name.Chr_scaffolds.fa.mod.EDTA.TElib.fa_* $name.Chr_scaffolds.fa.mod.EDTA.TElib.fa
