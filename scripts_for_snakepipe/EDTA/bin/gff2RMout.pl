#!/usr/bin/env perl
use strict;
use warnings;

# Shujun Ou (shujun.ou.1@gmail.com)
# Facilitated by ChatGPT
#
my $usage = "Usage: perl gff2RMout.pl genome.gff3 genome.out";
die "\n$usage\n\n" unless $ARGV[1];

# Open the GFF3 file for reading
open(my $gff3_file, "<", $ARGV[0]) or die "Cannot open input.gff3: $!";

# Open the RepeatMasker .out file for writing
open(my $rm_out_file, ">", $ARGV[1]) or die "Cannot open output.out: $!";

# Read and process the GFF3 file line by line
while (my $line = <$gff3_file>) {
    chomp $line;
    next if $line =~ /^#/;
    my @fields = split /\t/, $line;

    # Extract necessary information (adjust indices based on your GFF3 structure)
    my $seqName = $fields[0];
    my $source = $fields[1];
    my $featureType = $fields[2];
    my $start = $fields[3];
    my $end = $fields[4];
    my $score = $fields[5];
    my $strand = $fields[6];
    my $phase = $fields[7];
    my $attributes = $fields[8];
    my $family = $1 if $attributes =~ /Name=([^\s;]+);/;
    my $classification = $1 if $attributes =~ /classification=([^\s;]+);?/;

    # Format the output for RepeatMasker .out (modify according to the expected format)
    printf $rm_out_file "10000 0.001 0.001 0.001\t$seqName\t$start\t$end\tNA\t$strand\t$family\t$classification\tNA NA NA NA\n";
}

# Close the filehandles
close($gff3_file);
close($rm_out_file);

