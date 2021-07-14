#!/usr/bin/env perl
use strict;

@ARGV or die "Please provide a directory name containing the fastq format sequence files!\n";
my $classifer_output = $ARGV[0];

 my %tbhash = (
     "Mycobacterium tuberculosis complex" =>1,
     "Mycobacterium tuberculosis" => 1, 
     "Mycobacterium africanum" =>1, 
     "Mycobacterium bovis" =>1, 
     "Mycobacterium microti" =>1,
     "Mycobacterium canettii" =>1, 
     "Mycobacterium caprae" =>1,
     "Mycobacterium pinnipedii" =>1,
     "Mycobacterium suricattae" =>1,
     "Mycobacterium mungi" =>1,
     "Mycobacterium dassie" =>1,
     "Mycobacterium orygis" =>1
 );

 my %seqs = ();
open my $taxon_classifier, '<', $classifer_output;
while(<$taxon_classifier>){
    next if /^U/;
    my @l = split(/\t/, $_);
    my $taxon = $l[2];
    $taxon =~ s/\s+\(.*//;
    print STDERR "*******$taxon\n";
    if(exists $tbhash{$taxon}){
        $seqs{$l[1]} = 1;
        print "$l[1]\n";
    }
    

}

