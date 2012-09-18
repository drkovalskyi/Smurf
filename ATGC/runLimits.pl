#!/bin/env perl
use strict;
use warnings;
die "Usage:\n\t$0 <dir>\n\n" unless @ARGV==1;
my @cards = ();
foreach my $card(`cd $ARGV[0]; ls -1 atgc_*.card`){
    $card =~ s/\n//;
    my $name = $card;
    $name =~ s/\.card$//;
    print "Running card $ARGV[0]/$card\n";
    system("cd $ARGV[0]; combine -M Asymptotic $card -n $name --rMax 2 &> $name.log");
}
