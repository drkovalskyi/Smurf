#!/bin/env perl
use strict;
use warnings;
die "Usage:\n\t$0 <input file selection> <output directory>\n" unless @ARGV==2;
my $pattern = $ARGV[0];
$pattern .= "/*.root" if (-d $pattern);
foreach my $file (`ls -1 $pattern|grep \.root`){
    $file =~ s/\n//;
    my $filename = $file;
    $filename =~ s/^(.*?\/)([^\/]+)$/$2/;
    printf("skimming $file to $ARGV[1]/$filename\n");
    system("./skimtree $file $ARGV[1]/$filename");
}
exit
