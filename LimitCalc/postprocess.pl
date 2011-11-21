#!/bin/env perl
use strict;
use warnings;
my $usage = "Usage:\n\t$0 <path>/<file with a list of cards> <n-jets title> <type>\n\nValid types are: Bayesian, CLs, CLs-asymptotic\n\n";
die $usage unless @ARGV == 3;
my ($dir,$name) = ($ARGV[0] =~ /^(.+?)\/([^\/]+)$/);
die $usage if (!defined $dir || !defined $name); 
$name =~ s/.*?([^\/]+)$/$1/;
$name =~ s/\.[^\.]+$//;
my $type = $ARGV[2];
my $title = "H #rightarrow WW #rightarrow 2l2#nu + $ARGV[1]";
my $macro = `mktemp`;
$macro =~ s/\n//;
open(OUT, ">$macro") || die "Failed to write to a ROOT macro: $macro\n$!\n";
print OUT <<EOF;
{
  gSystem->Load("lands.so");
  gSystem->CompileMacro("PlotExpectedLimits.C","k");
  PlotExpectedLimits("$dir","$name","$title","$type");
}
EOF
close OUT;
die "Failed to post-process\n" if (system("root -q -b $macro")!=0);

system("epstopdf limits.eps");
system("mv limits.pdf $dir/$name-$type.pdf");
system("mv limits.gif $dir/$name-$type.gif");
exit
