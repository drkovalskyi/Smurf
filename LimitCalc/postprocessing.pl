#!/bin/env perl
use strict;
use warnings;
die "Usage:\n\t$0 <file with a list of cards> <njets>\n" unless @ARGV == 2;
my $name = $ARGV[0];
$name =~ s/.*?([^\/]+)$/$1/;
$name =~ s/\.[^\.]+$//;
my $title = "H #rightarrow WW #rightarrow 2l2#nu + $ARGV[1]";
system("root -q -b 'makePlots.cxx(\"$name\",\"$title\")'");
system("epstopdf limits.eps");
system("mv limits.pdf $name.pdf");
exit
