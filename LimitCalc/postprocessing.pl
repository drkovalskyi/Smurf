#!/bin/env perl
use strict;
use warnings;
my $usage = "Usage:\n\t$0 <path>/<file with a list of cards> <njets>\n";
die $usage unless @ARGV == 2;
my ($dir,$name) = ($ARGV[0] =~ /^(.+?)\/([^\/]+)$/);
die $usage if (!defined $dir || !defined $name); 
$name =~ s/.*?([^\/]+)$/$1/;
$name =~ s/\.[^\.]+$//;
my $title = "H #rightarrow WW #rightarrow 2l2#nu + $ARGV[1]";
system("root -q -b 'makePlots.cxx(\"$dir\",\"$name\",\"$title\")'");
system("epstopdf limits.eps");
system("mv limits.pdf $dir/$name.pdf");
system("mv limits.gif $dir/$name.gif");
exit
