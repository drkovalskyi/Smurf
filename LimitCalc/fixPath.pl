#!/bin/env perl
use strict;
use warnings;
die "Usage:\n\t $0 <directory>\n\nThis script will add absolute path to all root file names mentioned in cards\n"
    unless @ARGV==1;

foreach my $file (`find $ARGV[0] -type f`){
    next unless $file =~ /\.txt$/;
    $file =~ s/\n//;
    my $path = $file;
    if ( $path =~ /(.*?)\/[^\/]+$/ ){
	$path = $1;
    } else {
	$path = "./";
    }
    my $pwd = `cd $path; pwd`;
    $pwd =~ s/\n//;
    my $input = `cat $file`;
    my $output = $input;
    $output =~ s/(\S*?)([^\/\s]+\.root)/$pwd\/$2/gm;
    if ($input ne $output){
	print "patched $file\n";
	open(OUT,">$file")||die "Cannot write to file: $file\n$!\n";
	print OUT $output;
	close OUT;
    }
}

# perl -pi -e 's/\/smurf\/ceballos\/tmva\/inputLimits_V6/inputLimits/' inputs_0j_200pb_shape_allann/*.txt

exit;
