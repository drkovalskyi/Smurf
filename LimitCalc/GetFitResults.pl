#!/bin/env perl
# 
# Script to run multiple LandS calculators. 
#
#
use strict;
use warnings;
use threads;
use threads::shared;

my $lands = "../../LandS/test/lands.exe";

my $nMaxThreads = 8;

#########################################################################
my $usage = "Usage:\n \t$0 <dir>/<card name> \n\n";
die $usage unless @ARGV == 1;
die "LandS not found. Please fix it.\n" unless -e $lands;

my $type = "ProfiledLikelihood";
$type = $ARGV[1] if (@ARGV == 2);

my @threads = ();

my $command = "--bDumpFitResults -M ScanningMuFit --scanRs 1 -vR 0  --minuitSTRATEGY 1 --maximumFunctionCallsInAFit 50000";

my %cards = ();
my $inputCard = $ARGV[0];
my ($dir,$name) = ($inputCard =~ /^(.+?)\/([^\/]+)$/);
die $usage if (!defined $dir || !defined $name); 
$name =~ s/.*?([^\/]+)$/$1/;
$name =~ s/\.[^\.]+$//;
foreach my $line (`cat $inputCard`){
    $line =~ s/\n//;
    next if ($line =~ /^\s*\#/);
    next if ($line !~ /\S/);
    my @elements = split(/\s+/,$line);
    die "Wrong format of input file. Each line show be: <mass> <card1> [<card2>...]\n" 
	unless @elements>1;
    my $mass = shift @elements;
    $cards{$mass} = \@elements;
}
print "Number of mass points to process: ", scalar keys %cards,"\n";

my $path = `pwd`;
$path =~ s/\n//;

print "Processing cards:\n";
system("mkdir $dir/output") unless -d "$dir/output" || -l "$dir/output";
foreach my $mass(sort {$a<=>$b} keys %cards){
    while (threads->list(threads::running)>=$nMaxThreads){ sleep 1; }
    my $cards = "";
    foreach my $card(@{$cards{$mass}}){
	my $fullcard = "$path/$dir/$card";
	$cards .= " $fullcard";
	my $n = "$name-fit-$card";
	$n =~ s/\//-/g;
	push @threads, 
	threads->new(\&run_job, "$lands $command --name $dir/output/$n -d $fullcard>$dir/output/$n.log 2>&1");
    }
    my $n = "$name-fit-$mass-all";
    push @threads, threads->new(\&run_job, "$lands $command --name $dir/output/$n -d$cards>$dir/output/$n.log 2>&1");
}

sub run_job{
    my $command = shift;
    print ".";
    system($command);
}

while (threads->list(threads::running)){
    sleep 1;
}
foreach my $thread(@threads){
    $thread->join();
}
print "\nDone\n";

my $nMassPoints = scalar keys %cards;
my %results = ();
my $macro=`mktemp`;

foreach my $mass(sort {$a<=>$b} keys %cards){
    foreach my $card(@{$cards{$mass}}){
	my $fullcard = "$path/$dir/$card";
	my $n = "$name-fit-$card";
	my $n2 = "$name-fit-all-$card";
	my $nall = "$name-fit-$mass-all";
	$n =~ s/\//-/g;
	$n2 =~ s/\//-/g;
	if (my ($rootfile) = (`grep histo_Data $dir/$card` =~ /(\S+?\.root)\s+histo_Data/)){
	    if ($card =~ /hww(..)_(\d)j_shape/){
		open(OUT,">$macro")||die ("Cannot write to file $macro\n$!\n");
		print OUT '{ gROOT->LoadMacro("visualizeFitResults.C+");';
		print OUT "\nmakePlots(\"j$2$1\",\"$dir/output/$n\",\n\"$rootfile\");\n}\n"; 
		close OUT;
		system("root -b -q $macro");
		open(OUT,">$macro")||die ("Cannot write to file $macro\n$!\n");
		print OUT '{ gROOT->LoadMacro("visualizeFitResults.C+");';
		print OUT "\nmakePlots(\"j$2$1\",\"$dir/output/$nall\",\n\"$rootfile\",\n\"$dir/output/$n2\");\n}\n"; 
		close OUT;
		system("root -b -q $macro");
	    } else {
		die "failed to undersand what channel was used\n";
	    }
	} else {
	    die "failed to fine data shape input root file\n";
	}
    }
}


exit
