#!/bin/env perl
use strict;
use warnings;
my $usage = "Usage:\n\t$0 <path>/<file with a list of cards> <n-jets title> <type>\n\nValid types are: Bayesian, CLs, CLs-asymptotic\n\n";
die $usage unless @ARGV == 3;
my $inputCard = $ARGV[0];
my ($dir,$name) = ($inputCard =~ /^(.+?)\/([^\/]+)$/);
die $usage if (!defined $dir || !defined $name); 
$name =~ s/.*?([^\/]+)$/$1/;
$name =~ s/\.[^\.]+$//;
my $type = $ARGV[2];
my $title = $ARGV[1];
my $macro = `mktemp`;
$macro =~ s/\n//;

if ( $type eq "CLs-asymptotic" || $type eq "Bayesian" ){
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
}

if ( $type eq "CLs" ){
    my $pwd = `pwd`;
    $pwd =~ s/\n//;
    my %seeds = ();
    # merge multiple seeds into one m2lnQ file
    my %cards = ();
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
    my %results = ();
    foreach my $mass(sort {$a<=>$b} keys %cards){
	my $n = "$name-$mass-cls";
	my $filename = "$dir/output/$n.root";
	system("hadd $filename $dir/output/$n-*_m2lnQ.root")
	    unless -e $filename;
	# system("rm -v $dir/output/$n-*root");
	if (! -e "$filename.results"){
	    system("rm ../../LandS/test/plot_*gif");
	    system("cd ../../LandS/test/ ; root -b -q \'fitRvsCLs.C+(\"$pwd/$filename\",\"plot\")\'>$pwd/$filename.results");
	    system("mkdir -p $dir/log/$mass") unless -e "$dir/log/$mass";
	    system("cp -v ../../LandS/test/plot_*gif $dir/log/$mass");
	}
	next unless -e "$filename.results";
	my $result = `cat $filename.results | grep BANDS`;
	$result =~ s/.*?://; 
	$result =~ s/\+.*?\s/ /g;
	$result =~ s/^\s+//;
	$result =~ s/\s+$//;
	my @results = split(/\s+/,$result);
	die "Unepxected number of columns in the result: ", scalar @results, "\n$result\n" unless (@results == 6);
	$results{$mass} = $result;
    }

    open (OUT, ">$dir/output/$name-cls-report");
    foreach my $mass (sort {$a<=>$b} keys %results){
	my $result = $results{$mass};
	# $result =~ s/\s+/\t /g;
	print OUT "$mass $result\n";
    }
    close OUT;

    open(OUT, ">$macro") || die "Failed to write to a ROOT macro: $macro\n$!\n";
    print OUT <<EOF;
    {
	gSystem->Load("lands.so");
	gSystem->CompileMacro("PlotExpectedLimits.C","k");
	PlotExpectedLimits("$dir/output/$name-cls-report","$title");
    }
EOF
    close OUT;
    die "Failed to post-process\n" if (system("root -q -b $macro")!=0);
}

system("epstopdf limits.eps");
system("mv limits.pdf $dir/$name-$type.pdf");
system("mv limits.gif $dir/$name-$type.gif");
exit
