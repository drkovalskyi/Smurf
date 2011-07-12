#!/bin/env perl
# 
# Script to run multiple LandS calculators. It can be used on a local
# machine (24 uaf monsters are the best) or in a batch if you have CAF
# account at CERN.
#
#
use strict;
use warnings;


my $lands = "../../LandS/test/lands.exe";
my $expectLimitCommand = "-tB 1000 -tPB 30 -t 1000 -M Bayesian --doExpectation 1";
# my $expectLimitCommand = "-tB 1000 -tPB 30 -t 100 -M Bayesian --doExpectation 1";
my $observedLimitCommand = "-tB 100000 -M Bayesian";
my $useLSF = 1; 
my $doObserved = 1;
my $doExpected = 1;
my $nJobsPerMassPoint = 10;

#########################################################################
die "Usage:\n\t$0 <file with a list of cards>\n" unless @ARGV == 1;
die "LandS not found. Please fix it.\n" unless -e $lands;

sub seed{
    return int(rand()*10000);
}

my %cards = ();
my $name = $ARGV[0];
$name =~ s/.*?([^\/]+)$/$1/;
$name =~ s/\.[^\.]+$//;
foreach my $line (`cat $ARGV[0]`){
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

# make a wrapper script
my $wrapper = "lands.sh";
my $path = `pwd`;
$path =~ s/\n//;
if ( ! -e $wrapper ){
    open(OUT,">$wrapper") || die "Cannot create file $wrapper\n$!\n";
    print OUT << "EOF";
#!/bin/sh
source /afs/cern.ch/cms/sw/cmsset_default.sh
cd $path
eval `scram runtime -sh`
$lands \$*
EOF
    close OUT;
    system("chmod u+x $wrapper");
}

system("mkdir output") unless -d "output";
foreach my $mass(sort {$a<=>$b} keys %cards){
    foreach (1..$nJobsPerMassPoint){
	my $seed = seed();
	my $n = "$name-$mass-$seed";
	if ( $doExpected ) {
	    if ($useLSF){
		system("bsub -q cmscaf1nd -o output/$n.log \"$wrapper $expectLimitCommand --seed $seed --name output/$n -d ".join(" ",@{$cards{$mass}})."\"");
	    } else {
		system("nice $lands $expectLimitCommand --seed $seed --name output/$n -d ".join(" ",@{$cards{$mass}}).">output/$n.log 2>&1 &");
	    }
	}
    }
    if ( $doObserved ) {
	print "nice $lands $observedLimitCommand -d ".join(" ",@{$cards{$mass}}).">output/$name-$mass.observed 2>&1 &\n";
	system("nice $lands $observedLimitCommand -d ".join(" ",@{$cards{$mass}}).">output/$name-$mass.observed 2>&1 &");
    }
}

exit