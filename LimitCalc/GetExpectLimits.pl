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
# my $expectLimitCommand = "-M Hybrid --freq --ExpectationHints Asymptotic --scanRs 1 --freq --nToysForCLsb 1500 --nToysForCLb 500"; 
# my $expectLimitCommand = "-M Hybrid --freq --ExpectationHints Asymptotic --scanRs 1 --freq --nToysForCLsb 1000 --nToysForCLb 500"; 
# my $expectLimitCommand = "-M Hybrid --freq --ExpectationHints Asymptotic --scanRs 1 --freq --nToysForCLsb 4000 --nToysForCLb 2000"; 
# my $expectLimitCommand = "-M Hybrid --freq --ExpectationHints Bayesian --scanRs 1 --freq --nToysForCLsb 10000 --nToysForCLb 5000"; 
# my $expectLimitCommand = "-tB 1000 -tPB 30 -t 1000 -M Bayesian --doExpectation 1";
my $expectLimitCommand = "-tB 1000 -tPB 30 -t 100 -M Bayesian --doExpectation 1";
my $observedLimitCommand = "-tB 10000 -M Bayesian";
# my $observedLimitCommand = "-tB 100000 -M Hybrid";
my $useLSF = 1; 
my $doObserved = 0;
my $doExpected = 1;
my $nJobsPerMassPoint = 10;

#########################################################################
my $usage = "Usage:\n \t$0 <dir>/<card name>\n";
die $usage unless @ARGV == 1;
die "LandS not found. Please fix it.\n" unless -e $lands;

sub seed{
    return int(rand()*10000);
}

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

system("mkdir $dir/output") unless -d "$dir/output";
foreach my $job(1..$nJobsPerMassPoint){
    my $seed = seed();
    foreach my $mass(sort {$a<=>$b} keys %cards){
	my $cards = "";
	foreach my $card(@{$cards{$mass}}){
	    $cards .= " $dir/$card";
	}
	my $n = "$name-$mass-$seed";
	if ( $doExpected ) {
	    if ($useLSF){
		system("bsub -q cmscaf1nd -o $dir/output/$n.log \"$wrapper $expectLimitCommand --seed $seed --name $dir/output/$n -d$cards\"");
	    } else {
		system("nice $lands $expectLimitCommand --seed $seed --name $dir/output/$n -d$cards>$dir/output/$n.log 2>&1 &");
	    }
	}
	if ( $doObserved  && $job==1) {
	    print "nice $lands $observedLimitCommand -d$cards>$dir/output/$name-$mass.observed 2>&1 &\n";
	    system("nice $lands $observedLimitCommand -d$cards>$dir/output/$name-$mass.observed 2>&1 &");
	}
    }
}

exit
