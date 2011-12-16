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

my $useLSF = 1; 
my $doObserved = 0;
my $doExpected = 1;
my $nJobsPerMassPoint = 10;
# my $outputHostName = "cmswn1927.cern.ch";
my $outputHostName = "lxplus.cern.ch";
my $queue = "1nd"; # 1 day
# my $queue = "cmscaf1nd"; # 1 day
# my $queue = "cmscaf1nh"; # 1 hour
#########################################################################
my $usage = "Usage:\n \t$0 <dir>/<card name> <type>\n\nValid types are: Bayesian, CLs, CLs-asymptotic\n\n";
die $usage unless @ARGV == 2;
die "LandS not found. Please fix it.\n" unless -e $lands;

my $type = $ARGV[1];
my $jobname = "";
my $expectLimitCommand;
my $observedLimitCommand;

##
## Bayesian
##
if ( $type eq "Bayesian" ){
    $jobname = "bayesian";
    $expectLimitCommand = "-tB 1000 -tPB 30 -t 100 -M Bayesian --doExpectation 1";
    $observedLimitCommand = "-tB 10000 -M Bayesian";
}

##
## CLs - asymptotic
##
if ( $type eq "CLs-asymptotic" ){
    $jobname = "cls_asym";
    $expectLimitCommand = "-M Asymptotic -D asimov_b --minuitSTRATEGY 1 --maximumFunctionCallsInAFit 50000 -rMin 0.1 -rMax 100";
    $observedLimitCommand = "-M Asymptotic --minuitSTRATEGY 1 --maximumFunctionCallsInAFit 50000 -rMin 0.1 -rMax 100";
    $nJobsPerMassPoint = 1;
    $useLSF = 0;
}

##
## CLs
##
if ( $type eq "CLs" ){
    $nJobsPerMassPoint = 100;
    $jobname = "cls";
    # GOOD $expectLimitCommand   = "-M Hybrid --freq --ExpectationHints Asymptotic --scanRs 1 --nToysForCLsb 50 --nToysForCLb 25 --minuitSTRATEGY 1 --maximumFunctionCallsInAFit 50000 -rMin 0.01 -rMax 20"; 
    $expectLimitCommand   = "-M Hybrid --freq -vR \"[0.1,2,0.1] [2.2,5.0,0.2] [6.0,15.0,1]\" --scanRs 1 --nToysForCLsb 50 --nToysForCLb 30 --minuitSTRATEGY 1 --maximumFunctionCallsInAFit 50000"; 
    $observedLimitCommand = "-M Hybrid --freq --nToysForCLsb 10000 --nToysForCLb 5000 --minuitSTRATEGY 1 --maximumFunctionCallsInAFit 50000 -rMin 0.01 -rMax 20";
    $doObserved = 0;
}

die "Incorrect type $type\n" 
    unless (defined $expectLimitCommand && defined $observedLimitCommand);

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
my $wrapper = "process_job.sh";
my $path = `pwd`;
$path =~ s/\n//;
# if ( ! -e $wrapper ){
open(OUT,">$wrapper") || die "Cannot create file $wrapper\n$!\n";
print OUT << "EOF";
#!/bin/sh
# first argument is the output directory on the host
# second is log file
# the rest is for LandS
workdir=$path
remotehost="$outputHostName"
lands=$path/$lands
EOF
print OUT << 'EOF';
TMPDIR=`pwd`
SCRATCH=`mktemp -d`
remotedir=$1
logfile=$2
source /afs/cern.ch/cms/sw/cmsset_default.sh
cd $workdir
eval `scram runtime -sh`
cd $SCRATCH
$lands ${@:3} > $logfile 2>&1
gzip $logfile
scp * ${remotehost}:${remotedir} && rm *
EOF
close OUT;
system "chmod u+x $wrapper";

system("mkdir $dir/output") unless -d "$dir/output" || -l "$dir/output";
system("mkdir $dir/log") unless -e "$dir/log";
foreach my $job(1..$nJobsPerMassPoint){
    my $seed = seed();
    foreach my $mass(sort {$a<=>$b} keys %cards){
	my $cards = "";
	foreach my $card(@{$cards{$mass}}){
	    $cards .= " $path/$dir/$card";
	}
	my $n = "$name-$mass-$jobname";
	$n .= "-$seed" if ($nJobsPerMassPoint>1);
	if ( $doExpected ) {
	    if ($useLSF){
		system("bsub -q $queue -o $dir/log/$n.lsf \"$wrapper $path/$dir/output $n.log $expectLimitCommand --seed $seed --name $n -d$cards\"");
	    } else {
		system("nice $lands $expectLimitCommand --seed $seed --name $dir/output/$n -d$cards>$dir/output/$n.log 2>&1 &");
	    }
	}
	if ( $doObserved  && $job==1) {
	    # print "nice $lands $observedLimitCommand -d$cards>$dir/output/$name-$mass-$jobname.observed 2>&1 &\n";
	    $n = "$name-$mass-$jobname";
	    if ($useLSF){
		system("bsub -q $queue -o $dir/log/$n.lsf \"$wrapper $path/$dir/output $n.observed $observedLimitCommand -d$cards\"");
	    } else {
		system("nice $lands $observedLimitCommand -d$cards>$dir/output/$name-$mass-$jobname.observed 2>&1 &");
	    }
	}
    }
}

exit
