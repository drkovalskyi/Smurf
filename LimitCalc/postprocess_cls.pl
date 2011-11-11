#!/bin/env perl
use strict;
use warnings;
my $usage = "Usage:\n\t$0 <path>/<file with a list of cards> <njets>\n";
die $usage unless @ARGV == 2;
my ($dir,$name) = ($ARGV[0] =~ /^(.+?)\/([^\/]+)$/);
die $usage if (!defined $dir || !defined $name); 
$name =~ s/.*?([^\/]+)$/$1/;
$name =~ s/\.[^\.]+$//;
my $pwd = `pwd`;
$pwd =~ s/\n//;
my %seeds = ();
my $title = "H #rightarrow WW #rightarrow 2l2#nu + $ARGV[1]";
    foreach my $file (`find $dir/output|grep $name|grep m2lnQ.root`){
	$file =~ s/\n//;
	next if ($file !~ /\.root$/);
	if (my ($mass,$seed) = ($file =~ /(\d+)-(\d+)_m2lnQ.root/)){
	    system("cd ../../LandS/test/ ; root -b -q \'fitRvsCLs.C+(\"$pwd/$file\",\"plot\")\'>$pwd/$file.results")
		unless -e "$file.results";
	    next unless -e "$file.results";
	    my $result = `cat $file.results | grep BANDS`;
	    $result =~ s/.*?://; 
	    $result =~ s/\+.*?\s/ /g;
	    $result =~ s/^\s+//;
	    $result =~ s/\s+$//;
	    my @results = split(/\s+/,$result);
	    die "Unepxected number of columns in the result: ", scalar @results, "\n$result\n" unless (@results == 6);
	    $seeds{$seed}->{$mass} = $result;
	}
	
    }
foreach my $seed (sort keys %seeds){
    open (OUT, ">$dir/output/$name-report.$seed");
    foreach my $mass (sort {$a<=>$b} keys %{$seeds{$seed}}){
	my $result = $seeds{$seed}->{$mass};
	# $result =~ s/\s+/\t /g;
	print OUT "$mass $result\n";
    }
    close OUT;
}
# system("root -q -b 'makePlots.cxx(\"$dir\",\"$name\",\"$title\")'");
# system("epstopdf limits.eps");
# system("mv limits.pdf $dir/$name.pdf");
# system("mv limits.gif $dir/$name.gif");
exit
