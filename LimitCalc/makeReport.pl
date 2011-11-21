#!/bin/env perl
use strict;
use warnings;

my $usage = "Usage:\n \t$0 <dir>/<card name>\n";
die $usage unless @ARGV == 1;
my %cards = ();
my $inputCard = $ARGV[0];
my ($dir,$name) = ($inputCard =~ /^(.+?)\/([^\/]+)$/);
die $usage if (!defined $dir || !defined $name); 
$name =~ s/.*?([^\/]+)$/$1/;
$name =~ s/\.[^\.]+$//;
my $nCardsPerMassPoint = 0;
foreach my $line (`cat $inputCard`){
    $line =~ s/\n//;
    next if ($line =~ /^\s*\#/);
    next if ($line !~ /\S/);
    my @elements = split(/\s+/,$line);
    die "Wrong format of input file. Each line show be: <mass> <card1> [<card2>...]\n" 
	unless @elements>1;
    my $mass = shift @elements;
    $nCardsPerMassPoint = @elements if ($nCardsPerMassPoint == 0);
    die "Inconsistent number of cards per mass point. It should be the same number to parse it right\n"
	unless $nCardsPerMassPoint == @elements;
    $cards{$mass} = \@elements;
}
die "0 cards is found\n" unless $nCardsPerMassPoint>0;

sub FilterColumns{
    my @input = @_;
    if ($input[1] =~ /lnN/){
#	@input = @input[0,5..$#input];
	@input = @input[0,4..$#input];
    } else {
#	@input = @input[0,4..$#input];
	@input = @input[0,3..$#input];
    }
    return @input;
}


my $tableBegin = << 'EOF';
\begin{table}
{%\footnotesize
 \tiny
 \begin{center}
 \begin{tabular}{SEPARATORS}
 \hline
 TITLE
 \hline
EOF

my $tableEnd = << 'EOF';
\hline
\end{tabular}
\end{center}
}
\caption{TEXT}
\end{table}
EOF

open(OUT,">$dir/$name-report.tex");
print OUT << 'EOF';
\documentclass{report}
\usepackage{fullpage}
\begin{document}
EOF

# print "Number of mass points to process: ", scalar keys %cards,"\n";
foreach my $card(0..$nCardsPerMassPoint-1){
    my $title;
    my @names = ("") x $nCardsPerMassPoint;
    foreach my $mass( sort {$a<=>$b} keys %cards ){
	my $cardname = $cards{$mass}->[$card];
	my $file = "$dir/$cardname";
	$cardname =~ s/(\d+\/)//;
	$names[$card] = $cardname;
	# print "$file\n";
	my @columns;
	my @errors2;
	my $observation;
	foreach my $line(`cat $file`){
	    if ( !defined $title && ($line =~ /process/) ){
		$title = $line;
		$title =~ s/\n//;
		my @titleColumns = split(/\s+/,$title); 
		@titleColumns = FilterColumns(@titleColumns);
		push @titleColumns, "\$\\sum\$Bkg";
		push @titleColumns, "Data";
		my $separators = "l | c c | ".("c "x(@titleColumns - 5))." | c c";
		$title = join(" & ",@titleColumns);
		$title .= " \\\\";
		my $tab = $tableBegin;
		$tab =~ s/SEPARATORS/$separators/m;
		$tab =~ s/TITLE/$title/m;
		print OUT $tab;
            }
	    $line =~ s/\n//;
	    if ($line =~ /rate/){
		$line =~ s/rate/$mass/;
		@columns = FilterColumns(split(/\s+/,$line));
# 		foreach my $i(1..$#columns){
# 		    $columns[$i] = sprintf("%.1f", $columns[$i]);
# 		}
	    }
	    if ($line =~ /lnN/){
		my @errors = FilterColumns(split(/\s+/,$line));
		die "Errors and yeild don't match\n" unless (scalar @errors == scalar @columns);
		if (! @errors2){
		    @errors2 = @errors;
		    foreach my $i(1..$#errors2){
			my $err = $errors2[$i];
			$err = 1 if ($err eq "-");
			$errors2[$i] = ($err-1)*($err-1);
		    }
		} else {
		    foreach my $i(1..$#errors2){
			my $err = $errors[$i];
			$err = 1 if ($err eq "-");
			$errors2[$i] += ($err-1)*($err-1);
		    }
		}
	    }   
	    if ($line =~ /Observation\s+(\d+)/){
		$observation = $1;
	    }
	}
	my $sum = 0;
	my $sumErr2 = 0;
	foreach my $i(1..$#columns){
	    if ($i>2){ # skip first two columns, which are qqH and ggH
		$sum += $columns[$i];
		$sumErr2 += $errors2[$i]*$columns[$i]*$columns[$i];
	    }
	    $columns[$i] = sprintf("\$%0.1f\\pm%0.1f\$", $columns[$i], $columns[$i]*sqrt($errors2[$i]));
	}
	push @columns, sprintf("\$%0.1f\\pm%0.1f\$", $sum, sqrt($sumErr2));
	push @columns, "$observation";
	print OUT join(" & ",@columns)." \\\\\n";
    }
    my $tab = $tableEnd;
    my $label = $names[$card];
    $label =~ s/_/\\_/g;
    $tab =~ s/TEXT/Summary of card $label/m;
    print OUT $tab;
}
print OUT '\end{document}'."\n";
close OUT;
exit
