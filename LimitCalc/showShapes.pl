#!/bin/env perl
use strict;
use warnings;
die "Usage:\n\t$0 <card for shape analysis>\n\n" unless @ARGV==1;

my %colors = ("Data"   => "kBlack", 
	      "ZH"     => "kRed+1", 
	      "WH"     => "kRed+1",
	      "qqH"    => "kRed+1",
	      "ggH"    => "kRed+1",
	      "qqWW"   => "kAzure-9",
	      "ggWW"   => "kAzure-9",
	      "VV"     => "kAzure-2",
	      "Top"    => "kYellow", 
	      "Zjets"  => "kGreen+2",
	      "Wjets"  => "kGray+1", 
	      "Wgamma" => "kAzure-2");

my $path;
my %processes = ();
foreach my $line (`cat $ARGV[0]`){
    $line =~ s/\n//;
    $line =~ s/^\s*//;
    $line =~ s/\s*$//;
    if ( $line =~ /^process\s+[^\d\-]/){
	my @p = split(/\s+/,$line);
	shift @p;
	foreach my $process (@p){
	    $processes{$process} = [];
	}
    }
}

foreach my $line (`cat $ARGV[0]`){
    $line =~ s/\n//;
    $line =~ s/^\s*//;
    $line =~ s/\s*$//;
    if (my ($process,$channel,$file,$histogram) = ($line =~ /^\s*shapes\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)/) ){
	die "No support for channels. Contact the developer\n" unless ($channel eq "*");
	if ($process eq "data_obs"){
	    $processes{"Data"} = [$file, $histogram];
	} 
	if ($process eq "*"){
	    foreach my $process(sort keys %processes){
		next if ($process eq "Data");
		my $h = $histogram;
		$h =~ s/\$PROCESS/$process/;
		$processes{$process} = [$file, $h];
	    }
	}
    }
}

my @histograms = ();
foreach my $process(sort keys %processes){
    my ($file, $hist) = @{$processes{$process}}; 
    # print "$process \t$file \t$hist\n";
}

# make macro
my $macroFile = $ARGV[0].".plot.cxx";
my $plotFile = $ARGV[0].".png";
open (OUT,">$macroFile")||die"Cannot write to file $macroFile\n$!\n";
print OUT "{\n";
print OUT "gROOT->SetStyle(\"Plain\");\n";
print OUT "gStyle->SetOptTitle(0);\n";
print OUT "THStack* stack = new THStack(\"hs\",\"Predictions\");\n";
my ($fd, $hd) = @{$processes{"Data"}}; 
print OUT "TFile* f = TFile::Open(\"$fd\");\n";
print OUT "TH1* $hd = (TH1*)f->Get(\"$hd\");\n";
print OUT "TH1* hdata = $hd;\n";
print OUT "$hd->SetMarkerStyle(20);\n";
print OUT "TH1* hobs = (TH1*)$hd->Clone(\"hobs\");\n";
my $i=0;
foreach my $process(sort keys %processes){
    next if ($process eq "Data" || $process eq "ggH" || $process eq "qqH" || 
	     $process eq "WH" || $process eq "ZH");
    my ($file, $hist) = @{$processes{$process}}; 
    print OUT "TFile* f$i = TFile::Open(\"$file\");\n";
    print OUT "TH1* $hist = (TH1*)f$i->Get(\"$hist\");\n";
    print OUT "hobs->Add($hist,-1);\n";
    my $color = "kBlack";
    $color = $colors{$process} if (defined $colors{$process});
    print OUT "$hist->SetFillColor($color);\n";
    print OUT "stack->Add($hist);\n";
    $i++;
}   
print OUT "TH1* hsig(0);\n";
foreach my $process(sort keys %processes){
    next if ($process ne "ggH" && $process ne "qqH" && 
	     $process ne "WH" && $process ne "ZH");
    my ($file, $hist) = @{$processes{$process}}; 
    print OUT "TFile* f$i = TFile::Open(\"$file\");\n";
    print OUT "TH1* $hist = (TH1*)f$i->Get(\"$hist\");\n";
    print OUT "if(!hsig) {hsig = (TH1*)$hist->Clone(\"hsig\");} else {hsig->Add($hist);}\n";
    $i++;
}
print OUT "hsig->SetLineColor(kRed+1);hsig->SetFillColor(0); hsig->SetLineWidth(2);\n";
print OUT << "EOF";
  double maximum = stack->GetMaximum();
  if (maximum<histo_Data->GetMaximum()) maximum = histo_Data->GetMaximum();
  maximum *= 1.2;
  stack->SetMaximum(maximum);
  hdata->SetMaximum(maximum);
  new TCanvas("cc","cc",600,600);
  cc->cd();
  TPad *pad1 = new TPad("p_main", "p_main", 0.0, 0.3, 1.0, 1.0);
  pad1->SetLogy(0);
  pad1->SetBottomMargin(0.05);
  pad1->SetRightMargin(0.07);
  pad1->Draw();
  pad1->cd();
  stack->Draw("hist");
  hdata->Draw("e1 same");
  cc->cd();
  TPad *pad2 = new TPad("p_pull", "p_pull", 0.0, 0.0, 1.0, 0.3);
  pad2->Draw();
  pad2->SetLogy(0);
  pad2->SetGridy(1);
  pad2->cd();
  pad2->SetTopMargin(0.01);
  pad2->SetRightMargin(0.07);
  pad2->SetBottomMargin(0.4);
  hobs->SetStats(0);
  hobs->Draw("e1");
  hsig->SetStats(0);
  hsig->Draw("same hist");
EOF
print OUT "}";
close OUT;


exit
