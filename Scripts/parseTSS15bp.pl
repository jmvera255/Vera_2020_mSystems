#! /usr/bin/perl
use strict;
#### this script will take a TSS bed file and remove those TSS that are <= 15bp 
#### downstream of another TSS

if(@ARGV != 1){
	die "Must provide a bed file!\n";
}

my %bed;

open(BED, "< $ARGV[0]") || die "Cannot open bed file!\n";
while(my $line = <BED>){
	chomp($line);
	if($line !~ /^#/){
		my @tabs = split("\t", $line);
		push(@{$bed{$tabs[0]}{$tabs[5]}{$tabs[1]}}, $tabs[3]);
	}
}

foreach my $chr (sort {$a cmp $b} keys %bed){
	foreach my $strand (keys %{$bed{$chr}}){
		my @positions = (sort {$a <=>$b} keys %{$bed{$chr}{$strand}});
		if($strand =~ /-/){
			for(my $i == 0; $i < scalar(@positions) - 1; $i++){
				my $j = $i + 1;
				if($positions[$j] > $positions[$i] + 15){
					push(my @line, $chr, $positions[$i], $positions[$i] + 1, $bed{$chr}{$strand}{$positions[$i]}[0], "0", $strand);
					print join("\t", @line) . "\n";
				}
			}
			my $j = scalar(@positions) - 1;
			#my $i = $j - 1;
			#if($positions[$j] > $positions[$i] + 15){
				push(my @line, $chr, $positions[$j], $positions[$j] + 1, $bed{$chr}{$strand}{$positions[$j]}[0],"0", $strand);
				print join("\t", @line) . "\n";
			#}
		}
		else{
			for(my $i = scalar(@positions) - 1; $i > 0; $i--){
				my $j = $i - 1;
				if($positions[$j] < $positions[$i] - 15){
					push(my @line, $chr, $positions[$i], $positions[$i] + 1, $bed{$chr}{$strand}{$positions[$i]}[0], "0", $strand);
					print join("\t", @line) . "\n";
				}
			}
			#if($positions[0] < $positions[1] - 15){
				push(my @line, $chr, $positions[0], $positions[0] + 1, $bed{$chr}{$strand}{$positions[0]}[0],"0", $strand);
				print join("\t", @line) . "\n";
			#}
		}
	}
}
