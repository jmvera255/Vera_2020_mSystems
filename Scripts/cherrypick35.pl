#! /usr/bin/perl

### this script differs from cherrypick35.v4 in that there is an additional round of cherrypicking of
### -35 preliminary sequences after malign

use strict;
use DateTime;
my $dt = DateTime->now;
my $start = $dt->hms;
print STDERR "Start time = $start\n";

use List::Util qw( min max);
use List::MoreUtils qw(first_index);
use Getopt::Std;
my (%opts);
getopts("i:l:s:T:p:", \%opts);


if(@ARGV != 1){
	die "No file provided!

-i <int> = window start, 0-based count from left
-l <int> = window length; window of sequence to calculate information content
-s <int> = slide value
-T <int> = threshold for stopping malign
-p <int> = number of passes to make\n\n";
}

my(%passesHash, %testPasses, %slides);
my @alpha = ("A", "G", "T", "C");
	
my $fasta_ref = readFasta($ARGV[0]);
#my $Riw_ref = buildRiw($fasta_ref);
my($cherryHash_ref, $n, $cherryInfo) = cherrypick($fasta_ref);
my $dflogf_ref = malignPreCalc($n); ### is based on the # of sequences in alignment, only needs to be made once

##malign(cherryPickedSequences, precalc dflogf, slide, windowStart, windowLength
my $malignHash_ref = malign($cherryHash_ref, $dflogf_ref, $cherryInfo);

### return the best shuffle
my %bestAlign;
my @bestPasses = sort({$b <=> $a} keys %passesHash);
my $bestPass = $passesHash{$bestPasses[0]};
print STDERR "pass $bestPass had the best alignment w/an average information content of $bestPasses[0]\n";
#my %hash = %$malignHash_ref;
my %hash = %$cherryHash_ref;
foreach my $TSS (keys %hash){
	my @tabs = split("", $hash{$TSS});
	my $s = $slides{$TSS}[$bestPass];
	if($s ==0){
	}
	elsif($s < 0){
		###shift sequence in @tabs to the left with shift and push
		for(my $z = $s; $z < 0; $z++){
			my $temp = shift(@tabs);
			push(@tabs, $temp);
		}
	}
	else{
		### shift sequence in @tabs to the right with pop and unshift
		for(my $z = $s; $z > 0; $z--){
			my $temp = pop(@tabs);
			unshift(@tabs, $temp);
		}
	}
	#print join("", @tabs) . "\n";
	$bestAlign{$TSS} = join("", @tabs);
}
###cherrypick best shuffle
my $fasta_ref = \%bestAlign;
my($cherryHash_ref, $n, $cherryInfo) = cherrypick($fasta_ref);
my %bestAlign = %$cherryHash_ref;
foreach my $TSS (keys %bestAlign){
	print "$bestAlign{$TSS}\n";
}

my $dt = DateTime->now;
my $start = $dt->hms;
print STDERR "End time = $start\n";

sub malign{

	### define variables passed to sub
	my $hash_ref = $_[0];
	my %cherryFasta = %$hash_ref;
	my $hash_ref = $_[1];
	my %dflogf = %$hash_ref;


	my $testHash_ref = \%cherryFasta;	
	#print STDERR "Begining malign, test information content = $newInfo\n";
	print STDERR "Begining malign. Starting average information content = $_[2]\n";
	my $newInfo = $$_[2];

	###prep variables for first shuffle of first pass
	my %tempCherry = %cherryFasta; ###tempCherry is equivalent to the original expanded cherryFasta
	$testPasses{0} = \%cherryFasta;
	for (my $pass = 1; $pass <= $opts{p}; $pass++){	
		foreach my $seq (keys %cherryFasta){
			#my @tabs = split("", $cherryFasta{$seq});
			my @tabs = split("", $tempCherry{$seq});
			$slides{$seq}[0] = 0;
			my $previousPass = $pass - 1;
			#print STDERR "$seq $tempCherry{$seq} for pass# $pass, previous pass = $previousPass ";
			#print STDERR "which is a slide of $slides{$seq}[$previousPass]\n";
			#my %tempCherry = %cherryFasta;
			delete($tempCherry{$seq});  ###remove current sequence for expanded cherry
			my $tempCherryHash_ref = \%tempCherry;
			#print STDERR "pass $pass, after $seq removal " . scalar(keys %tempCherry)." seqs remain\n";
			### create new expanded counts table i.e. n'(b,l)
			my $tempCounts_ref = countsTable($tempCherryHash_ref);
			#print_table($tempCounts_ref);
			my %tempCounts = %$tempCounts_ref;
			my(@dH, @dHslide);
			### calc dH for each slide
			for(my $i=-1*$opts{s}; $i<=$opts{s}; $i++){
				my $adjustedI = $i - $slides{$seq}[$previousPass];
				#if(abs($adjustedI) <= $opts{s}){
					my @tempdh;
					my $sum = 0;
					### because the sequence array is expanded, adjustments need to be made
					### to the for loop values
					### l = 2,3,4,5,6,7
					for(my $l = $opts{i} + $opts{s}; $l < $opts{i} + $opts{s} + $opts{l}; $l++){
						# k = position of seq to evaluate after shift $i
						#my $k = $l + $i;
						my $k = $l - $adjustedI;
						# d = base identity of position k
						my $d = $tabs[$k];
						#print "For window position $l, \$k = $k\n";
						my $count = $tempCounts{$l}{$d};
						my $H = $dflogf{$count};
						push(@tempdh, $H);
					}
					$sum += $_ for @tempdh;
					#push(@dHslide, $i);
					push(@dHslide, $adjustedI);
					push(@dH, $sum);
					#print STDERR "$seq adjustedI = $adjustedI, dH = $sum\n";
				#}
			}
	
			### find slide position with min dH
			my $min = min(@dH);
			#my $index = first_index {/$min/} @dH;
			#push(my @bestSlide, $dhslide[$index]);
			my @bestSlide; ### hold adjustedI for lowest dH
			for (my $z = 0; $z < scalar(@dH); $z++){
				my $test = $dH[$z];
				if($test == $min){
					push(@bestSlide, $dHslide[$z]);
				}
			}
			my $best = $bestSlide[0];
			if(scalar(@bestSlide) > 1){
				my $r = rand(scalar(@bestSlide));
				$best = $bestSlide[$r];
			}
			#my $adjustedI = $i - $slides{$seq}[$previousPass];
			my $finalSlide = $best + $slides{$seq}[$previousPass];
			push(@{$slides{$seq}}, $finalSlide);
			#print STDERR "$seq pass $pass best slide = $best, previous = $slides{$seq}[$previousPass], finalSlide = $finalSlide\n";
			### create new, full fasta hash with current seq shifted
			if($best ==0){
			}
			elsif($best < 0){
				###shift sequence in @tabs to the left with shift and push
				for(my $z = $best; $z < 0; $z++){
					my $temp = shift(@tabs);
					push(@tabs, $temp);
				}
			}
			else{
				### shift sequence in @tabs to the right with pop and unshift
				for(my $z = $best; $z > 0; $z--){
					my $temp = pop(@tabs);
					unshift(@tabs, $temp);
				}
			}
			###add shuffled sequence back to tempCherry
			$tempCherry{$seq} = join("",@tabs);
			#print STDERR "$seq after pass $pass is " . join("", @tabs) . "\n";
		} ### end of all shuffles, i.e. end of foreach $seq loop
	
		### evaluate info content of shuffled fasta
		#foreach my $seq2 (sort {$a cmp $b} keys %tempCherry){
			##print ">$seq\n$tempCherry{$seq}\n";
		#	print "$pass\t$seq2\t$tempCherry{$seq2}\n";
		#}
		my $shuffledFasta_ref = \%tempCherry;
		my $newInfo = calcInfoContent($shuffledFasta_ref);
		print STDERR "After pass $pass, new information content = $newInfo\n";
		
		$passesHash{$newInfo} = $pass;
	}### end of passes
	return \%cherryFasta
}### end of sub

sub readFasta {
	my(%fasta, $seq, $name);
	my $file = $_[0];
	open(FILE, "< $file") || die "cannot open $file!\n";
	while(my $line = <FILE>){
		chomp($line);
		if($line =~/^>(\S+)/){
			$name = $1;
			$seq = "";
		}
		else{
			$seq = $seq . $line;
			$fasta{$name} = $seq;
			#print "$name\t$seq\n";
		}
	}
	close(FILE);
	
	### create expanded sequence hash
	foreach my $TSS (keys %fasta){
		for(my $i = 0; $i < $opts{s}; $i++){
			$fasta{$TSS} = "-" . $fasta{$TSS} . "-";
		}
	}

	return \%fasta;
}

### this sub used for cherrypicking
sub buildRiw {
	my $hash_ref = $_[0];
	my %fasta = %$hash_ref;
	#my $array_ref = $_[1];
	#my @alpha = @$array_ref;
	my(%counts, %freq, %Riw);
	my $n = scalar(keys %fasta);
	my $En = (1/log(2))*((4-1)/(2*$n));

	foreach my $seq (keys %fasta){
		my @tabs =  split("", $fasta{$seq});
		for(my $i = 2*$opts{s}; $i < 2*$opts{s} + $opts{l}; $i++){
		#for(my $i = $opts{i}; $i < $opts{i} + $opts{l}; $i++){
			$counts{$i}{$tabs[$i]}++;
		}
	}

	### create frequency table
	foreach my $base (@alpha){
		foreach my $pos (sort {$a <=> $b} keys %counts){
			$freq{$pos}{$base} = $counts{$pos}{$base}/$n
		}
	}
	
	### create Riw, Riw keys = 1,2,3,4,5,6
	foreach my $base (@alpha){
		foreach my $pos (sort {$a <=> $b} keys %freq){
			if($freq{$pos}{$base} == 0){
				my $newFreq = 1/($n+2);
				$Riw{$pos}{$base} = 2 - (-1 * (log($newFreq)/log(2)) + $En);
			}
			else{
				$Riw{$pos}{$base} = 2 - (-1 * (log($freq{$pos}{$base})/log(2)) + $En);
			}
		}
	}
	return \%Riw;
}

sub print_table {
	my $hash_ref = $_[0];
	my %table = %$hash_ref;
	foreach my $base (@alpha){
		my @temp;
		foreach my $pos (sort {$a <=> $b} keys %table){
			push(@temp, $table{$pos}{$base});
		}
		my $line = join("\t", @temp);
		print "$base\t$line\n";
	}
}

sub cherrypick{
	my $fasta_ref = $_[0];
	my %cherryFasta = %$fasta_ref;

	my(@bits, $info, $n);
	my $oldCount = scalar(keys %cherryFasta);
	my $test = 0;
	print STDERR "old count = $oldCount, test = $test\n";

	until($oldCount == $test){
		#print "Next pass of cherry picking\n";
		my $Riw_ref = buildRiw($fasta_ref);
		my %Riw = %$Riw_ref;
		foreach my $seq (keys %cherryFasta){
			my @bitArray;
			my $sum = 0;
			my @tabs =  split("", $cherryFasta{$seq});
			###positions will = 2,3,4,5,6,7
			for(my $i = 2*$opts{s}; $i < 2*$opts{s} + $opts{l}; $i++){
				#my $j = $i - $opts{i};
				push(@bitArray, $Riw{$i}{$tabs[$i]});
			}
			$sum += $_ for @bitArray;
			if($sum < 0){
				delete($cherryFasta{$seq});
			}
			else{
				push(@bits, $sum);
				#print "$seq\t$sum\n";
			}
		}
		$n = scalar(keys %cherryFasta);
		print STDERR "$n sequences > 0 bits\n";
		my $sum = 0;
		$sum += $_ for @bits;
		$info = $sum/$n;
		print STDERR "Average information content of remaining sequences = $info\n";
		$test = $oldCount;
		$oldCount = scalar(keys %cherryFasta);	
		print STDERR "after cherrypicking, old count = $oldCount, test = $test, n = $n\n";
		$fasta_ref = \%cherryFasta;
		@bits = ();
	}
	return(\%cherryFasta, $n, $info);
}

###precalc dflogf
sub malignPreCalc{
	my $n = $_[0];
	my %dflogf;
	
	### enter i = 0
	my $val = (-1*((1)/$n)*(log((1/$n))/log(2)));
	$dflogf{0} = $val;

	for(my $i = 1; $i < $n; $i++){
		my $val = (-1*(($i+1)/$n)*(log(($i+1/$n))/log(2))) - (-1*($i/$n)*(log($i/$n)/log(2)));
		$dflogf{$i} = $val;
	}
	return \%dflogf;
}

#countsTable($fastaHash_ref)
#creats countsTable across window of expanded fasta hash 
#windowing parameters take slide spacing into account 
sub countsTable {
	my %counts;
	my $fasta_ref = $_[0];
	my %fasta = %$fasta_ref;
	for(my $i = $opts{i} + $opts{s}; $i < $opts{i} + $opts{l} + $opts{s}; $i++){
		foreach my $b (@alpha){
			$counts{$i}{$b} = 0;
		}
	}
	foreach my $seq (keys %fasta){
		my @tabs =  split("", $fasta{$seq});
		### counts positions = 2,3,4,5,6,7
		for(my $i = $opts{i} + $opts{s}; $i < $opts{i} + $opts{l} + $opts{s}; $i++){
			$counts{$i}{$tabs[$i]}++;
		}
	}
	return \%counts; ###expanded counts table
}

###calculate info content of inputted fasta hash, assumes an expanded hash
### info = calcInfoContent(expandedShiftedFasta_ref, alphaArray_ref)
sub calcInfoContent{
	my $fasta_ref = $_[0];
	my %fasta = %$fasta_ref;
	#my $n = scalar(keys %fasta);
	#my $alpha_ref = $_[1];
	#my @alpha = @$alpha_ref;
	my $counts_ref = countsTable($fasta_ref);
	my %counts = %$counts_ref;
	my(%freq, %Riw, @bits);
	my $n = scalar(keys %fasta);
	my $En = (1/log(2))*((4-1)/(2*$n));
	
	### create frequency table
	foreach my $base (@alpha){
		foreach my $pos (sort {$a <=> $b} keys %counts){
			if($counts{$pos}{$base} == 0){
				my $newFreq = 1/($n+2);
				$freq{$pos}{$base} = $newFreq;
			}
			else{
				$freq{$pos}{$base} = $counts{$pos}{$base}/$n
			}
		}
	}

	### create %Riw
	foreach my $base (@alpha){
		foreach my $pos (sort {$a <=> $b} keys %freq){
			if($freq{$pos}{$base} == 0){
				my $newFreq = 1/($n+2);
				$Riw{$pos}{$base} = 2 - (-1 * (log($newFreq)/log(2)) + $En);
			}
			else{
				$Riw{$pos}{$base} = 2 - (-1 * (log($freq{$pos}{$base})/log(2)) + $En);
			}
		}
	}
	#print "The Riw matrix positions are: ";
	#foreach my $pos (sort {$a <=> $b} keys %freq){
	#	print "$pos, ";
	#}
	#print "\n";

	###calculate per sequence bits
	foreach my $seq (keys %fasta){
		my @bitArray;
		my $sum = 0;
		my @tabs =  split("", $fasta{$seq});
		###calculate bit per position of given sequence
		#for(my $i = $opts{i}; $i < $opts{i} + $opts{l}; $i++){
		### $i = 2,3,4,5,6,7
		for(my $i = $opts{i} + $opts{s}; $i < $opts{i} + $opts{l} + $opts{s}; $i++){
			push(@bitArray, $Riw{$i}{$tabs[$i]});
		}
		### sum bits per position for whole sequence
		$sum += $_ for @bitArray;
		push(@bits, $sum);
	}
	my $n = scalar(keys %fasta);
	#print STDERR "$n sequences > 0 bits\n";
	my $sum = 0;
	$sum += $_ for @bits;
	my $info = $sum/$n;
	return $info;
}
