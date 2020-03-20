#! /usr/bin/perl
####this version of multiscan will remove -35 sequences with Ri < 0bits in cherrypcik35 subroutine 
use strict;
use Math::Trig ':pi';
use DateTime;
my $dt = DateTime->now;
my $start = $dt->hms;
print STDERR "Start time = $start\n";

use List::Util qw( min max);
use List::MoreUtils qw(first_index);
use Getopt::Std;
my (%opts);
getopts("1:3:T:", \%opts);


if(@ARGV != 1){
	die "No file provided!

-1 <file.txt> = tab-delimited file of -10 malign results with seqID-tab-sequence-tab-spacing
	        file must contain a header marked with #<int1>,<int2>; use 0-based
		where <int1> = position of start of -10 window in aligned sequences, include \"-\"
		and <int2> = window length
-3 <file.txt> = file of -35 malign results, no sequence IDs, with header like -10 file
-T <int> = position of TSS in fasta sequences, 1-based\n\n";
}

### define additional global variables
my @alpha = ("A", "G", "T", "C");
my(%malign10, %malign10initial, %fasta, %Riw35initial, %seq10, %hits35);
my(@window10, %flexible, @window35);

### run subroutines
readFasta($ARGV[0]);
read10();
read35();
scan35();
cherrypick35();
Flexible10(); ### Ri(-10) = $flexible{$TSS}[0]
Flexible35(); ### Ri(-35) = $flexible{$TSS}[1]
calcGS();     ### GS = $flexible{$TSS}[2]
evalFlex();

my $dt = DateTime->now;
my $end = $dt->hms;
print STDERR "End time = $end\n";

### subroutines
### calculate total flexibility score for each sequence in %flexibility
### return hits with flex score > 0
### return TSS name, -10 seq, -10 spacing, -35 seq, -35 spacing, Ri(-10), Ri(-35), GS, Ri(total)
sub evalFlex {
	print "#TSS\t-35\tRi(-35)\tSpacing\t-10\tRi(-10)\tGS\tRi(total)\t-10spacing\n";
	foreach my $TSS (keys %flexible){
		my $Ritotal = 0;
		$Ritotal += $_ for @{$flexible{$TSS}};
		if($Ritotal > 0){
			print "$TSS\t$hits35{$TSS}[1]\t$flexible{$TSS}[1]\t$hits35{$TSS}[0]\t";
			my @seq = split("",$seq10{$TSS}[1]);
			my @outSeq;
			for(my $i = $window10[0]; $i < $window10[0] + $window10[1]; $i++){
				push(@outSeq, $seq[$i]);
			}
			my $out10 = join("",@outSeq);
			print "$out10\t$flexible{$TSS}[0]\t$flexible{$TSS}[2]\t$Ritotal\t";
			print "$seq10{$TSS}[0]\n";
		}
	}
}


### use the spacing stored in %hits35 to calculate GS
### add GS to flexible
sub calcGS {
	### %accessibility = n, d = range = 22-27, center = 24
	my (%accessibility, %GS, @Nsum);
	my $N = 0;

	###populate %accessibility across range
	#for (my $d = 22; $d <=27; $d++){
	for (my $d = 21; $d <=26; $d++){
		my $val = (pi2/10.6)*($d - 23);
		my $nd = 1 + cos($val);
		$accessibility{$d} = $nd;
		push(@Nsum, $nd);
	}
	$N += $_ for @Nsum;

	###populate %GS across range
	for (my $d = 21; $d <=26; $d++){
		my $val = $accessibility{$d}/$N;
		my $g = -1 * (log($val)/log(2));
		$GS{$d} = $g;
	}

	###add GS to Flexible
	foreach my $TSS (keys %hits35){
		#my $d = $hits35{$TSS}[0];
		my $d = $hits35{$TSS}[0];
		push(@{$flexible{$TSS}}, $GS{$d});
	}
}

### use the sequences from %hits35 to build a new -35 Riw matrix
### add -35 score to %flexible
sub Flexible35 {
	my(%counts, %freq, %Riw);
	my $n = scalar(keys %hits35);
	foreach my $TSS (keys %hits35){
		my @seq = split("",$hits35{$TSS}[1]);
		for(my $i = 0; $i < 6; $i++){
			$counts{$i}{$seq[$i]}++;
		}
	}
	### create frequency table, keys/pos = 0,1,2,3,4,5
	foreach my $base (@alpha){
		foreach my $pos (sort {$a <=> $b} keys %counts){
			$freq{$pos}{$base} = $counts{$pos}{$base}/$n;
		}
	}
				
	### create final -35 Riw, keys/pos = 0,1,2,3,4,5
	foreach my $base (@alpha){
		foreach my $pos (sort {$a <=> $b} keys %freq){
			if($freq{$pos}{$base} == 0){
				my $newFreq = 1/($n+2);
				$Riw{$pos}{$base} = 2 - (-1 * (log($newFreq)/log(2)));
			}
			else{
				$Riw{$pos}{$base} = 2 - (-1 * (log($freq{$pos}{$base})/log(2)));
			}
		}
	}
	
	### calculate individual -35 score, add to flexible
	foreach my $TSS (keys %hits35){
		my $sum = 0;
		my @bits;
		my @seq = split("",$hits35{$TSS}[1]);
		for(my $i = 0; $i < 6; $i++){
			push(@bits, $Riw{$i}{$seq[$i]});
		}
		$sum += $_ for @bits;
		push(@{$flexible{$TSS}}, $sum);
	}
}

sub Flexible10 {
	my(%counts, %freq, %Riw);
	my $n = scalar(keys %hits35);
	foreach my $TSS (keys %hits35){
		my @seq = split("",$seq10{$TSS}[1]);
		for(my $i = $window10[0]; $i < $window10[0] + $window10[1]; $i++){
			$counts{$i}{$seq[$i]}++;
		}
	}
	### create frequency table
	foreach my $base (@alpha){
		foreach my $pos (sort {$a <=> $b} keys %counts){
			$freq{$pos}{$base} = $counts{$pos}{$base}/$n;
		}
	}
				
	### create -10 Riw
	foreach my $base (@alpha){
		foreach my $pos (sort {$a <=> $b} keys %freq){
			if($freq{$pos}{$base} == 0){
				my $newFreq = 1/($n+2);
				$Riw{$pos}{$base} = 2 - (-1 * (log($newFreq)/log(2)));
			}
			else{
				$Riw{$pos}{$base} = 2 - (-1 * (log($freq{$pos}{$base})/log(2)));
			}
		}
	}
	
	### calculate individual -10 score, add to flexible
	foreach my $TSS (keys %hits35){
		my $sum = 0;
		my @bits;
		my @seq = split("",$seq10{$TSS}[1]);
		for(my $i = $window10[0]; $i < $window10[0] + $window10[1]; $i++){
			push(@bits, $Riw{$i}{$seq[$i]});
		}
		$sum += $_ for @bits;
		push(@{$flexible{$TSS}}, $sum);
	}
}

### use Riw35initial to scan for best >0 bit match in list of -10 sequences
### scans with score > 0 are stored in $hits{sequenceName} = spacing
sub scan35 {
	foreach my $s (keys %seq10){
		my @sequence = split("", $fasta{$s});
		my(@bitsArray, @spacing);
		for(my $i = 21; $i<=26; $i++){ ### Schneider -35 range is 21-26 but only to 2nd base in 35 6mer
		### so adding 1 to this spacing gets me to the 1st base in the 35 6mer
			my(@bits);
			my $sum = 0;
			for(my $j = 0; $j < 6; $j++){
				my $pos = $opts{T} - $seq10{$s}[0] - $i + $j - 2;
				my $k = $j + $window35[0];  ###k = Riw35 postion
				push(@bits, $Riw35initial{$k}{$sequence[$pos]});
			}
			$sum += $_ for @bits;
			push(@bitsArray, $sum);
			push(@spacing, $i);
		}

		###find spacing with best -35 score
		my $max = max(@bitsArray);
		my @bestScan;
		for (my $z = 0; $z < scalar(@bitsArray); $z++){
			my $test = $bitsArray[$z];
			if($test == $max){
				push(@bestScan, $spacing[$z]);
			}
		}
		my $best = $bestScan[0];
		if(scalar(@bestScan) > 1){
			my $r = rand(scalar(@bestScan));
			$best = $bestScan[$r];
		}

		if($max > 0){
			my(@best35);
			for(my $j = 0; $j < 6; $j++){
				my $pos = $opts{T} - $seq10{$s}[0] - $best + $j - 2;
				push(@best35, $sequence[$pos]);
			}
			my $out35 = join("", @best35);

			#push(@{$hits35{$s}}, $max, $best);
			#$hits35{$s} = $best; ###$hit35{TSS} = best hit spacing 
			push(@{$hits35{$s}}, $best, $out35); ### $hits35{TSS} -> [0] = spacing, [1] = -35 seq
		}
	}
}

### parse -35 malign results
sub read35 {
	
	my(%counts, %freq);

	### parse -35 alignment file and create counts
	open(FILE, " < $opts{3}" ) || die "Cannot open -35 malign file $opts{3}!!\n";
	my $lineCount = 0;
	while(my $line = <FILE>){
		chomp($line);
		if($line =~ /^#/){
			if($lineCount == 0){
				$lineCount++;
				$line =~ /^#(\d+),(\d+)/;
				push(@window35, $1, $2);
			}
		}
		else{
			if($lineCount != 0){
				$lineCount++;
				my @seq = split("", $line);
				### count/freq/Riw35 postions = window35[0] to window[0] +[1]
				for(my $i = $window35[0]; $i < $window35[0] + $window35[1]; $i++){
					$counts{$i}{$seq[$i]}++;
				}
			}
			else{
				die "No header in -35 malign file! File must include header\n";
			}
		}
	}
	### print counts
	foreach my $base (@alpha){
		my @temp;
		foreach my $pos (sort {$a <=> $b} keys %counts){
		#	print STDERR "$pos\t";
			push(@temp, $counts{$pos}{$base});
		}
		my $line = join("\t", @temp);
	}
	### create frequency table
	foreach my $base (@alpha){
		foreach my $pos (sort {$a <=> $b} keys %counts){
			$freq{$pos}{$base} = $counts{$pos}{$base}/$lineCount;
		}
	}
				
	### create Riw35initial
	foreach my $base (@alpha){
		foreach my $pos (sort {$a <=> $b} keys %freq){
			if($freq{$pos}{$base} == 0){
				my $newFreq = 1/($lineCount+2);
				$Riw35initial{$pos}{$base} = 2 - (-1 * (log($newFreq)/log(2)));
			}
			else{
				$Riw35initial{$pos}{$base} = 2 - (-1 * (log($freq{$pos}{$base})/log(2)));
			}
		}
	}
}


sub read10 {
	
	my(%counts);

	### parse -10 alignment file and create counts
	open(FILE, " < $opts{1}" ) || die "Cannot open -10 malign file $opts{1}!!\n";
	my $lineCount = 0;
	while(my $line = <FILE>){
		chomp($line);
		if($line =~ /^#/){
			if($lineCount==0){
				$line =~ /^#(\d+),(\d+)/;
				push(@window10, $1, $2);
				$lineCount++;
			}
		}
		else{
			if($lineCount!=0){
				$lineCount++;
				my @tabs = split("\t", $line);
				push(@{$seq10{$tabs[0]}}, $tabs[2], $tabs[1]); ###[0] = spacing, [1] = sequence
			}
			else{
				die "No header in -10 malign file! File must include header\n";
			}
		}
	}
}



sub readFasta {
	my($seq, $name);
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
		}
	}
	close(FILE);
	return \%fasta;
}

sub cherrypick35{

	my(@bits, $info, $n);
	my $oldCount = scalar(keys %hits35);
	my $test = 0;
	print STDERR "old count = $oldCount, test = $test\n";

	until($oldCount == $test){
		my $hits35_ref = \%hits35;
		my $Riw_ref = buildRiw($hits35_ref);
		my %Riw = %$Riw_ref;
		foreach my $seq (keys %hits35){
			my @bitArray;
			my $sum = 0;
			my @tabs =  split("", $hits35{$seq}[1]);
			for(my $i = 0; $i < $window35[1]; $i++){
				push(@bitArray, $Riw{$i}{$tabs[$i]});
			}
			$sum += $_ for @bitArray;
			if($sum < 0){
				delete($hits35{$seq});
			}
			else{
				push(@bits, $sum);
			}
		}
		$n = scalar(keys %hits35);
		print STDERR "$n sequences > 0 bits\n";
		my $sum = 0;
		$sum += $_ for @bits;
		$info = $sum/$n;
		print STDERR "Average information content of remaining sequences = $info\n";
		$test = $oldCount;
		$oldCount = scalar(keys %hits35);	
		print STDERR "after cherrypicking, old count = $oldCount, test = $test, n = $n\n";
		@bits = ();
	}
}

sub buildRiw {
	my $hash_ref = $_[0];
	my %fasta = %$hash_ref;
	#my $array_ref = $_[1];
	#my @alpha = @$array_ref;
	my(%counts, %freq, %Riw);
	my $n = scalar(keys %fasta);
	
	foreach my $seq (keys %fasta){
		my @tabs =  split("", $fasta{$seq}[1]);
		for(my $i = 0; $i < $window35[1]; $i++){
			$counts{$i}{$tabs[$i]}++;
		}
	}

	### create frequency table
	foreach my $base (@alpha){
		foreach my $pos (sort {$a <=> $b} keys %counts){
			$freq{$pos}{$base} = $counts{$pos}{$base}/$n
		}
	}

	### create Riw, Riw keys = 0,1,2,3,4,5
	foreach my $base (@alpha){
		foreach my $pos (sort {$a <=> $b} keys %freq){
			if($freq{$pos}{$base} == 0){
				my $newFreq = 1/($n+2);
				$Riw{$pos}{$base} = 2 - (-1 * (log($newFreq)/log(2)));
			}
			else{
				$Riw{$pos}{$base} = 2 - (-1 * (log($freq{$pos}{$base})/log(2)));
			}
		}
	}
	return \%Riw;
}
