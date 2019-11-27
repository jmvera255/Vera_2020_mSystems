#! /usr/bin/perl

use strict;
use Getopt::Std;

my (%opts);
getopts("s:f:a:T:Z:A:", \%opts);

my $usage = "********
Usage: parse-10align.pl <options> > STDOUT
********
prints preliminary -35 sequences

where:
	-f <sequences.fa>
	-a <align.txt>
	-s <align.seqIDs.txt>
	-T <int> TSS position in sequences.fa (1-based)
	-Z <int> position in alignment for -10 zero base (i.e. -11 A) (1-based)
	-A <seqID and sequence.txt>
\n";

if(not defined($opts{Z})){
	die "$usage";
}

#### build hash of input sequences from sequences.fa
my(%fastaInput, $seq, $name);
open(FILE, "< $opts{f}") || die "cannot open fasta file: $opts{f}!\n";
while(my $line=<FILE>){
	chomp($line);
	if($line =~/^>(\S+)/){
		$name = $1;
		$seq = "";
	}
	else{
		$seq = $seq . $line;
		$fastaInput{$name} = $seq;
	}
}
close(FILE);

### parse alignment sequence IDs
my(@seqIDs);
if(defined($opts{s})){
	open(FILE, "< $opts{s}") || die "Cannot open seqIDs at $opts{s}!\n";
	while(my $line = <FILE>){
		chomp($line);
		push(@seqIDs, $line);
	}
close(FILE);
}

### parse alignment, get spacing histogram, get -35 start position
my(%histogram, %startsHash);
my $lineCount = 0;
open(OUT, "> prelim10.spacing.txt");
open(FILE, "< $opts{a}") || die "Cannot open alignment file $opts{a}!\n";
while(my $line =<FILE>){
	chomp($line);
	my @sequence = split("", $line);
	my $left = 3;
	my $i = 0;
	while($sequence[$i] =~ /-/){
		$left--;
		$i++;
	}
	my $spacing = scalar(@sequence) - $opts{Z} - $left;
	#print "$seqIDs[$lineCount] length of aligment = " . scalar(@sequence) . ", Z = " . $opts{Z} . ", left = $left, spacing = $spacing\n";
	$histogram{$spacing}++;
	print OUT "$spacing\n";
	#my $start35 = -1*$spacing - 24 - 3;
	my $start35 = -1*$spacing - 24;
	#print $seqIDs[$lineCount] . " -35 start = " . $start35 . "\n";
	$startsHash{$seqIDs[$lineCount]} = $start35;
	$lineCount++;
}

### pull out -35 region from %fastaInput
for(my $i = 0; $i < scalar(@seqIDs); $i++){
	my $id = $seqIDs[$i];
	print ">$seqIDs[$i]\n";
	my @seq = split("", $fastaInput{$id});
	
	### redefine start35 to print out starting 1bp upstream of presumed -35
	my $start = $opts{T} + $startsHash{$id} - 1 - 1;
	#my $start = $opts{T} + $startsHash{$id} - 1;
	### print out an 8mer region
	for(my $j = $start ;$j <= $start + 7; $j++){
	
	### redefine start35 to print out starting 5bp upstream of presumed -35
	#my $start = $opts{T} + $startsHash{$id} - 5 - 1;
	### print out an 15mer region
	#for(my $j = $start ;$j <= $start + 14; $j++){
		#print "j = $j\n";
		print"$seq[$j]";
	}
	print "\n";
}
	

