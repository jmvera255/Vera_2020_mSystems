#! /usr/bin/perl

use strict;
use Getopt::Std;

my (%opts, $ColNum);
getopts("g:m:cA:", \%opts);

my $usage = "*********
Usage: filterBed.pl <options> <File.bed> <limit> > STDOUT
*********
This script will execute various filters or modifications (denoted with the -m option)
to a provided bed file where limit is an integer provided by the user to denote the 
cut off value for the \"m\" options, e.g. -m 1 with limit = 100 will only return
bed features >= 100bp.

Where:
	-m <int> selects what method to apply to bed filter:
		1 = length, greater than: print features >= limit
		2 = add - add length to each end of feature
		3 = add5 - add length only to 5' end of features
		4 = add3 - add length only to 3' end of features
		5 = sub5 - subtract length only from 5' end
		6 = sub3 - subtract length only from 3' end
		7 = test - test error in start stop coordinates, print offending lines
		8 = scoreG - print feature with score >= limit
		9 = scoreL - print feature with score <= limit
		10 = random - print <limit> # of lines at random
		11 = up5 - return region of length <int>  upstream of feature
		12 = down3 - return region of length <int> downstream of feature
		13 = sub - subtract length from each end of feature
		14 = test - test error in start stop coordinates, print non-offending lines
		15 = length, less than: print features < limit
		16 = shift feature in upstream direction by <int> bases
		17 = shift feature in downstream direction by <int> bases
		18 = return the first <int> bases from the 5' end of feature
		19 = return the last <int> bases from the 3' end of feature
		20 = ScoreE - print features with score == <int>
		21 = middle - return the middle <int> bases of each feature
	
	-g <genome.dict> 
		this options is not required but should be use when adding length
		to a bed feature or when selected the up/down stream regions
		so that the new coordinates do not exceed the size of the actuall genome

	-c	Use this option in combination with the -g option for bacterial 
		circular chromosomes or plasmids so that when new coordinates 
		exceed the actual genome coordinates they will wrap/span the genome start/stop
	
	-A <str> Add string to feature name

This script was written by Jessica M. Vera, for questions please contact her.\n\n";

#############_usage statmentes_############
if(@ARGV != 2){
	die "\nInappropriate number of arguments provided!\n\n$usage";
}

if(defined $opts{c} && not defined $opts{g}){
	die "\nYou must provide a genome dictionary when selecting the -c option!\n\n$usage";
}


######_collect @ARGV_######
my $file = shift;
my $length = shift;

###########_print_STDERR_message_##############
my %report = ("6", "User has requested to subtract $length bp from the 3' end",
"12", "User has requested downstream regions only (down3)",
"11", "User has requested upstream regions only (up5)",
"4", "User has requested to add $length bp to the 3' end",
"3", "User has requested to add $length bp to the 5' end",
"2", "User has requested to add $length bp to both ends",
"5", "User has requested to subtract $length bp from the 5' end", 
"13", "User has requested to subtract $length bp from both ends",
"16", "User has requested to shift all features in upstream direction by $length bases",
"17", "User has requested to shift all features in downstream direction by $length bases", 
"18", "User has requested the first $length bp of each feature",
"19", "User has requested the last $length bp of each feature",
"1", "User has selected annotations >= $length bp",
"21", "User has selected to return the middle $length bases of each feature");

my $message = $report{$opts{m}};
print STDERR "$message\n";

#####_parse genome.dict_#####
my(%dict);

if(defined $opts{g}){
	open(DICT, "< $opts{g}") || die "cannot open $opts{g}!\n";
	while(my $line = <DICT>){
		chomp($line);
		my @tabs = split("\t", $line);
		$dict{$tabs[0]} = $tabs[1];
	}
}

######_parse bed file and apply method_#####
my(@lines);
open(FILE, "< $file") || die "Cannot open file $file";
while(my $line = <FILE>){
	chomp($line);
	if($line =~ /^#/){
	}
	else{
		push(@lines, $line);
		my @tabs = split("\t", $line);
		if(defined $opts{A}){
			$tabs[3] = $tabs[3] . $opts{A};
		}
		if($opts{m} == 21){ ### return middle
			my $halfFeature = int((($tabs[2] - $tabs[1])/2) + 0.5) + $tabs[1];
			my $halfLen = int(($length/2) + 0.5);
			$tabs[1] = $halfFeature - $halfLen;
			$tabs[2] = $tabs[1] + $length;
			my $newLine = join("\t", @tabs);
			print "$newLine\n";
		}
		elsif($opts{m} == 1){ ###filter by length
			my $test = ($tabs[2] - $tabs[1]);
			if($test >= $length){
				print "$line\n";
			}
		}
		elsif($opts{m} == 17){ ###shift upstream
			$tabs[1] = $tabs[1] + $length;
			$tabs[2] = $tabs[2] + $length;
			my $newLine = join("\t", @tabs);
			print "$newLine\n";
		}
		elsif($opts{m} == 16){ ###shift upstream
			$tabs[1] = $tabs[1] - $length;
			$tabs[2] = $tabs[2] - $length;
			my $newLine = join("\t", @tabs);
			print "$newLine\n";
		}
		elsif($opts{m} == 15){ ###filter by length
			my $test = ($tabs[2] - $tabs[1]);
			if($test < $length){
				print "$line\n";
			}
		}
		elsif($opts{m} == 13){   ###subtract length from both ends
			$tabs[1] = $tabs[1] + $length;
			$tabs[2] = $tabs[2] - $length;
			my $newLine = join("\t", @tabs);
			print "$newLine\n";
		}		
		elsif($opts{m} == 2){   ###add length to both ends
			if(defined $opts{g}){
				dict(1, $line, $length);
			}
			else{
				if($tabs[1] - $length < 0){
					$tabs[1] = 0;
				}
				else{
					$tabs[1] = $tabs[1] - $length;
				}
				$tabs[2] = $tabs[2] + $length;
				my $newLine = join("\t", @tabs);
				print "$newLine\n";
			}
		}
		elsif($opts{m} == 5){
			if($tabs[5] =~ /-/){  ###subtract from 5' end
				$tabs[2] = $tabs[2] - $length;
			}
			else{
				$tabs[1] =  $tabs[1] + $length;
			}
			my $newLine = join("\t", @tabs);
			print "$newLine\n";
		}
		elsif($opts{m} == 6){   ###subtract from 3' end
			if($tabs[5] =~ /-/){
				$tabs[1] = $tabs[1] + $length;
			}
			else{
				$tabs[2] =  $tabs[2] - $length;
			}
			my $newLine = join("\t", @tabs);
			print "$newLine\n";
		}
		elsif($opts{m} == 3){   ###add to 5' end
			if(defined $opts{g}){
				dict(2,$line,$length);
			}
			else{
				if($tabs[5] =~ /-/){
					$tabs[2] = $tabs[2] + $length;
				}
				else{
					if($tabs[1] < $length){
						$tabs[1] = 0;
					}
					else{
						$tabs[1] =  $tabs[1] - $length;
					}
				}
				my $newLine = join("\t", @tabs);
				print "$newLine\n";
			}
		}
		elsif($opts{m} == 4){   ###add to 3' end
			if(defined $opts{g}){
				dict(5,$line,$length);
			}
			else{
				if($tabs[5] =~ /-/){
					if($tabs[1] < $length){
						$tabs[1] = 0;
					}
					else{
						$tabs[1] = $tabs[1] - $length;
					}
				}
				else{
					$tabs[2] =  $tabs[2] + $length;
				}
				my $newLine = join("\t", @tabs);
				print "$newLine\n";
			}
		}
		elsif($opts{m} == 7){  ###test coordinates
			if($tabs[2] < $tabs[1]){
				print "$line\n";
			}
		}
		elsif($opts{m} == 14){  ###test coordinates
			if($tabs[2] > $tabs[1]){
				print "$line\n";
			}
		}
		elsif($opts{m} == 20){    ###print with score >=
			if($tabs[4] == $length){
				print "$line\n";
			}
		}
		elsif($opts{m} == 8){    ###print with score >=
			if($tabs[4] >= $length){
				print "$line\n";
			}
		}
		elsif($opts{m} == 9){    ####print with score <=
			if($tabs[4] <= $length){
				print "$line\n";
			}
		}
		elsif($opts{m} == 11){    ### return upstream region (up5)
			if(defined $opts{g}){
				dict(3,$line,$length);
			}
			else{
				if($tabs[5] =~ /-/){
					$tabs[1] = $tabs[2];
					$tabs[2] = $tabs[1] + $length;
					my $newLine = join("\t", @tabs);
					print "$newLine\n";
				}
				else{
					$tabs[2] = $tabs[1];
					if($tabs[1] < $length){
						$tabs[1] = 0;
					}
					else{
						$tabs[1] = $tabs[2] - $length;
					}
					my $newLine = join("\t", @tabs);
					print "$newLine\n";
				}
			}
		}
		elsif($opts{m} == 12){   ####return downstream region (down3)
			if(defined $opts{g}){
				dict(4,$line,$length);
			}
			else{
				if($tabs[5] =~ /-/){
					$tabs[2] = $tabs[1];
					if($tabs[1] < $length){
						$tabs[1] = 0;
					}
					else{
						$tabs[1] = $tabs[1] - $length;
					}
					my $newLine = join("\t", @tabs);
					print "$newLine\n";
				}
				else{
					$tabs[1] = $tabs[2];
					$tabs[2] = $tabs[2] + $length;
					my $newLine = join("\t", @tabs);
					print "$newLine\n";
				}
			}
		}
		elsif($opts{m} == 18){
			if($tabs[5] =~ /-/){
				$tabs[1] = $tabs[2] - $length;
				my $newLine = join("\t", @tabs);
				print "$newLine\n";
			}
			else{
				$tabs[2] = $tabs[1] + $length;
				my $newLine = join("\t", @tabs);
				print "$newLine\n";
			}
		}
		elsif($opts{m} == 19){
			if($tabs[5] =~ /-/){
				$tabs[2] = $tabs[1] + $length;
				my $newLine = join("\t", @tabs);
				print "$newLine\n";
			}
			else{
				$tabs[1] = $tabs[2] - $length;
				my $newLine = join("\t", @tabs);
				print "$newLine\n";
			}
		}
	}
}
	
if($opts{m} == 10){   ###print lines are random
	my $a = scalar(@lines);
	$a = $a - 1;
	for(my $i =0; $i < $length; $i++){
		my $t = int(rand($a));
		print "$lines[$t]\n";
	}
}

###########_subroutines_############
sub dict {
	if($_[0] == 1){   ###add length to both ends
		my $length = $_[2];
		my @cols = split("\t", $_[1]); 
		my @alt = @cols;
		if($cols[1] - $length < 0){
			if(defined $opts{c}){
				$alt[1] = $dict{$cols[0]} + ($cols[1] - $length);
				$alt[2] = $dict{$cols[0]};
				my $altLine = join("\t", @alt);
				print "$altLine\n";
			}
			$cols[1] = 0;
		}
		else{
			$cols[1] = $cols[1] - $length;
		}
		if($cols[2] + $length > $dict{$cols[0]}){
			if(defined $opts{c}){
				$alt[1] = 0;
				$alt[2] = $cols[2] + $length - $dict{$cols[0]};
				my $altLine = join("\t", @alt);
				print "$altLine\n";
			}
			$cols[2] = $dict{$cols[0]};
		}
		else{
			$cols[2] = $cols[2] + $length;
		}
		my $newLine = join("\t", @cols);
		print "$newLine\n";
	}
	elsif($_[0] == 2){   ###add to 5' end
		my $length = $_[2];
		my @cols = split("\t", $_[1]); 
		my @alt = @cols;
		if($cols[5] =~ /-/){
			if($cols[2] + $length > $dict{$cols[0]}){
				$cols[2] = $dict{$cols[0]};
			}
			else{
				$cols[2] = $cols[2] + $length;
			}
		}
		else{
			if($cols[1] < $length){
				if(defined $opts{c}){
					$alt[1] = $dict{$cols[0]} + ($cols[1] - $length);
					$alt[2] = $dict{$cols[0]};
					my $altLine = join("\t", @alt);
					print "$altLine\n";
				}
				$cols[1] = 0;
			}
			else{
				$cols[1] =  $cols[1] - $length;
			}
		}
		my $newLine = join("\t", @cols);
		print "$newLine\n";
		$cols[1] = 0;
	}
	elsif($_[0] == 3){   ###up5
		my $length = $_[2];
		my @cols = split("\t", $_[1]); 
		my @alt = @cols;
		if($cols[5] =~ /-/){
			$cols[1] = $cols[2];
			if($cols[2] + $length > $dict{$cols[0]}){
				if(defined $opts{c}){
					$alt[1] = 0;
					$alt[2] = $cols[2] + $length - $dict{$cols[0]};
					my $altLine = join("\t", @alt);
					print "$altLine\n";
				}
				$cols[2] = $dict{$cols[0]};
			}
			else{
				$cols[2] = $cols[1] + $length;
			}
		}
		else{
			$cols[2] = $cols[1];
			if($cols[1] < $length){
				if(defined $opts{c}){
					$alt[1] = $dict{$cols[0]} + ($cols[1] - $length);
					$alt[2] = $dict{$cols[0]};
					my $altLine = join("\t", @alt);
					print "$altLine\n";
				}
				$cols[1] = 0;
			}
			else{
				$cols[1] = $cols[2] - $length;
			}
		}
		my $newLine = join("\t", @cols);
		print "$newLine\n";
	}
	elsif($_[0] == 4){   ##down3
		my $length = $_[2];
		my @cols = split("\t", $_[1]); 
		my @alt = @cols;
		if($cols[5] =~ /-/){
			$cols[2] = $cols[1];
			if($cols[1] < $length){
				if(defined $opts{c}){
					$alt[1] = $dict{$cols[0]} + ($cols[1] - $length);
					$alt[2] = $dict{$cols[0]};
					my $altLine = join("\t", @alt);
					print "$altLine\n";
				}
				$cols[1] = 0;
			}
			else{
				$cols[1] = $cols[1] - $length;
			}
		}
		else{
			$cols[1] = $cols[2];
			if($cols[2] + $length > $dict{$cols[0]}){
				if(defined $opts{c}){
					$alt[1] = 0;
					$alt[2] = $cols[2] + $length - $dict{$cols[0]};
					my $altLine = join("\t", @alt);
					print "$altLine\n";
				}
				$cols[2] = $dict{$cols[0]};
			}
			else{
				$cols[2] = $cols[2] + $length;
			}
		}
		my $newLine = join("\t", @cols);
		print "$newLine\n";
	}
	elsif($_[0] == 5){   ##add to 3' end
		my $length = $_[2];
		my @cols = split("\t", $_[1]);
		my @alt = @cols; 
		if($cols[5] =~ /-/){
			if($cols[1] < $length){
				if(defined $opts{c}){
					$alt[1] = $dict{$cols[0]} + ($cols[1] - $length);
					$alt[2] = $dict{$cols[0]};
					my $altLine = join("\t", @alt);
					print "$altLine\n";
				}
				$cols[1] = 0;
			}
			else{
				$cols[1] = $cols[1] - $length;
			}
		}
		else{
			if($cols[2] + $length > $dict{$cols[0]}){
				$cols[2] = $dict{$cols[0]};
				if(defined $opts{c}){
					$alt[1] = 0;
					$alt[2] = $cols[2] + $length - $dict{$cols[0]};
					my $altLine = join("\t", @alt);
					print "$altLine\n";
				}
			}
			else{
				$cols[2] =  $cols[2] + $length;
			}
		}
		my $newLine = join("\t", @cols);
		print "$newLine\n";
	}
}
