#! /usr/bin/bash

### this is the master bash wrapper script that will execute the pipeline
### to identify sig70 promoter elements (results written to multiscan.out)
### starting with a set of promoter sequences as input

# Usage statement
[ $# -eq 0 ] && { echo "
Usage: $0 <promoters.fa> <headerString>

where: 
promoters.fa is a fasta file containing the DNA sequence of your promoters
starting with the TSS base and continuing an additional 100bp upstream

headerString is a string that will be appended as a header to various output files; 
currently this argument is required.
"; exit 1; }


python get-10Consensus.py $1 $2 > get-10Consensus.out \
2> get-10Consensus.err

perl parse-10align.pl -f $1 -a $2\.-10malign.txt \
-s $2\.-10malign.seqIDs.txt -T 101 -Z 10 > $2\.prelim35.fa

mv prelim10.spacing.txt $2\.prelim10.spacing.txt
  
perl cherrypick35.pl -i1 -l6 -s1 -p5 $2\.prelim35.fa > $2\.cherrypick35.txt \
2> $2\.cherrypick35.err

sed -i '1s/^/#2,6\n/' $2\.cherrypick35.txt

paste $2\.-10malign.seqIDs.txt $2\.-10malign.txt $2\.prelim10.spacing.txt > $2\.multiscan10.txt

sed -i '1s/^/#8,6\n/' $2\.multiscan10.txt

perl multiscan.pl -1 $2\.multiscan10.txt -3 $2\.cherrypick35.txt -T101 \
$1 > $2\.multiscan.out 2> $2\.multiscan.err
