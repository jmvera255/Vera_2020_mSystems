#! /usr/bin/bash

### this is the master bash wrapper script that will execute the pipeline
### to identify sig70 promoter elements (results written to  multiscan.out)
### starting with a set of TSS annotations (in BED format) as input

# Usage statement
[ $# -eq 0 ] && { echo "
Usage: $0 <TSS.bed> <genome.fasta> <headerString>

where: 
TSS.bed is a BED annotation file of single nucleotide TSS positions

genome.fasta is the genome sequence for you species of interest

headerString is a string that will be appended as a header to various output files; 
currently this argument is required.
"; exit 1; }


# pre-process TSS BED files to get promoter sequences
perl parseTSS15bp.pl $1 > temp.bed

perl filterBed.pl -m3 temp.bed 100 > $3.up100.bed

perl /opt/Perl_scripts/bedsfetch.pl $3.bed $2 \
> $3.up100.fa

# peform sigma70 promoter element analysis
python get-10Consensus.py $1 $3 > get-10Consensus.out \
2> get-10Consensus.err

perl parse-10align.pl -f $1 -a $3.-10malign.txt \
-s $3.-10malign.seqIDs.txt -T 101 -Z 10 > $3.prelim35.fa

mv prelim10.spacing.txt $3.prelim10.spacing.txt
  
perl cherrypick35.pl -i1 -l6 -s1 -p5 $3.prelim35.fa > $3.cherrypick35.txt \
2> $3.cherrypick35.err

sed -i '1s/^/#2,6\n/' $3.cherrypick35.txt

paste $3.-10malign.seqIDs.txt $3.-10malign.txt $3.prelim10.spacing.txt > $3.multiscan10.txt

sed -i '1s/^/#8,6\n/' $3.multiscan10.txt

perl multiscan.pl -1 $2\.multiscan10.txt -3 $2\.cherrypick35.txt -T101 $1 \
> $3.multiscan.out 2> $3.multiscan.err
