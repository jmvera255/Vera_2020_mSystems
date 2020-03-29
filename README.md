## Genome-scale transcription–translation mapping reveals features of Zymomonas mobilis transcription units and promoters

The scripts and files available in this repository were used for the work presented in
**\<citation TBD\>**

## Content

#### DataFiles/
>NC_000913.2.fa: _Escherichia coli_ str. K-12 substr. MG1655, complete genome     
>NC_011916.fa: *Caulobacter crescentus* NA1000, complete genome      
>ZM4.fa: compilation of *Z. mobilis* ZM4 chromosome and plasmid sequences [1]    
>Thomason_2015.EcoTSS.bed: *E. coli* TSS annotations [2]   
>Zhou_2015.CcrTSS.bed: *C. crescentus* TSS annotations [3]     
>Zymo.up100.promoters.fa: 3080 *Z. mobilis* TSS annotations from this study 

#### Scripts/
>bedsfetch.pl: convert a BED annotation to a fasta sequence   
>cherrypick35.pl: filter and malign preliminary -35 elements    
>defined_functions.py: library of user defined functions for get-10Consensus.py      
>get-10Consensus.py: perform malign on inputted promoter sequences, return preliminary -10 element results      
>filterBed.pl: manipulate BED files   
>multiscan.pl: produce final Sigma A/70 model results    
>parse-10align.pl: parse preliminary -10 element results from get-10Consensus.py, return discriminator lengths, generate preliminary set of -35 sequences     
>parseTSS15bp.pl: remove TSS within 15bp of an upstream TSS in the same orientation     
>get_sig70_from_fasta.sh: wrapper script for pipeline starting with 100bp promoter sequences    
>get_sig70_from_TSS.sh: warpper script for pipeline starting with TSS annotations (used for *E. coli* and *C. crescentus* analysis)

### Citations
[1] Biotechnol Biofuels. 2018 May 2;11:125. doi:
[10.1186/s13068-018-1116-x](https://biotechnologyforbiofuels.biomedcentral.com/articles/10.1186/s13068-018-1116-x)   
[2] J Bacteriol. 2015 Jan 1;197(1):18-28. doi: [10.1128/JB.02096-14](https://doi.org/10.1128/JB.02096-14)  
[3] PLoS Genet. 2015 Jan 8;11(1):e1004831. doi: [10.1371/journal.pgen.1004831](https://doi.org/10.1371/journal.pgen.1004831)   
