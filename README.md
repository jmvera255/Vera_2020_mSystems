## Please note this repo is still underconstruction

The scripts available in this repository were used for the work presented in
**\<citation TBD\>**

## Content

The content available in this *repository* is as follows:

#### DataFiles/
>NC_000913.2.fa: _Escherichia coli_ str. K-12 substr. MG1655, complete genome     
>NC_011916.fa: *Caulobacter crescentus* NA1000, complete genome      
>ZM4.fa: compilation of *Z. mobilis* ZM4 chromosome and plasmid sequences    
>Thomason_2015.EcoTSS.bed: *E. coli* TSS annotations used as input for this study[1]   
>Zhou_2015.CcrTSS.bed: *C. crescentus* TSS annotations used as input for this study[2]     
>Zymo.up100.promoters.fa: 3080 *Z. mobilis*   

#### Scripts/
>bedsfetch.pl: convert a BED annotation to a fasta sequence   
>cherrypick35.pl: filter and malign preliminary -35 elements    
>defined_functions.py: library of user defined functions for TBD   
>filterBed.pl: manipulate BED files   
>multiscan.pl: produce final Sigma A/70 model results    
>parse-10align.pl: parse preliminary -10 element results from TBD, return discriminator lengths, generate preliminary set of -35 sequences     
>parseTSS15bp.pl: remove TSS within 15bp of an upstream TSS in the same orientation     





### Citations
[1]: J Bacteriol. 2015 Jan 1;197(1):18-28. doi: [10.1128/JB.02096-14](https://doi.org/10.1128/JB.02096-14).  
[2]: PLoS Genet. 2015 Jan 8;11(1):e1004831. doi: [10.1371/journal.pgen.1004831](https://doi.org/10.1371/journal.pgen.1004831)
