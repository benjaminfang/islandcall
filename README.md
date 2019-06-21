# islandcall
## Introduction
Islandcall is a utility to find genome islands(GIs) whitin bacteria genome. It search GIs which contain a gene(flag gene) by coparison approach.
## Install
### Requirements
* NCBI blast
* nucmer
### Install
This utility is writen using Python3.  
>  islandcall.py [-h] [-flag_seq_type {aa,nc}] [-res_dir RES_DIR]  
>                 [-blast_filter_coverage BLAST_FILTER_COVERAGE]  
>                 [-blast_filter_identity BLAST_FILTER_IDENTITY]  
>                 [-expand_len EXPAND_LEN] flag_seq_file dataset  
                 
## usage
* flag_seq_file  
flag gene sequence
* dataset  
data set in which you want to search GIs
* -res_dir
directroy where to put all result
* blast_filter_coverage  
cutoff on coverage of blast result  
* blast_filter_identity  
cutoff on identity of blast result
* expand_len  
legth to expand flag gene match sequence

## result
all result will putinto res_dir. and in this directroy two file will created.  
* flag_gene_blast_res  
* island_calling_res

