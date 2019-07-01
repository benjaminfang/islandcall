# islandcall

## Introduction
Islandcall is a utility to find genome islands(GIs) whitin bacteria genome. It search GIs which contain a gene(flag gene) by coparison approach.
## Install
### Requirements
* NCBI blast
* nucmer  

### Install  
This utility is writen using Python3.  
>  python3 islandcall.py [-h] [-flag_seq_type {aa,nc}] [-res_dir RES_DIR]  
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
the format is like below:  

```
GENUS -genus name-  
    SPECIES -species name-  
         STRAIN -strain name with path-  
         query  subject  coverage_of_query   identity  e_value   score   query_start query_end   subject_start   subject_end subject_match_sequence  [expand_start:expand_end]|expand_start:expand_end|[expand_start:expand_end] orientation_of_expand_sequence         
```  

* island_calling_res
the format is like below:  

```
GENUS -genus name-  
    SPECIES -species name-
        STRAIN -strain name with path-
        MATCH -match number-
        BLANK_STRAIN -blank stain name with path which expanded sequence compared with-
        query(flag gene)    subject  coverage    identity    subject_start  subject_end expand_seq_location oritation    blank_strain_subject_contig_id    flag_expanded_gene_query_left_match_block    coverage    identity    sub_s   sub_e   expand_seq_position sub_blank   query_left_block    query_right_block   sub_left_block  sub_right_block merge_len   interval_len   island_start(by STRAIN coordinate with merge length)    island_end(by STRAIN coordinate with merge length)  island_sequence(with merge_len added)
```
