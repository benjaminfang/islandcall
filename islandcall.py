#! /usr/bin/env python3
"""
This utility extract island of a flag gene within a interesting database of
bacteria genomes.

flag gene: a gene which consereved in a gene island, and represent the
posibility of the existence of the gene island of such type. utille now a amino
acide flag gene was acceptable.

dataset: a dataset which contain all data which is needed of all process. the
dataset was organised by structure like this:"dataset/genue/species/strain
/seq.fasta,seq.fasta.blastdb". NOTE that, this utility just appliable to
complete genome.

requirement: tblastn, nucmer
"""
import argparse
import os

def get_args():
    args=argparse.ArgumentParser(description='this a utility extracting gene \
                                              island from bacteria genome')
    args.add_argument('flag_seq_file',type=str,help='the fasta \
                       file contain flag gene which identify a gene island.')
    args.add_argument('-flag_seq_type',type=str,default='aa',choices=['aa','nc'], \
                        help='flag sequence type, untill now onle amino acid \
                              was acceptable.')
    args.add_argument('dataset',type=str,help='dataset in which the island will \
                       will be extracted.')
    args.add_argument('-res_dir',type=str,default='./islandcall_res',help='result\
                        directory in chich all result will stored.')
    args.add_argument('-blast_filter_coverage',type=float,default=0.9,help='blast\
                        result coverage filter factor.')
    args.add_argument('-blast_filter_identity',type=float,default=0.3,help='blast \
                        result identity filter factor.')
    args=args.parse_args()
    flag_seq_file,flag_seq_type,dataset,res_dir=args.flag_seq_file,args.flag_seq_type,\
                                        args.dataset,args.res_dir
    blast_coverage_cutoff=args.blast_filter_coverage
    blast_identity_cutoff=args.blast_filter_identity
    return (flag_seq_file,flag_seq_type,dataset,res_dir,
            blast_coverage_cutoff,blast_identity_cutoff)

def blast_flag_gene(flag_seq_file,flag_seq_type,dataset,blast_coverage_cutoff,
                    blast_identity_cutoff,res_dir):
    def species_taxonomy(start_dir):
        taxonomy={}
        for f in os.listdir(start_dir):
            genes=f.split('_')[0]
            if genes not in taxonomy: taxonomy[genes]={}
            species=f
            if species not in taxonomy[genes]: taxonomy[genes][species]=None
            strain=[os.path.join(start_dir,f,ff) for ff in os.listdir(os.path.join(start_dir,f))]
            taxonomy[genes][species]=strain
        return taxonomy

    def rm_noncomplete_strain(data_in):
        data_out={}
        for genus in data_in:
            for species in data_in[genus]:
                for strain in data_in[genus][species]:
                    strain_name=os.path.basename(strain)
                    if strain_name[-9:]=='_complete':
                        if genus not in data_out:
                            data_out[genus]={species:[]}
                        else:
                            if species not in data_out[genus]:
                                data_out[genus][species]=[]
                        data_out[genus][species].append(strain)
        return data_out

    f_out=open(os.path.join(res_dir,'blast_res'),'w')
    taxonomy=species_taxonomy(dataset)
    taxonomy=rm_noncomplete_strain(taxonomy)

    return 0




def main():
    (flag_seq_file,flag_seq_type,dataset,res_dir,
     blast_coverage_cutoff,blast_identity_cutoff)=get_args()
    blast_flag_gene(flag_seq_file,flag_seq_type,dataset,blast_coverage_cutoff,
                    blast_identity_cutoff,res_dir)
    return 0














if __name__=='__main__':
    res=main()
