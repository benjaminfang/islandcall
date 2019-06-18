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
import fbio.fparse


def get_args():
    args = argparse.ArgumentParser(description='this a utility extracting gene \
                                              island from bacteria genome')
    args.add_argument('flag_seq_file', type=str, help='the fasta \
                       file contain flag gene which identify a gene island.')
    args.add_argument('-flag_seq_type', type=str, default='aa',
                      choices=['aa', 'nc'], help='flag sequence type, untill \
                      now onle amino acid was acceptable.')
    args.add_argument('dataset', type=str, help='dataset in which the island will \
                       will be extracted.')
    args.add_argument('-res_dir', type=str, default='./islandcall_res', help='result\
                        directory in chich all result will stored.')
    args.add_argument('-blast_filter_coverage', type=float, default=0.9, help='blast\
                        result coverage filter factor.')
    args.add_argument('-blast_filter_identity', type=float, default=0.3, help='blast \
                        result identity filter factor.')
    args = args.parse_args()
    flag_seq_file, flag_seq_type, dataset, res_dir = \
        args.flag_seq_file, args.flag_seq_type, args.dataset, args.res_dir
    blast_coverage_cutoff = args.blast_filter_coverage
    blast_identity_cutoff = args.blast_filter_identity
    return (flag_seq_file, flag_seq_type, dataset, res_dir,
            blast_coverage_cutoff, blast_identity_cutoff)


def blast_flag_gene(flag_seq_file, flag_seq_type, dataset, blast_coverage_cutoff,
                    blast_identity_cutoff, res_dir):

    def species_taxonomy(start_dir):
        taxonomy = {}
        for f in os.listdir(start_dir):
            genus = f.split('_')[0]
            if genus not in taxonomy: taxonomy[genus] = {}
            species = f
            if species not in taxonomy[genus]: taxonomy[genus][species] = None
            strain = [os.path.join(start_dir, f, ff) for ff in os.listdir(
                os.path.join(start_dir, f))]
            taxonomy[genus][species] = strain
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

    def expand_seq(fast_dt):
        pass

    def blast_cmd(query,db):
        cmd_name='tblastn'
        query='-query '+query
        db='-db '+db
        outfmt='-outfmt 7'
        out='-out tblastn.tmp_complete_genome'
        num_threads='-num_threads 16'
        cmd=' '.join([cmd_name,query,db,outfmt,out,num_threads])
        os.system(cmd)
        return 'tblastn.tmp_complete_genome'

    def parse_tblastn_res(tblastn_res_f, fasta_f, coverage_cutoff, identity_cutoff):
        expand_len=50000
        data_out=[]
        query_len=400
        fasta_dt=fbio.fparse.Fasta_parse(fasta_f)
        fasta_dt.join_lines()
        blast_dt=[line.rstrip().split() for line in open(tblastn_res_f) if line[0]!='#']
        for record in blast_dt:
            print(record)
            query,subject,identity,alig_len,mismatch,gaps,q_s,q_e,s_s,s_e,e_value,score=record
            #identity=float(identity)/100
            #coverage=int(alig_len)/query_len
            #q_s,q_e,s_s,s_e=int(q_s),int(q_e),int(s_s),int(s_e)
            #if coverage>=coverage_cutoff and identity>=identity_cutoff:
            #    print('>>>',fasta_f,subject,s_s,s_e)
            #    seq,s_s_expand,s_e_expand,seq_expand=expand_seq(fasta_dt,subject,s_s,s_e,expand_len)
            #    data_out.append([query,subject,str(identity),str(coverage),e_value,score,\
            #    str(q_s),str(q_e),str(s_s),str(s_e),str(seq),s_s_expand,s_e_expand,seq_expand])
        return data_out

    if not os.path.exists(res_dir): os.mkdir(res_dir)

    f_out=open(os.path.join(res_dir,'blast_res'),'w')
    taxonomy=species_taxonomy(dataset)
    taxonomy=rm_noncomplete_strain(taxonomy)

    for genus in taxonomy:
        print('GENUS', genus, file=f_out);
        for species in taxonomy[genus]:
            print('    SPECIES',species,file=f_out)
            for strain in taxonomy[genus][species]:
                print('        STRAIN',strain,file=f_out)
                db=os.path.join(strain,os.path.basename(strain)+'.fasta')
                fasta_f=db
                blast_res=blast_cmd(flag_seq_file,db)
                blast_parse_res=parse_tblastn_res(blast_res,fasta_f,blast_coverage_cutoff,blast_identity_cutoff)
                for match in blast_parse_res:
                    #print(match)
                    print('        ','\t'.join(match),file=f_out,sep='')
        break
    f_out.close()

    return 0


def main():
    (flag_seq_file,flag_seq_type,dataset,res_dir,
     blast_coverage_cutoff,blast_identity_cutoff)=get_args()
    blast_flag_gene(flag_seq_file,flag_seq_type,dataset,blast_coverage_cutoff,
                    blast_identity_cutoff,res_dir)
    return 0


if __name__=='__main__':
    res=main()
