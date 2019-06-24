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

Two result file will generated after the finish of this utility.
the fist: flag_gene_blast_res.
the second: island_calling_res.
"""

import argparse
import os
import re
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
    args.add_argument('-expand_len', type=int, default=50000, help='legth when expand \
        flag match seq.')
    args = args.parse_args()
    flag_seq_file, flag_seq_type, dataset, res_dir = \
        args.flag_seq_file, args.flag_seq_type, args.dataset, args.res_dir
    blast_coverage_cutoff = args.blast_filter_coverage
    blast_identity_cutoff = args.blast_filter_identity
    expand_len = args.expand_len
    return (flag_seq_file, flag_seq_type, dataset, res_dir,
            blast_coverage_cutoff, blast_identity_cutoff, expand_len)


def blast_flag_gene(flag_seq_file, flag_seq_type, dataset, blast_coverage_cutoff,
                    blast_identity_cutoff, res_dir, expand_len):

    def species_taxonomy(start_dir):
        taxonomy = {}
        for f in os.listdir(start_dir):
            genus = f.split('_')[0]
            if genus not in taxonomy: taxonomy[genus] = {}
            species = f
            if species not in taxonomy[genus]: taxonomy[genus][species] = None
            strain = [os.path.join(os.path.abspath(start_dir), f, ff) for ff in os.listdir(
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

    def blast_cmd(query, db, blast_tmp_f):
        cmd_name='tblastn'
        query='-query '+query
        db='-db '+db
        outfmt='-outfmt 7'
        out='-out '+blast_tmp_f
        num_threads='-num_threads 16'
        cmd=' '.join([cmd_name,query,db,outfmt,out,num_threads])
        os.system(cmd)
        return 0

    def extract_seq(start, end, oritation, seq):
        base_dic={'A':'T', 'T':'A', 'G':'C', 'C':'G', 'a':'t', 't':'a', 'g':'c', 'c':'g'}
        if oritation == '+':
            if start-end<=0:
                return seq[start-1:end]
            else:
                raise Exception('start should less than or equal end.')
        elif oritation == '-':
            if start-end>=0:
                return ''.join([base_dic.get(b,'N') for b in seq[end-1:start][::-1]])
        else:
            raise Exception('oritation should marked as "+" or "-".')

    def get_expanded_seq(seq_expand_range, oritation, seq):
        expand_seq = ''
        for piece in seq_expand_range:
            if None not in piece:
                expand_seq += extract_seq(piece[0], piece[1], oritation, seq)
        return expand_seq

    def expand_seq_complete_genome(fast_dt, subject, s_s, s_e, expand_len):
        seq_len = None
        seq_expand_range = [[None,None],[None,None],[None,None]]
        oritation = None
        seq_subject = None
        seq_expanded = None
        for seqhead in fast_dt.data:
            if seqhead.split()[0]==subject:
                seq=fast_dt.data[seqhead]
        seq_len = len(seq)
        if s_s<s_e:
            oritation = '+'
            seq_subject = extract_seq(s_s,s_e,oritation,seq)
            if s_s-expand_len<=0:
                seq_expand_range[0][0] = seq_len+(s_s-expand_len)
                seq_expand_range[0][1] = seq_len
                seq_expand_range[1][0] = 1
            else:
                seq_expand_range[1][0] = s_s-expand_len
            if s_e+expand_len>seq_len:
                seq_expand_range[1][1] = seq_len
                seq_expand_range[2][0] = 1
                seq_expand_range[2][1] = s_e+expand_len-seq_len
            else:
                seq_expand_range[1][1] = s_e+expand_len
        elif s_s>s_e:
            oritation = '-'
            seq_subject = extract_seq(s_s,s_e,oritation,seq)
            if s_s+expand_len>seq_len:
                seq_expand_range[0][0] = s_s+expand_len-seq_len
                seq_expand_range[0][1] = 1
                seq_expand_range[1][0] = seq_len
            else:
                seq_expand_range[1][0] = s_s+expand_len
            if s_e-expand_len<=0:
                seq_expand_range[1][1] = 1
                seq_expand_range[2][0] = seq_len
                seq_expand_range[2][1] = seq_len+(s_e-expand_len)
            else:
                seq_expand_range[1][1] = s_e-expand_len
        else:
            raise Exception('opps...s_s can not equal s_e. under blast result context')
        seq_expanded=get_expanded_seq(seq_expand_range, oritation, seq)
        return seq_subject, seq_expand_range, oritation, seq_expanded

    def parse_tblastn_res(tblastn_res_f, fasta_f, coverage_cutoff, identity_cutoff, expand_len):
        #expand_len = 50000
        query_len = 400
        data_out = []
        fasta_dt = fbio.fparse.Fasta_parse(fasta_f)
        fasta_dt.join_lines()
        blast_dt = [line.rstrip().split() for line in open(tblastn_res_f) if line[0]!='#']
        for record in blast_dt:
            print('       ',record)
            query,subject,identity,alig_len,mismatch,gaps,q_s,q_e,s_s,s_e,e_value,score=record
            identity = round(float(identity)/100,5)
            coverage = int(alig_len)/query_len
            if coverage>1: coverage=1
            q_s,q_e,s_s,s_e=int(q_s),int(q_e),int(s_s),int(s_e)
            if coverage>=coverage_cutoff and identity>=identity_cutoff:
                # DEBUG:
                print('        ^ This is ok')
                seq_subject, seq_expand_range, oritation, seq_expanded= \
                    expand_seq_complete_genome(fasta_dt, subject, s_s, s_e, expand_len)
                seq_expand_range_str=[]
                for piece in seq_expand_range:
                    if None in piece:
                        seq_expand_range_str.append(':')
                    else:
                        seq_expand_range_str.append( str(piece[0]) + ':' + str(piece[1]) )
                seq_expand_range_str = '|'.join(seq_expand_range_str)

                data_out.append([query, subject, str(coverage), str(identity), \
                    e_value, score, str(q_s), str(q_e), str(s_s), str(s_e), \
                    str(seq_subject), seq_expand_range_str, oritation, seq_expanded])
        return data_out

    if not os.path.exists(res_dir): os.mkdir(res_dir)
    f_out_file_name=os.path.join(res_dir, 'flag_gene_blast_res')
    blast_tmp_f=os.path.join(res_dir, 'blast_tmp_file')
    f_out=open(f_out_file_name, 'w')
    taxonomy=species_taxonomy(dataset)
    taxonomy=rm_noncomplete_strain(taxonomy)
    print('+++++++++++blasting++++++++++')
    for genus in taxonomy:
        print('>', genus)
        print('GENUS', genus, file=f_out);
        for species in taxonomy[genus]:
            print('    >>', species)
            print('    SPECIES',species,file=f_out)
            for strain in taxonomy[genus][species]:
                print('        >>>',strain)
                print('        STRAIN',strain,file=f_out)
                db = os.path.join(strain, os.path.basename(strain)+'.fasta')
                fasta_f = db
                blast_cmd(flag_seq_file, db, blast_tmp_f)
                blast_parse_res = parse_tblastn_res(blast_tmp_f, fasta_f, \
                    blast_coverage_cutoff, blast_identity_cutoff, expand_len)
                for match in blast_parse_res:
                    print('        ','\t'.join(match), file=f_out, sep='')
    f_out.close()
    return f_out_file_name


def callisland(blast_res_file, res_dir, expand_len):

    def struc_blast_res(blast_res_file):
        data_out={}
        re1=re.compile('^GENUS\s')
        re2=re.compile('^\s{4}SPECIES\s')
        re3=re.compile('^\s{8}STRAIN\s')
        for line in open(blast_res_file):
            line=line.rstrip()
            if re1.match(line):
                genus=line.split()[1]
            elif re2.match(line):
                species=line.split()[1]
            elif re3.match(line):
                strain=line.split()[1]
                if genus not in data_out:
                    data_out[genus]={species:{}}
                elif species not in data_out[genus]:
                    data_out[genus][species]={}
                data_out[genus][species][strain]=[]
            else:
                data_out[genus][species][strain].append(line)
        return data_out

    def filter_empty_genus(data_in):
        data_out = {}
        for genus in data_in:
            i = 0
            for species in data_in[genus]:
                for strain in data_in[genus][species]:
                    if len(data_in[genus][species][strain])>0:
                        i += 1
            if i>0:
                data_out[genus] = data_in[genus]
        return data_out

    def split_strain_have_flag(data_in):
        data_out = {'Y':{},'N':{}}
        for genus in data_in:
            genus_Y = {}
            genus_N = {}
            for species in data_in[genus]:
                for strain in data_in[genus][species]:
                    if len(data_in[genus][species][strain])==0:
                        if species not in genus_N: genus_N[species]={}
                        genus_N[species][strain] = None
                    else:
                        if species not in genus_Y: genus_Y[species]={}
                        genus_Y[species][strain] = data_in[genus][species][strain]
            if len(genus_Y)>0 and len(genus_N)>0:
                data_out['Y'][genus] = genus_Y
                data_out['N'][genus] = genus_N
        return data_out

    def get_species_N_strain(genus, species, blast_dt_splited):
        data_out=[]
        if species in blast_dt_splited['N'][genus]:
            for strain in blast_dt_splited['N'][genus][species]:
                #print(strain)
                data_out.append(strain)
        return data_out

    def get_genus_N_strain(genus, blast_dt_splited):
        data_out=[]
        for species in blast_dt_splited['N'][genus]:
            for strain in blast_dt_splited['N'][genus][species]:
                data_out.append(strain)
        return data_out

    def run_nucmer(ref_file, query_file, out_file):
        cmd_name = 'nucmer'
        delta = '--delta ' + out_file
        cmd = ' '.join([cmd_name, delta, ref_file, query_file])
        print('        nucmer command:',cmd)
        os.system(cmd)
        return out_file

    def struc_nucmer_res(nucmer_res, expand_len):

        def parse_block(block):
            new_block = []
            alig_distribution = []
            alig_distribution_new = []
            for line in block:
                if len(line)>1:
                    new_block.append(line)
            ref_seq_name, ref_seq_len, query_seq_len = new_block[0]
            ref_seq_name = ref_seq_name[2:]
            ref_seq_len = int(ref_seq_len)
            query_seq_len = int(query_seq_len)
            for line in new_block[1:]:
                alig_distribution.append([(int(line[0]), int(line[1])), \
                    (int(line[2]), int(line[3]))])
            for alig in alig_distribution:
                ref_s, ref_e = alig[0]
                query_s, query_e = alig[1]
                if query_s < query_e:
                    match_consistance = 'cis'
                else:
                    match_consistance = 'trans'
                    query_s, query_e = query_e, query_s
                    ref_s, ref_e = ref_e, ref_s
                alig_distribution_new.append([(query_s, query_e), \
                    (ref_s, ref_e), query_e-query_s, match_consistance])
            alig_distribution_new.sort(key=lambda x:x[0][0])

            return ref_seq_name, ref_seq_len, query_seq_len,alig_distribution_new

        data_out = {}
        a_line = [line.rstrip().split() for line in open(nucmer_res)]
        block = []
        i = 0
        j = 0
        for line in a_line:
            if line[0][0]=='>':
                if i==0: i = 1
                if j==1:
                    ref_seq_name, ref_seq_len, query_seq_len, alig_distribution_sorted = \
                        parse_block(block)
                    flag_gene_position = [expand_len+1, query_seq_len-expand_len]
                    data_out[ref_seq_name] = {'ref_seq_len':       ref_seq_len,
                                              'query_seq_len':     query_seq_len,
                                              'flag_gene_position':flag_gene_position,
                                              'alig_distribution_sorted':alig_distribution_sorted}
                if j==0: j = 1
                block = []
            if i==1:
                block.append(line)
        if i==1:
            ref_seq_name, ref_seq_len, query_seq_len, alig_distribution_sorted = \
                parse_block(block)
            flag_gene_position = [expand_len+1, query_seq_len-expand_len]
            data_out[ref_seq_name] = {'ref_seq_len':       ref_seq_len,
                                      'query_seq_len':     query_seq_len,
                                      'flag_gene_position':flag_gene_position,
                                      'alig_distribution_sorted':alig_distribution_sorted}
        return data_out

    def find_flag_left_right_match(nucmer_dt, match_block_distance):

        def judge_island(left_block, right_block):
            merge_len = None
            interval_len = None
            have_island = None
            left_block_oritation = left_block[3]
            right_block_oritaion = right_block[3]
            left_range = left_block[1]
            left_range.sort()
            right_range = right_block[1]
            right_range.sort()
            sequence_pattren = [(left_range[0],'ll'),(left_range[1],'lr'), \
                (right_range[0],'rl'),(right_range[1],'rr')]
            sequence_pattren.sort(key=lambda x:x[0])
            sequence_pattren = [ele[1] for ele in sequence_pattren]
            print('        suject sequence_pattren:',sequence_pattren)
            # DEBUG: Here bugs need to fix. because different suject block oritation
            # DEBUG: cis, cis; cis, trans; trans, cis; trans, trans.
            if left_block_oritation == 'cis' and right_block_oritaion == 'trans':
                have_island = 'N'
            elif left_block_oritation == 'trans' and right_block_oritaion == 'cis':
                have_island = 'N'
            else:
                if sequence_pattren==['ll', 'lr', 'rl', 'rr']:
                    interval_len = right_range[0] - left_range[1]
                    if interval_len < match_block_distance:
                        have_island = 'Y'
                    else:
                        have_island = 'N'
                elif sequence_pattren == ['ll', 'rl', 'lr', 'rr']:
                    merge_len = left_range[1] - right_range[0]
                    if merge_len > 100: merge_len = 100
                    have_island = 'Y'
                elif sequence_pattren == ['rl', 'll', 'lr', 'rr']:
                    if (left_range[1] - right_range[1]) < 100 or (right_range[0] - left_range[0]) < 100:
                        have_island = 'Y'
                    else:
                        have_island = 'N'
                elif sequence_pattren == ['ll', 'rl', 'rr', 'lr']:
                    if (right_range[1] - left_range[1]) < 100 or (left_range[0] - right_range[0]) < 100:
                        have_island = 'Y'
                    else:
                        have_island = 'N'
                elif sequence_pattren == ['rl', 'll', 'rr', 'lr']:
                    have_island = 'N'
                elif sequence_pattren == ['rl', 'rr', 'll', 'lr']:
                    have_island = 'N'
            # DEBUG: print(have_island, merge_len, interval_len)
            return have_island, merge_len, interval_len

        data_out = {}
        match_block_min_len = 500
        for ref_contig in nucmer_dt:
            flag_left = nucmer_dt[ref_contig]['flag_gene_position'][0]
            flag_right = nucmer_dt[ref_contig]['flag_gene_position'][1]
            query_match_series = {}
            for ele in nucmer_dt[ref_contig]['alig_distribution_sorted']:
                # DEBUG: print(ele[2])
                if ele[2]>match_block_min_len:
                    query_match_series[ele[0]]=ele
            # DEBUG: print(query_match_series)
            # DEBUG: print(flag_left,flag_right)
            tmp_left = []
            tmp_right = []
            for ele in query_match_series:
                if ele[1] < flag_left:
                    tmp_left.append(ele)
                if ele[0] > flag_right:
                    tmp_right.append(ele)
            # DEBUG: print(tmp_left)
            # DEBUG: print(tmp_right)
            if len(tmp_left)>0:
                tmp_left.sort(key=lambda x:x[1])
                left_block = query_match_series[tmp_left[-1]]
            else:
                left_block = None
            if len(tmp_right)>0:
                tmp_right.sort(key=lambda x:x[0])
                right_block = query_match_series[tmp_right[0]]
            else:
                right_block = None
            if left_block and right_block:
                print('       left and right block:',left_block,right_block)
                have_island, merge_len, interval_len = judge_island(left_block, right_block)
                if have_island == 'Y':
                    data_out[ref_contig] = [nucmer_dt[ref_contig]['flag_gene_position'], \
                        left_block[0], right_block[0], left_block[1], right_block[1], \
                        merge_len, interval_len]

        if len(data_out)>0:
            return 'Y', data_out
        else:
            return 'N', data_out

    def parse_nucmer_res(nucmer_res, ref_file, query_file, expand_len):
        ref_dt = fbio.fparse.Fasta_parse(ref_file)
        ref_dt.join_lines()
        query_dt = fbio.fparse.Fasta_parse(query_file)
        query_dt.join_lines()

        nucmer_dt = struc_nucmer_res(nucmer_res, expand_len)
        # DEBUG: for head in nucmer_dt:
        # DEBUG:     print(head)
        print('        nucmer_res_dt:',nucmer_dt)
        match_block_distance=100
        if_have_island, island_info = find_flag_left_right_match(nucmer_dt, \
            match_block_distance)

        return if_have_island, island_info

    def print_island_info(query, suject, s_s, s_e, island_info, subject_expanded_seq, file_out):
        for island_contig in island_info:
            flag_gene_position_start, flag_gene_position_end, \
                query_left_block, query_right_block, sub_leaf_block, sub_right_block, \
                merge_len, interval_len = island_info[island_contig][0][0], \
                island_info[island_contig][0][1], island_info[island_contig][1], \
                island_info[island_contig][2], island_info[island_contig][3], island_info[island_contig][4], \
                island_info[island_contig][5], island_info[island_contig][6]
            island_left = query_left_block[1]
            island_right = query_right_block[0]
            if merge_len:
                island_left = island_left - merge_len
                island_right = island_right + merge_len
            island_seq = subject_expanded_seq[island_left:island_right]
            # trans coordinate to int contain strain seq.
            island_left_traned = int(s_s) - (flag_gene_position_start - island_left)
            island_right_traned = int(s_e) + (island_right - flag_gene_position_end)

            print('       ', '\t'.join([query, subject, s_s, s_e, island_contig,
                str(island_left_traned),str(island_right_traned),
                ':'.join([str(query_left_block[0]), str(query_left_block[1])]),
                ':'.join([str(query_right_block[0]), str(query_right_block[1])]),
                ':'.join([str(sub_leaf_block[0]), str(sub_leaf_block[1])]),
                ':'.join([str(sub_right_block[0]), str(sub_right_block[1])]),
                 str(merge_len), str(interval_len), island_seq]), file=file_out)

        return 0

    print('\n++++++++++islandcalling+++++++++')
    blast_dt = struc_blast_res(blast_res_file)
    # DEBUG: print(blast_dt)
    blast_dt = filter_empty_genus(blast_dt)
    blast_dt_splited = split_strain_have_flag(blast_dt)

    # DEBUG:>-----------------
    debug_blast_dt={}
    for c in blast_dt_splited:
        debug_blast_dt[c]={}
        for genus in blast_dt_splited[c]:
            debug_blast_dt[c][genus]={}
            for species in blast_dt_splited[c][genus]:
                debug_blast_dt[c][genus][species]={}
                for strain in blast_dt_splited[c][genus][species]:
                    if blast_dt_splited[c][genus][species][strain]:
                        debug_blast_dt[c][genus][species][os.path.basename(strain)] = 'Y'
                    else:
                        debug_blast_dt[c][genus][species][os.path.basename(strain)] = 'N'
    print(debug_blast_dt)
    # DEBUG: <----------------

    nucmer_out_file = os.path.join(res_dir, 'nucmer_res')
    flag_segment_file = os.path.join(res_dir, 'flag_gene_segment.fasta')
    island_res_file=os.path.join(res_dir, 'islandcall_calling_res')
    island_res_f_out = open(island_res_file, 'w')
    for genus in blast_dt_splited['Y']:
        print('>',genus)
        print('GENUS', genus, file=island_res_f_out)
        genus_N_strain=get_genus_N_strain(genus, blast_dt_splited)
        for species in blast_dt_splited['Y'][genus]:
            print('    >>',species)
            print('    SPECIES', species, file=island_res_f_out)
            species_N_strain = get_species_N_strain(genus, species, blast_dt_splited)
            genus_N_strain = set(genus_N_strain) - set(species_N_strain)
            # DEBUG: print('genus_N_strain', genus_N_strain)
            for strain in blast_dt_splited['Y'][genus][species]:
                print('        >>>',strain)
                print('        STRAIN', strain, file=island_res_f_out)
                match_count = 1
                print('        match:', match_count)
                print('        MATCH', match_count, file=island_res_f_out)
                for match in blast_dt_splited['Y'][genus][species][strain]:
                    find_island = None
                    island_info = None
                    match = match.split()
                    query, subject, s_s, s_e, subject_expanded_seq = match[0], match[1], \
                        match[8], match[9], match[-1]
                    f_out_flag = open(flag_segment_file, 'w')
                    head = ' '.join(['>', os.path.basename(strain), str(match_count)])
                    print(head, file=f_out_flag)
                    print(subject_expanded_seq, file=f_out_flag)
                    f_out_flag.close()
                    query_file = flag_segment_file
                    for N_strain in species_N_strain:
                        ref_file = os.path.join(N_strain, os.path.basename(N_strain) + '.fasta')
                        numcer_res = run_nucmer(ref_file, query_file, nucmer_out_file)
                        if_have_island, island_info = parse_nucmer_res(numcer_res, ref_file, \
                            query_file, expand_len)
                        if if_have_island=='Y':
                            print('yaohho, find island in species blank strain...')
                            print('        BLANK_STARIN', ref_file, file=island_res_f_out)
                            print_island_info(query, subject, s_s, s_e, island_info,
                                subject_expanded_seq, island_res_f_out)
                            find_island = True
                            break
                    if not find_island:
                        print('contiune, looking in genus blank strain...')
                        for N_strain in genus_N_strain:
                            ref_file = os.path.join(N_strain, os.path.basename(N_strain) + '.fasta')
                            nucmer_res = run_nucmer(ref_file, query_file, nucmer_out_file)
                            if_have_island, island_info = parse_nucmer_res(nucmer_res, \
                                ref_file, query_file, expand_len)
                            if if_have_island=='Y':
                                print('yaohho, find island in genus blank strain...')
                                print('        BLANK_STARIN', ref_file, file=island_res_f_out)
                                print_island_info(query, subject, s_s, s_e, island_info,
                                    subject_expanded_seq, island_res_f_out)
                                find_island = True
                                break
                    match_count += 1
                    if find_island:
                        print(island_info)
    island_res_f_out.close()
    return 0


def main():
    (flag_seq_file, flag_seq_type, dataset, res_dir, blast_coverage_cutoff, \
        blast_identity_cutoff, expand_len)=get_args()
    f_out_file_name=blast_flag_gene(flag_seq_file,flag_seq_type,dataset,blast_coverage_cutoff,
        blast_identity_cutoff,res_dir, expand_len)
    callisland(f_out_file_name, res_dir, expand_len)
    return 0


if __name__=='__main__':
    main()





