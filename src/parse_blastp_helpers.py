
import re
import sys
import pandas as pd
import pycurl
from io import BytesIO
import time
import xml.etree.ElementTree as ET
import numpy as np
import concurrent.futures
import urllib.request
import glob
import os
from sourmash import MinHash
import itertools
import scipy
from scipy import spatial
import random


# for getting complementary DNA sequence
def complement(seq):
    dna_tab = str.maketrans("actg", "tgac")
    return seq.translate(dna_tab)

def filter_blasttable(input_file, filter_params= None):
    blastout = pd.read_csv(input_file, sep='\t', low_memory = False)
    print('Editing blastp output file')
    #blastout[["filename","assembly","accession","locus_tag","old_locus_tag","name","biosample","protein_name","coordinates","protein_id","pseudogene"]] = blastout['stitle'].str.split("!!", expand = True)
    blastout[["filename","assembly","accession","locus_tag","old_locus_tag","name","biosample","protein_name","coordinates","protein_id","pseudogene"]] = blastout['sacc'].str.split("!!", expand = True)
    blastout['direction']= [-1 if re.match('complement',str(c)) else 1 for c in blastout['coordinates']]
    blastout['join_status'] = [1 if re.match('join',str(c)) else 0 for c in blastout['coordinates']]
    #extracts start and end coordinates from the coordinates column
    blastout['start_coord'] = [re.search('\d+?(?=\.\.(\d|\>))',str(c)).group(0) for c in blastout['coordinates'] ]
    blastout['start_coord'] = [re.sub('join|complement|join\(|complement|>|<|\)|\(|,\d+',"",str(c)) for c in blastout['start_coord'] ]
    blastout['start_coord'] = blastout['start_coord'].astype(int)
    blastout['end_coord'] = [re.search('(?<=\.(\.|\>))\d+',str(c)).group(0) for c in blastout['coordinates'] ]
    blastout['end_coord'] = [re.sub('>|<|\)|\(',"",c) for c in blastout['end_coord'] ]
    blastout['end_coord'] = blastout['end_coord'].astype(int)
    blastout['middle_coord'] = (  blastout['start_coord'] + blastout['end_coord'] )/2
    print('Filtering blastp output file')
    #new code 12.2.21
    #blastout_filtered = blastout[(blastout["qseqid"].isin(["IacA","IacB","IacC","IacD","IacE","IacH","IacI"])) &
    #        (blastout['evalue'] <= 1e-10) & (blastout['bitscore']>=50)]
    prot_names=set(blastout['qseqid'])
    raw_blast_stats=[]
    for i in prot_names:
        genome_count = len(list(set(blastout[blastout['qseqid']== i]['assembly'] )))
        total_count = len(blastout[blastout['qseqid']== i])
        raw_blast_stats.append('Genomes with homologue to '+i+" "+str(genome_count)+"\n")
        raw_blast_stats.append('Total homologue count for '+i+" "+str(total_count)+"\n")
    with open("raw_blast_stats.txt", 'w') as f:
        f.writelines(raw_blast_stats)
    if filter_params is not None:
        blastout_filtered = df.query(filter_params)
        # sample query for blast purposes
        #"qseqid == 'IacA' and bitscore >= 250 or \\
        # qseqid == 'IacB' and bitscore >= 480 or \\
        # qseqid == 'IacC' and bitscore >= 90 or \\
        # qseqid == 'IacD' and bitscore >= 140 or \\
        # qseqid == 'IacE' and bitscore >= 240"
        blastout_filtered.to_csv('blastout_filtered.csv')
        #print raw blast out statistics:
        filtered_blast_stats=[]
        for i in prot_names:
            genome_count = len(list(set(blastout_filtered[blastout_filtered['qseqid']== i]['assembly'] )))
            total_count = len(blastout_filtered[blastout_filtered['qseqid']== i])
            filtered_blast_stats.append('Genomes with homologue to '+i+" "+str(genome_count)+"\n")
            filtered_blast_stats.append('Total homologue count for '+i+" "+str(total_count)+"\n")
        with open("filtered_blast_stats.txt", 'w') as f:
            f.writelines(filtered_blast_stats)
        return(blastout_filtered)
    else:
        return(blastout)

def parseblastout2(accession, blastout_dataframe, min_cluster_number=1, max_gene_dist= 10000, gene_color_dict= None):
    print("parsing proteins from: "+ accession )
    genome_match= blastout_filtered[blastout_filtered['accession'] ==accession]
    prot_names=set(blastout_filtered["qseqid"])
    #if len(genome_match) > 4:
    #the following filters the genomes if they have 6 of 7 of the iacABCDEH or I
    #this filter takes into account having multiples of a hit and filters based on presence or absence
    if sum(gene in list(genome_match['qseqid']) for gene in prot_names) >= min_cluster_number:
        print("found matches in"+accession)
        genome_match_sorted = genome_match.sort_values(by='start_coord')
        middle_coords = genome_match_sorted['middle_coord']
        middle_coords_np = np.array(middle_coords)
        genome_match_sorted['groups']=list(np.cumsum([0] + list(1*(middle_coords_np[1:] - middle_coords_np[0:-1] >= max_gene_dist))) + 1)
        if gene_color_dict is not None:
            gene_color_dict= gene_color_dict
        else
            list_of_clusters = [hex(random.randrange(0, 2**24)) for c in range(0,len(prot_names))]
            gene_color_dict = dict(zip(prot, col) for prot, col in prot_names, )
        for cluster_number, df in genome_match_sorted.groupby('groups'):
            if len(df) >= 1:
                hit_list = list(df['locus_tag'])
                old_locus_hit_list = list(df['old_locus_tag'])
                protein_name_list = list(df['protein_name'])
                protein_id_list = list(df['protein_id'])
                query_list = list(df['qseqid'])
                coord_list = list(zip(df['start_coord'], df['end_coord'],df['direction'],df['qseqid']))
                ncbi_graphics=''
                for i in coord_list:
                    ncbi_graphics= ncbi_graphics+ str(i[0])+':'+str(i[1])+"|"+i[3]+"|"+gene_color_dict[i[3]]+","
                ncbi_graphics= ncbi_graphics[0:-1]
                if sum(df['direction']) < 0:
                    df['actual_start_tmp'] = df['start_coord']
                    df['start_coord']= df['end_coord'] * -1
                    df['end_coord']= df['actual_start_tmp']* -1
                    df['direction'] = df['direction'] *-1
                    df.sort_values(by='start_coord',inplace=True)
                order = [ ("| " + gene['qseqid'] + " 〉") if gene['direction'] == 1 else ("〈 " + gene['qseqid'] + " |") for index,gene in df.iterrows()  ]
                order_pident = [ ("| " + gene['qseqid']+':'+ str(round(gene['pident'],2)) + " 〉") if gene['direction'] == 1 else ("〈 " +  gene['qseqid']+':'+ str(round(gene['pident'],2)) + " |") for index,gene in df.iterrows()  ]
                order_evalue = [ ("| " + gene['qseqid']+':'+ str(gene['evalue']) + " 〉") if gene['direction'] == 1 else ("〈 " +  gene['qseqid']+':'+ str(gene['evalue']) + " |") for index,gene in df.iterrows()  ]
                order_bitscore = [ ("| " + gene['qseqid']+':'+ str(gene['bitscore']) + " 〉") if gene['direction'] == 1 else ("〈 " +  gene['qseqid']+':'+ str(gene['bitscore']) + " |") for index,gene in df.iterrows()  ]
                order_score = [ ("| " + gene['qseqid']+':'+ str(gene['score']) + " 〉") if gene['direction'] == 1 else ("〈 " +  gene['qseqid']+':'+ str(gene['score']) + " |") for index,gene in df.iterrows()  ]
                order_length = [ ("| " + gene['qseqid']+':'+ str(gene['length']) + " 〉") if gene['direction'] == 1 else ("〈 " +  gene['qseqid']+':'+ str(gene['length']) + " |") for index,gene in df.iterrows()  ]
                dist = list(np.array(df['start_coord'][1:]) -  np.array(df['end_coord'][:-1]))
                if len(dist) > 1:
                    dist_list = ["-"+str(d)+"-" for d in dist]
                else:
                    dist_list = str(0)
                #obtains "| A 〉-23-| B 〉-23-| C 〉"
                if len(dist) > 1:
                    synteny_dir_dist=''.join(sum(zip(order, dist_list+[0]), ())[:-1])
                else:
                    synteny_dir_dist = order
                #obtains "| A 〉| B 〉| C 〉"
                synteny_dir =''.join(order)
                #obtains "| A:23.23 〉| B:23.23〉| C:23.23 〉"
                synteny_dir_pident =''.join(order_pident)
                #obtains "| A:bitscore 〉| B:bitscore〉| C:bitscore 〉"
                synteny_dir_bitscore =''.join(order_bitscore)
                #obtains "| A:score 〉| B:score〉| C:score 〉"
                synteny_dir_score =''.join(order_score)
                #obtains "| A:length 〉| B:length〉| C:length 〉"
                synteny_dir_length =''.join(order_length)
                #obtains "| A:evalue 〉| B:evalue〉| C:evalue 〉"
                synteny_dir_evalue =''.join(order_evalue)
                #obtains "A-B-C"
                synteny= re.sub("\n" ,"-", df['qseqid'].to_string(index=False))
                filename = list(set(df['filename']))[0]
                biosample = list(set(df['biosample']))[0]
                name = list(set(df['name']))[0]
                assembly= re.sub("\{|\}|\'","", str(set(df['assembly'])) )
                accession = re.sub("\{|\}|\'","", str(set(df['accession'])) )
                #sgi = re.sub("\{|\}|\'","", str(set(df['sgi'])) )
                cluster_len= max(df['end_coord']) - min(df['start_coord'])
                number_of_hits = len(df)
                pseudogene_list= list(df['pseudogene'])
                list_of_clusters.append([filename, biosample, number_of_hits , cluster_len, synteny, synteny_dir_dist, synteny_dir_pident, synteny_dir_evalue, synteny_dir_bitscore, synteny_dir_score, synteny_dir_length, synteny_dir, assembly, accession, name, hit_list, old_locus_hit_list, protein_name_list, protein_id_list, pseudogene_list, query_list, coord_list, cluster_number, ncbi_graphics])
        return(list_of_clusters)

def fetch_accessions(input_dataframe):
    parse_blastp_input= list(set(blastout_filtered['accession']))
    return()


def parseblastout2_parallel(list_of_accessions,blastout_dataframe, n_cpus=1, max_gene_dist=10000, min_cluster_number=1, gene_color_dict=None):
    result_list=[]
    with concurrent.futures.ProcessPoolExecutor(max_workers=n_cpus) as executor:
        for i in executor.map(parseblastout2, list_of_accessions, itertools.repeat(blastout_dataframe), itertools.repeat(min_cluster_number), itertools.repeat(max_gene_dist), itertools.repeat(gene_color_dict):
            result_list.append(i)
            pass
    cluster_positive = list(filter(None, result_list))
    cluster_positive_flat = [item for sublist in cluster_positive for item in sublist]
    cluster_positive_df= pd.DataFrame(cluster_positive_flat, columns=('filename', 'biosample', 'number_of_hits' ,'cluster_len', 'synteny', 'synteny_dir_dist', 'synteny_dir_pident', 'synteny_dir_evalue', 'synteny_dir_bitscore', 'synteny_dir_score', 'synteny_dir_length', 'synteny_dir', 'assembly', 'accession', 'name', 'hit_list', 'old_locus_hit_list', 'protein_name_list', 'protein_id_list', 'pseudogene_list', 'query_list', 'coord_list', 'cluster_number', 'ncbi_graphics'))
    #filter so that x # or more hits must be clustered together
    cluster_positive_df=cluster_positive_df[[len(set(i))>=1 for i in cluster_positive_df['synteny'].str.findall('\w')]]
    #
    cluster_positive_df['contig'] = cluster_positive_df['name'].str.contains('supercont|ctg|node|contig|scaffold|contigs',case=False)
    cluster_positive_df['complete_genome']= cluster_positive_df['name'].str.contains('complete',case=False)
    cluster_positive_df['plasmid'] = cluster_positive_df['name'].str.contains('plasmid')
    cluster_positive_df['has_overlap'] = cluster_positive_df['synteny_dir_dist'].str.contains("--")
    cluster_positive_df['duplicated']= cluster_positive_df.duplicated(subset="filename")
    return(cluster_positive_df)
