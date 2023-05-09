
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
#from sourmash import MinHash
import itertools
import scipy
from scipy import spatial
import random


# for getting complementary DNA sequence
def complement(seq):
    dna_tab = str.maketrans("actg", "tgac")
    return seq.translate(dna_tab)

def filter_blasttable(input_file, output_directory="results", filter_params= None):
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
    with open(output_directory+"/"+"raw_blast_stats.txt", 'w') as f:
        f.writelines(raw_blast_stats)
    if filter_params is not None:
        blastout_filtered = blastout.query(filter_params)
        # sample query for blast purposes
        #"qseqid == 'IacA' and bitscore >= 250 or \\
        # qseqid == 'IacB' and bitscore >= 480 or \\
        # qseqid == 'IacC' and bitscore >= 90 or \\
        # qseqid == 'IacD' and bitscore >= 140 or \\
        # qseqid == 'IacE' and bitscore >= 240"
        blastout_filtered.to_csv(output_directory+"/"+'blastout_filtered.csv')
        #print raw blast out statistics:
        filtered_blast_stats=[]
        for i in prot_names:
            genome_count = len(list(set(blastout_filtered[blastout_filtered['qseqid']== i]['assembly'] )))
            total_count = len(blastout_filtered[blastout_filtered['qseqid']== i])
            filtered_blast_stats.append('Genomes with homologue to '+i+" "+str(genome_count)+"\n")
            filtered_blast_stats.append('Total homologue count for '+i+" "+str(total_count)+"\n")
        with open(output_directory+"/"+"filtered_blast_stats.txt", 'w') as f:
            f.writelines(filtered_blast_stats)
        return(blastout_filtered)
    else:
        blastout.to_csv(output_directory+"/"+'blastout.csv')
        return(blastout)

def parseblastout2(accession, blastout_dataframe, min_cluster_number=1, max_gene_dist= 10000, gene_color_dict= None):
    list_of_clusters=[]
    print("parsing proteins from: "+ accession )
    genome_match= blastout_dataframe[blastout_dataframe['accession'] ==accession]
    prot_names=set(blastout_dataframe["qseqid"])
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
        else:
            list_of_colors = [ "#"+str(hex(random.randint(0,16777215)))[2:] for c in range(0,len(prot_names))]
            gene_color_dict = dict(zip(prot_names,list_of_colors ))
            print(gene_color_dict)
        for cluster_number, df in genome_match_sorted.groupby('groups'):
            if len(df) >= min_cluster_number:
                hit_list = list(df['locus_tag'])
                old_locus_hit_list = list(df['old_locus_tag'])
                protein_name_list = list(df['protein_name'])
                protein_id_list = list(df['protein_id'])
                query_list = list(df['qseqid'])
                coord_list = list(zip(df['start_coord'], df['end_coord'],df['direction'],df['qseqid']))
                min_coord= min(df['start_coord'])
                max_coord= max(df['end_coord'])
                flip = sum(df['direction'])
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
                url_direction = "1" if flip < 0 else "0"
                start_url='https://www.ncbi.nlm.nih.gov/nuccore/'+ accession + '?report=graph&tracks=[key:sequence_track,name:Sequence,display_name:Sequence,id:STD649220238,annots:Sequence,ShowLabel:false,ColorGaps:false,shown:true,order:1][key:gene_model_track,name:Genes,display_name:Genes,annots:Unnamed,Options:MergeAll,CDSProductFeats:false,NtRuler:true,AaRuler:true,HighlightMode:2,ShowLabel:true,shown:true,order:6]&assm_context='+ assembly+'&mk='
                end_url= '&label=1&decor=1&spacing=2&v='+ str(min_coord)+':'+ str(max_coord)+'&gflip='+url_direction+'&c=000000&select=null&slim=0'
                ncbi_graphics_url = start_url + ncbi_graphics +end_url
                #sgi = re.sub("\{|\}|\'","", str(set(df['sgi'])) )
                cluster_len= max(df['end_coord']) - min(df['start_coord'])
                number_of_hits = len(df)
                pseudogene_list= list(df['pseudogene'])
                list_of_clusters.append([filename, biosample, number_of_hits , cluster_len, synteny, synteny_dir_dist, synteny_dir_pident, synteny_dir_evalue, synteny_dir_bitscore, synteny_dir_score, synteny_dir_length, synteny_dir, assembly, accession, name, hit_list, old_locus_hit_list, protein_name_list, protein_id_list, pseudogene_list, query_list, coord_list, cluster_number, ncbi_graphics_url])
        return(list_of_clusters)

def fetch_accessions(input_file):
    df=pd.read_csv(input_file)
    parse_blastp_input= list(set(df['accession']))
    return(parse_blastp_input)

def parseblastout2_parallel(list_of_accessions,blastout_file, n_cpus=1, max_gene_dist=10000, min_cluster_number=1, gene_color_dict=None):
    result_list=[]
    blastout_dataframe=pd.read_csv(blastout_file)
    with concurrent.futures.ProcessPoolExecutor(max_workers=n_cpus) as executor:
        for i in executor.map(parseblastout2, list_of_accessions, itertools.repeat(blastout_dataframe), itertools.repeat(min_cluster_number), itertools.repeat(max_gene_dist), itertools.repeat(gene_color_dict)):
            result_list.append(i)
            pass
    cluster_positive = list(filter(None, result_list))
    cluster_positive_flat = [item for sublist in cluster_positive for item in sublist]
    cluster_positive_df= pd.DataFrame(cluster_positive_flat, columns=('filename', 'biosample', 'number_of_hits' ,'cluster_len', 'synteny', 'synteny_dir_dist', 'synteny_dir_pident', 'synteny_dir_evalue', 'synteny_dir_bitscore', 'synteny_dir_score', 'synteny_dir_length', 'synteny_dir', 'assembly', 'accession', 'name', 'hit_list', 'old_locus_hit_list', 'protein_name_list', 'protein_id_list', 'pseudogene_list', 'query_list', 'coord_list', 'cluster_number', 'ncbi_graphics'))
    #filter so that x # or more hits must be clustered together
    cluster_positive_df=cluster_positive_df[[len(set(i))>=1 for i in cluster_positive_df['synteny'].str.findall('\w')]]
    cluster_positive_df['contig'] = cluster_positive_df['name'].str.contains('supercont|ctg|node|contig|scaffold|contigs',case=False)
    cluster_positive_df['complete_genome']= cluster_positive_df['name'].str.contains('complete',case=False)
    cluster_positive_df['plasmid'] = cluster_positive_df['name'].str.contains('plasmid')
    cluster_positive_df['has_overlap'] = cluster_positive_df['synteny_dir_dist'].str.contains("--")
    cluster_positive_df['duplicated']= cluster_positive_df.duplicated(subset="filename")
    cluster_positive_df.to_csv("results"+"/"+"cluster_hits.csv", index=False)
    return(cluster_positive_df)

def fetch_metadata(cluster_positive_file, ncbi_api_key):
    print('Fetching genome metadata')
    cluster_positive_df=pd.read_csv(cluster_positive_file)
    all_accs=list(set(cluster_positive_df['accession']))
    #fetch biosample id from accession
    # this works better than getting the metadata from the biosample name
    # there were errors
    batch=150
    accs_batches = [all_accs[i * batch:(i + 1) * batch] for i in range((len(all_accs) + batch - 1) // batch )]
    accs_batch_str=[]
    for i in accs_batches:
        id_list=re.sub('\n',',',str(i))
        id_list=re.sub(' ','',id_list)
        id_list=re.sub("NA\',|\'",'',id_list)
        id_list=re.sub("\[|\]","'",id_list)
        id_list=re.sub(",","&id=",id_list)
        id_list=re.sub("'","",id_list)
        accs_batch_str.append(id_list)
    gi_dict = {}
    biosample_ids_dict={}
    for i in accs_batch_str:
        accs_list = re.findall("(?<=\=)\w+|^\w+",i)
        url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?dbfrom=nuccore&db=biosample&id=" + i + "&api_key="+ncbi_api_key
        buffer = BytesIO()
        c = pycurl.Curl()
        c.setopt(c.URL, url)
        c.setopt(c.WRITEDATA, buffer)
        c.perform()
        c.close()
        body = buffer.getvalue()
        out=body.decode('iso-8859-1')
        root = ET.fromstring(out)
        for j in range(len(root.findall('./LinkSet'))):
            sgi= root.findall('./LinkSet')[j].findall('./IdList/Id')
            bsid= root.findall('./LinkSet')[j].findall('./LinkSetDb/Link/Id')
            accs_list[j]
            if bsid:
                biosample_ids_dict[str(accs_list[j])] = bsid[0].text
            else:
                biosample_ids_dict[str(accs_list[j])] = "NA"
            gi_dict[accs_list[j]] = sgi[0].text
    cluster_positive_df['biosample_id'] = cluster_positive_df['accession'].map(biosample_ids_dict)
    cluster_positive_df['gi'] = cluster_positive_df['accession'].map(gi_dict)
    all_ids=list(cluster_positive_df['biosample_id'])
    batch=150
    id_batches = [all_ids[i * batch:(i + 1) * batch] for i in range((len(all_ids) + batch - 1) // batch )]
    id_batch_str=[]
    for i in id_batches:
        id_list=re.sub('\n',',',str(i))
        id_list=re.sub(' ','',id_list)
        id_list=re.sub("NA\',|\'",'',id_list)
        id_list=re.sub("\[|\]","'",id_list)
        id_batch_str.append(id_list)
    isolation_source_dict={}
    isolation_source_dict["NA"] = "NA"
    env_biome_dict={}
    env_biome_dict["NA"] = "NA"
    env_feature_dict={}
    env_feature_dict["NA"] = "NA"
    host_dict={}
    host_dict["NA"] = "NA"
    host_sciname_dict={}
    host_sciname_dict["NA"] = "NA"
    strain_dict={}
    strain_dict["NA"] = "NA"
    all_meta_dict={}
    all_meta_dict["NA"]="NA"
    for i in id_batch_str:
        url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=biosample&id=" + i + "&api_key="+ncbi_api_key
        print(url)
        buffer = BytesIO()
        c = pycurl.Curl()
        c.setopt(c.URL, url)
        c.setopt(c.WRITEDATA, buffer)
        c.perform()
        c.close()
        body = buffer.getvalue()
        out=body.decode('iso-8859-1')
        root = ET.fromstring(out)
        for i in root.findall('./BioSample'):
            bs_text = ET.tostring(i).decode('ISO-8859-1')
            bs_id=re.findall('(?<=\sid\=")\w+',bs_text)[0]
            isolation_source= i.findall('./Attributes/Attribute[@attribute_name="isolation_source"]')
            if isolation_source:
                isolation_source_dict[bs_id] = isolation_source[0].text
            else:
                isolation_source_dict[bs_id] = "NA"
            env_biome= i.findall('./Attributes/Attribute[@attribute_name="env_biome"]')
            if env_biome:
                env_biome_dict[bs_id] = env_biome[0].text
            else:
                env_biome_dict[bs_id] = "NA"
            env_feature= i.findall('./Attributes/Attribute[@attribute_name="env_feature"]')
            if env_feature:
                env_feature_dict[bs_id] = env_feature[0].text
            else:
                env_feature_dict[bs_id] = "NA"
            host= i.findall('./Attributes/Attribute[@attribute_name="host"]')
            if host:
                host_dict[bs_id] = host[0].text
            else:
                host_dict[bs_id] = "NA"
            host_sciname= i.findall('./Attributes/Attribute[@attribute_name="host scientific name"]')
            if host_sciname:
                host_sciname_dict[bs_id] = host_sciname[0].text
            else:
                host_sciname_dict[bs_id] = "NA"
            strain= i.findall('./Attributes/Attribute[@attribute_name="strain"]')
            if strain:
                strain_dict[bs_id] = strain[0].text
            else:
                strain_dict[bs_id] = "NA"
            all_metadata = i.findall('./Attributes/Attribute[@attribute_name]')
            if all_metadata:
                all_metadata_list=[]
                for m in all_metadata:
                    all_metadata_list.append( (m.attrib['attribute_name'], m.text ))
                all_meta_dict[bs_id]= all_metadata_list
            else:
                all_meta_dict[bs_id] = "NA"
    #print(isolation_source_dict)
    cluster_positive_df['isolation_src'] = cluster_positive_df['biosample_id'].map(isolation_source_dict)
    cluster_positive_df['env_biome'] = cluster_positive_df['biosample_id'].map(env_biome_dict)
    cluster_positive_df['env_feature'] = cluster_positive_df['biosample_id'].map(env_feature_dict)
    cluster_positive_df['host'] = cluster_positive_df['biosample_id'].map(host_dict)
    cluster_positive_df['host_sciname'] = cluster_positive_df['biosample_id'].map(host_sciname_dict)
    cluster_positive_df['strain'] = cluster_positive_df['biosample_id'].map(strain_dict)
    cluster_positive_df['all_biosample_metadata'] = cluster_positive_df['biosample_id'].map(all_meta_dict)
    #match to gtdb independant;y of assembly version
    if not os.path.exists("gtdb/"):
        os.mkdir("gtdb/")
    if not os.path.exists("gtdb/bac120_metadata_r89.tsv"):
        print("downloading gtdb data")
        urllib.request.urlretrieve("https://data.ace.uq.edu.au/public/gtdb/data/releases/release89/89.0/bac120_metadata_r89.tsv",'gtdb/bac120_metadata_r89.tsv')
    cluster_positive_df['assembly_base']= cluster_positive_df['assembly'].str.split(pat='\.',expand=True)[0]
    all_assemblies= ['RS_'+ i for i in list(set(cluster_positive_df['assembly_base']))]
    gtdb_metadata = pd.read_csv('gtdb/bac120_metadata_r89.tsv', sep='\t')
    gtdb_metadata['accession_base'] = gtdb_metadata['accession'].str.split(pat='\.',expand=True)[0]
    gtdb_matches= gtdb_metadata[gtdb_metadata['accession_base'].isin(all_assemblies)]
    gtdb_matches['accession_base_'] = gtdb_matches['accession_base'].str.replace('RS_','')
    #gtdb_dict = gtdb_matches[['accession_base_','gtdb_taxonomy']].to_dict()
    gtdb_dict = dict(zip(gtdb_matches['accession_base_'], gtdb_matches['gtdb_taxonomy']))
    ncbi_dict = dict(zip(gtdb_matches['accession_base_'], gtdb_matches['ncbi_taxonomy']))
    cluster_positive_df['gtdb_tax'] = cluster_positive_df['assembly_base'].map(gtdb_dict)
    cluster_positive_df['ncbi_tax'] = cluster_positive_df['assembly_base'].map(ncbi_dict)
    cluster_positive_df['same_taxonomy'] = cluster_positive_df['gtdb_tax'] == cluster_positive_df['ncbi_tax']
    cluster_positive_df[['domain_gtdb','phylum_gtdb','class_gtdb','order_gtdb','family_gtdb','genus_gtdb','species_gtdb']]=cluster_positive_df['gtdb_tax'].str.split(";", expand = True)
    cluster_positive_df.to_csv(cluster_positive_file, index=False)
    return(cluster_positive_df)

def make_indexprot(file,fasta_file_directory):
    #header_re = re.compile('^>.+', re.MULTILINE)
    prot_re = re.compile('^>[\s|\S]+?(?=\n\s*\n)', re.MULTILINE)
    if not os.path.exists("index_files"):
        os.mkdir("index_files")
    print("attempting indexing for " + file)
    if os.path.exists(fasta_file_directory+'/'+file):
        filename = fasta_file_directory+'/'+file
        #outfilename = 'index_files/'+ file + ".index"
        outfilename = 'index_files/'+ file + ".indexprot"
        if os.path.exists(outfilename):
            print(outfilename + " already indexed")
        else:
            with open(filename, 'r') as file:
                file_contents = file.read()
                proteins = prot_re.findall(file_contents)
                out_prots = []
                for p in proteins:
                    p1 = re.sub('NULL$','NULL!!',p)
                    p2 = re.sub('PSEUDOGENE$','PSEUDOGENE!!',p1)
                    p3 = p2.replace('\n','!!',1)
                    p4 = p3.replace('\n','??')
                    out_prots.append(p4)
            with open(outfilename, 'w+') as outfile:
                outfile.write('\n'.join(out_prots) + '\n')
            print('made index for: ' + filename)
    else:
        print('No match for' + file)

def make_indexprot_parallel(cluster_positive_file, fasta_file_directory="fasta_files"):
    cluster_positive_df=pd.read_csv(cluster_positive_file)
    inputs_indexprot = [re.sub('.gbff','_proteins.fa', i) for i in list(set(cluster_positive_df['filename']))]
    if not os.path.exists("index_files"):
        os.mkdir("index_files")
    with concurrent.futures.ThreadPoolExecutor(max_workers=None) as executor:
        for i in executor.map(make_indexprot, inputs_indexprot, itertools.repeat(fasta_file_directory)):
            pass


def fetch_neighborhood2(index,cluster_positive_file,features_upstream = 0,features_downstream = 0):
    cluster_positive_df=pd.read_csv(cluster_positive_file)
    cluster = cluster_positive_df.iloc[index,:]
    acc= cluster['accession']
    assembly = re.sub('.gbff','_proteins.fa.indexprot', cluster['filename'])
    #make the genome database from the .fa.index file
    assembly_index_file ='index_files/'+assembly
    print(assembly_index_file)
    db = pd.read_csv(assembly_index_file, sep = "!!" ,header = None, engine='python')
    #db.columns = ["filename","assembly","accession","locus_tag","old_locus_tag","name","biosample","protein_name","coordinates","protein_id"]
    db.columns = ["filename","assembly","accession","locus_tag","old_locus_tag","name","biosample","protein_name","coordinates","protein_id","pseudogene","protein_seq"]
    db['direction']= [-1 if re.match('complement',c) else 1 for c in db['coordinates']]
    db['start_coord'] = [re.search('\d+?(?=\.\.(\d|\>))',str(c)).group(0) for c in db['coordinates'] ]
    db['start_coord'] = [re.sub('complement|>|<|\)|\(',"",c) for c in db['start_coord'] ]
    db['start_coord'] = db['start_coord'].astype(int)
    db['end_coord'] = [re.search('(?<=\.(\.|\>))\d+',str(c)).group(0) for c in db['coordinates'] ]
    db['end_coord'] = [re.sub('>|<|\)|\(',"",c) for c in db['end_coord'] ]
    db['end_coord'] = db['end_coord'].astype(int)
    hit_list = cluster['hit_list']
    query_list = cluster['query_list']
    cluster_number = cluster['cluster_number']
    hit_dict = dict(zip(hit_list,query_list))
    genome = db.loc[db['accession'] == acc].copy()
    start = genome[genome['locus_tag'] == hit_list[0]].index.values.astype(int)[0] - features_upstream
    stop = genome[genome['locus_tag'] == hit_list[-1]].index.values.astype(int)[0] + features_downstream
    neighborhood = genome.loc[start:stop,].copy()
    neighborhood['query_match'] = neighborhood['locus_tag'].map(hit_dict)
    coord_list = list(zip(neighborhood['start_coord'], neighborhood['end_coord'],neighborhood['direction'],neighborhood['query_match']))
    #function to find GC content of cluster vs genome
    gbff_str=str(db['filename'][0][1:])
    with open("gbff_files_unzipped/"+gbff_str) as file:
        gbff_file = file.read()
    genome_seq = "".join(re.findall("(?<=ORIGIN)[\s+\S+]+?(?=\/\/)",gbff_file))
    genome_seq = re.sub('\s|\d|\n','',genome_seq)
    Gg=genome_seq.count("g")
    Gc=genome_seq.count("c")
    Ga=genome_seq.count("a")
    Gt=genome_seq.count("t")
    genomeGC= (Gg+Gc)/(Gg+Gc+Ga+Gt)
    start=min(coord_list)[0]
    end=max(coord_list)[1]
    regex_str=acc+"[\s+\S+]+\/\/"
    all_cluster_fasta = re.findall(regex_str,gbff_file)[0]
    all_cluster_fasta = re.findall("(?<=ORIGIN)[\s+\S+]+(?=\/\/)",all_cluster_fasta)[0]
    all_cluster_fasta = re.sub(" |\d|\n","",all_cluster_fasta)
    cluster_seq = all_cluster_fasta[start-1:end-1]
    g=cluster_seq.count("g")
    c=cluster_seq.count("c")
    a=cluster_seq.count("a")
    t=cluster_seq.count("t")
    clusterGC = (g+c)/(g+c+a+t)
    diffGC = abs(clusterGC - genomeGC)
    #compare minhash values between cluster and genome
    #kmer_size=5
    #n=0
    #sc=1
    #cluster_minhash= MinHash(n=n, ksize=kmer_size,scaled=sc)
    #cluster_minhash.add_sequence(cluster_seq,force=True)
    #cluster_minhash.add_sequence(complement(cluster_seq),force=True)
    #
    #genome_minhash= MinHash(n=n, ksize=kmer_size,scaled=sc)
    #genome_minhash.add_sequence(genome_seq,force=True)
    #genome_minhash.add_sequence(complement(genome_seq),force=True)
    #minhash_sim=cluster_minhash.similarity(genome_minhash)
    # genome_minus_cluster=re.sub(cluster_seq,'',genome_seq)
    # #print(len(genome_seq)-len(genome_minus_cluster))
    # genome_minus_cluster_minhash=MinHash(n=n, ksize=kmer_size,scaled=sc)
    # genome_minus_cluster_minhash.add_sequence(genome_minus_cluster,force=True)
    # genome_minus_cluster_minhash.add_sequence(complement(genome_minus_cluster),force=True)
    # minhash_sim_minus_cluster=cluster_minhash.similarity(genome_minus_cluster_minhash)
    #print(minhash_sim)
    #compare tetranucleotide frequency between cluster and genomes
    bases=['a','t','g','c']
    four_mers=[''.join(p) for p in itertools.product(bases, repeat=4)]
    four_mer_count_genome= np.add([genome_seq.count(i) for i in four_mers], [complement(genome_seq).count(i) for i in four_mers])
    four_mer_freq_genome = [i/sum(four_mer_count_genome) for i in four_mer_count_genome]
    four_mer_count_cluster=np.add([cluster_seq.count(i) for i in four_mers], [complement(cluster_seq).count(i) for i in four_mers])
    four_mer_freq_cluster = [i/sum(four_mer_count_cluster) for i in four_mer_count_cluster]
    four_mer_distance=scipy.spatial.distance.cityblock(four_mer_freq_cluster,four_mer_freq_genome)
    ####
    if sum( neighborhood[neighborhood['query_match'].notnull()]['direction'] ) < 0:
            neighborhood['actual_start_tmp'] = neighborhood['start_coord']
            neighborhood['start_coord']= neighborhood['end_coord'] * -1
            neighborhood['end_coord']= neighborhood['actual_start_tmp']* -1
            neighborhood['direction'] = neighborhood['direction'] *-1
            neighborhood = neighborhood.sort_values(by='start_coord')
    neighborhood['query_match'] = neighborhood['query_match'].replace(np.nan,"x")
    nhbrhood_hit_list= list(neighborhood['query_match'])
    nhbrhood_locus_tags= list(neighborhood['locus_tag'])
    nhbrhood_old_locus_tags= list(neighborhood['old_locus_tag'])
    nhbrhood_prot_ids= list(neighborhood['protein_id'])
    nhbrhood_prot_name= list(neighborhood['protein_name'])
    nhbrhood_prot_seq= list(neighborhood['protein_seq'])
    order = [ ("| " + gene['query_match'] + " 〉") if gene['direction'] == 1 else ("〈 " + gene['query_match'] + " |") for index,gene in neighborhood.iterrows()  ]
    dist = list(np.array(neighborhood['start_coord'][1:]) -  np.array(neighborhood['end_coord'][:-1]))
    dist = ["-"+str(d)+"-" for d in dist]
    adj_coord_list = list(zip(neighborhood['start_coord'], neighborhood['end_coord'],neighborhood['direction'],neighborhood['query_match']))
    if min(neighborhood['start_coord']) <0:
        tare_value = abs(min(neighborhood['start_coord']))
        tared_adj_coord_list = list(zip([v + tare_value for v in neighborhood['start_coord']], [v + tare_value for v in neighborhood['end_coord']],neighborhood['direction'],neighborhood['query_match']))
    else:
        tare_value = min(neighborhood['start_coord'])
        tared_adj_coord_list = list(zip([v - tare_value for v in neighborhood['start_coord']], [v - tare_value for v in neighborhood['end_coord']],neighborhood['direction'],neighborhood['query_match']))
    # making an ITOL compatible string
    #gene_color_dict={ 'IacA':'#ff5969',
                        'IacB':'#2db34e',
                        'IacC':'#fb77e0',
                        'IacD':'#00bc7e',
                        'IacE':'#8d006e',
                        'IacF':'#cfdd63',
                        'IacG':'#0060d0',
                        'IacR':'#bb7b00',
                        'IacH':'#7c2c29',
                        'IacI':'#f1d17a',
                        'x':'#d1d1d1'}
    #max_len = tared_adj_coord_list[-1][1]
    #itol_diagram=[]
    #for g in tared_adj_coord_list:
    #    gene_string=[]
    #    gene_length=g[1]-g[0]
    #    if g[2] > 0:
    #        gene_string.append('RE')
    #        gene_string.append(str(g[0]))
    #        gene_string.append(str(g[1]-(0.1*gene_length)))
    #        #gene_string.append('#34b4eb')
    #        gene_string.append(gene_color_dict[g[3]])
    #        gene_string.append(str(g[3]))
    #        gene_string.append(',')
    #        gene_string.append('TR')
    #        gene_string.append(str(g[1]-(0.1*gene_length)))
    #        gene_string.append(str(g[1]))
    #        #gene_string.append('#34b4eb')
    #        gene_string.append(gene_color_dict[g[3]])
    #        gene_string.append('')
    #    else:
    #        gene_string.append('TL')
    #        gene_string.append(str(g[0]))
    #        gene_string.append(str(g[0]+(0.1*gene_length)))
    #        #gene_string.append('#34b4eb')
    #        gene_string.append(gene_color_dict[g[3]])
    #        gene_string.append('')
    #        gene_string.append(',')
    #        gene_string.append('RE')
    #        gene_string.append(str(g[0]+(0.1*gene_length)))
    #        gene_string.append(str(g[1]))
            #gene_string.append('#34b4eb')
    #        gene_string.append(gene_color_dict[g[3]])
    #        gene_string.append(str(g[3]))
    #    itol_gene='|'.join(gene_string)
    #    itol_diagram.append(itol_gene)
    #itol_diagram_joined=",".join(map(str,itol_diagram))
    #itol_diagram_string= str(max_len)+','+itol_diagram_joined
    #itol_diagram_string = re.sub(',\|',',',itol_diagram_string)
    #obtains "| A 〉-23-| B 〉-23-| C 〉"
    synteny_dir_dist=''.join(sum(zip(order, dist+[0]), ())[:-1])
    synteny_dir_dist = re.sub("Iac" ,"", synteny_dir_dist)
    #obtains "| A 〉| B 〉| C 〉"
    synteny_dir =''.join(order)
    synteny_dir = re.sub("Iac" ,"", synteny_dir)
    #obtains "| A:23.23 〉| B:23.23〉| C:23.23 〉"
    #synteny_dir_pident =''.join(order_pident)
    #synteny_dir_pident = re.sub("Iac" ,"", synteny_dir_pident)
    #obtains "A-B-C"
    synteny= re.sub("\n" ,"-", neighborhood['query_match'].to_string(index=False))
    synteny= re.sub("Iac| " ,"", synteny)
    synteny_alphabet = "".join([ gene['query_match'].replace("Iac","").upper() if gene['direction'] == 1 else gene['query_match'].replace("Iac","").lower() for index,gene in neighborhood.iterrows()  ])
    cluster_len= max(neighborhood['end_coord']) - min(neighborhood ['start_coord'])
    assembly= re.sub("\{|\}|\'|>","", str(set(neighborhood['assembly'])) )
    accession = re.sub("\{|\}|\'","", str(set(neighborhood['accession'])) )
    title= re.sub("\{|\}|\'", "",str(set(neighborhood['name'])) )
    print(assembly_index_file + " successfully used")
    #return([accession, assembly, title, len(neighborhood), cluster_len, synteny,synteny_alphabet, synteny_dir_dist, synteny_dir, #cluster_number,coord_list,adj_coord_list,tared_adj_coord_list,itol_diagram_string, #nhbrhood_hit_list,nhbrhood_locus_tags,nhbrhood_old_locus_tags,nhbrhood_prot_ids,nhbrhood_prot_name,nhbrhood_prot_seq, clusterGC, #genomeGC,diffGC,minhash_sim, four_mer_freq_cluster,four_mer_freq_genome, four_mer_distance, cluster_seq])
    return([accession, assembly, title, len(neighborhood), cluster_len, synteny,synteny_alphabet, synteny_dir_dist, synteny_dir, cluster_number,coord_list,adj_coord_list,tared_adj_coord_list, nhbrhood_hit_list,nhbrhood_locus_tags,nhbrhood_old_locus_tags,nhbrhood_prot_ids,nhbrhood_prot_name,nhbrhood_prot_seq, clusterGC, genomeGC,diffGC, four_mer_freq_cluster,four_mer_freq_genome, four_mer_distance, cluster_seq])


def helper_fetchneighborhood2(index):
    return fetchneighborhood2(index,features_upstream = 0,features_downstream = 0)
#
def fetch_neighborhood_parallel():
    fetchneighborhood_columns=['accession','assembly', 'title', 'feature_count_nhbr', 'cluster_len_nhbr', 'synteny_nhbr','synteny_alphabet_nhbr', 'synteny_dir_dist_nhbr', 'synteny_dir_nhbr','cluster_number','coord_list','adj_coord_list','tared_adj_coord_list', 'nhbrhood_hit_list','nhbrhood_locus_tags','nhbrhood_old_locus_tags','nhbrhood_prot_ids','nhbrhood_prot_name','nhbrhood_prot_seq', 'clusterGC','genomeGC','diffGC','four_mer_freq_cluster','four_mer_freq_genome','four_mer_distance','cluster_seq']
    inputs_fetchneighborhood = list(range(0,len(iac_positive_df)))
    outputs_fetchneighborhood=[]
    with concurrent.futures.ProcessPoolExecutor(max_workers=None) as executor:
        for i in executor.map(fetchneighborhood2, inputs_fetchneighborhood):
            outputs_fetchneighborhood.append(i)
            pass
