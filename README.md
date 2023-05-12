

# sample workflow:

The following is a sample workflow that searches for the well characterized phenylacetic acid gene cluster within **Pseudomonas putida** genomes in the RefSeq database.

### Downloading gbff files
The first step is to download files from RefSeq.
This can be done manually via a web interface if desired however the ```download_gbff_files_parallel()``` function enables doing this from the command line. The function first downloads the RefSeq assembly summary file and then gets the ftp paths for each genome for subsequent download. If you are running the workflow on an HPC or with lots of storage you can certainly download all RefSeq genomes. However, the function takes a parameter ```filter_params``` which lets you subset the entire RefSeq dataset to only include geomes of interest (for example just from a particular genus or species). For the example here the filter will allow only complete genomes which contain the word "Pseudomonas putida" and that are the latest version of the genome entry. You can filter based on any column of the assembly summary file using the  

```
download_gbff_files_parallel(filter_params="assembly_level == 'Complete Genome' & version_status == 'latest' & organism_name.str.contains('Pseudomonas putida')")
```

The files will be downloadeed to a folder named ```gbff_files```

The files will then be unzipped to a folder named ```gbff_unzipped_files```

### Extracting protein information from gbff files

After downloading and extracting the gbff files the next step is to parse the files to construct fasta protein files for a subsequent blast search.
Each protein within the gbff file gets extracted along with relevant information which is stored in the fasta header.
A sample extracted protein is listed below. Each component of the header is separated by a delimeter of two exclamation points ```!!``` Interestingly it was only one exclamation point but the parser got confused because there was a genome name that contained an exclamation point. Some people give their genome entries strange names.
The header entries are ```>filename!!Assembly_id!!Nucleotide_accession!!Genome_name!!Locus_tag!!Old_locus_tag!!Biosample_id!!Protein_name!!DNA_coordinates!!Protein_id!!Pseudogene```

The last entry in the header will be NULL unless the gene was classified as a pseudogene by the RefSeq pipeline.
```
>GCF_000016865.1_ASM1686v1_genomic.gbff!!GCF_000016865.1!!NC_009512!!PPUT_RS00005!!Pput_0001!!Pseudomonas_putida_F1__complete_genome.!!SAMN00623058!!chromosomal_replication_initiator_protein_DnaA!!385..1902!!WP_003253158.1!!NULL
MSVELWQQCVELLRDELPAQQFNTWIRPLQVEAEGDELRVYAPNRFVLDWVNEKYLGRLLELLGENGSGI
APALSLLIGSRRSSAPRAAPNAPVSAAVAASLAQTQAHKTAPAAAVEPVAVAAAEPVLVETSSRDSFDAM
AEPAAAPPSGGRAEQRTVQVEGALKHTSYLNRTFTFDTFVEGKSNQLARAAAWQVADNPKHGYNPLFLYG
GVGLGKTHLMHAVGNHLLKKNPNAKVVYLHSERFVADMVKALQLNAINEFKRFYRSVDALLIDDIQFFAR
KERSQEEFFHTFNALLEGGQQVILTSDRYPKEIEGLEERLKSRFGWGLTVAVEPPELETRVAILMKKADQ
AKVELPHDAAFFIAQRIRSNVRELEGALKRVIAHSHFMGRDITIELIRESLKDLLALQDKLVSVDNIQRT
VAEYYKIKISDLLSKRRSRSVARPRQVAMALSKELTNHSLPEIGDMFGGRDHTTVLHACRKINELKESDA
DIREDYKNLLRTLTT

```

In order to extract the information in the gbff files you can run the following command:
```
parse_gbff_parallel(n_cpus=1, input_directory="gbff_files_unzipped", output_directory="fasta_files")
```

### Constructing a blast database

```
make_blastdb.make_blastdb_parallel()
```

### Performing a blast search of the database

```
cluster_blast.cluster_blast(query_file = "query_cluster.fa", database_file= "fasta_files_combined/fasta_files_combined_master")
```

### Editing and filtering the tabular blast output
```
parse_blastp_helpers.filter_blasttable("results/cluster_blast_out.txt")
```

While the default function parses and separates certain values stored in the header, you can also filter any of the columns using a query string.
For example to filter the table to only include hits for PaaA that have a bitscore greater than 250 and hits for PaaB greater than 200 you can use the following command:
```
parse_blastp_helpers.filter_blasttable("results/cluster_blast_out.txt", filter_params="qseqid == 'PaaA' and bitscore >= 250 or qseqid == 'PaaB' and bitscore >= 200")
```
Any of the columns of the blast output file can be filtered in a similar manner

### Parsing the filtered tabular blast output to identify gene clusters

First obtain all unique genome accessions to be examined for gene clusters
```
accs_list = parse_blastp_helpers.fetch_accessions("results/blastout.csv")
```

Next, input those accessions to the following function:
```
parse_blastp_helpers.parseblastout2_parallel(accs_list,"results/blastout.csv", n_cpus=1, max_gene_dist=10000, min_cluster_number=1, gene_color_dict=None)
```
The ```max_gene_dist``` parameter sets a threshold of how far a gene must be away from a neighboring gene in a cluster to be excluded from that particular cluster.

The ```min_cluster_number``` parameter sets the minimum number of unique genes that must be present for the genes to be considered within a clusters

### Fetching metadata for the gene clusters
```
parse_blastp_helpers.fetch_metadata("results/cluster_hits.csv", ncbi_api_key="52553dfd9c090cfba1c3b28a45d8a648fd09")
```

### Making index files and extracting gene neighborhoods
While the results of the ```parse_blastout_parallel()``` function identifies the query genes within each genome it does not include information on the neighboring or interspersed genes that do not match to one of the query proteins. Further indexing and parsing is needed to obtain the gene neighborhoods.

First index files must be created using the following function:
```
parse_blastp_helpers.make_indexprot_parallel("results/cluster_hits.csv",fasta_file_directory= "fasta_files")
```

Next, the neighborhoods can be obtained by parsing the index files and original gbff files. The parameters ```features_upstream``` and ```features_downstream``` determine how many gene features upstream the first gene in the cluster, and downstream of the last gene in the cluster to fetch if available. In some cases there will be a limited amount of features returned due to the cluster being located towards the end of a contig.
```
parse_blastp_helpers.fetch_neighborhood_parallel("results/cluster_hits.csv", features_upstream=10,features_downstream=10)
```

The result is a large tab separated file that contains the following fields:

accession	assembly	title	feature_count_nhbr	cluster_len_nhbr	synteny_nhbr	synteny_alphabet_nhbr	synteny_dir_dist_nhbr	synteny_dir_nhbr	cluster_number	adj_coord_list	tared_adj_coord_list	nhbrhood_hit_list	nhbrhood_locus_tags	nhbrhood_old_locus_tags	nhbrhood_prot_ids	nhbrhood_prot_name	nhbrhood_prot_seq	clusterGC	genomeGC	diffGC	four_mer_freq_cluster	four_mer_freq_genome	four_mer_distance	cluster_seq	filename	biosample	number_of_hits	cluster_len	synteny	synteny_dir_dist	synteny_dir_pident	synteny_dir_evalue	synteny_dir_bitscore	synteny_dir_score	synteny_dir_length	synteny_dir	name	hit_list	old_locus_hit_list	protein_name_list	protein_id_list	pseudogene_list	query_list	coord_list	ncbi_graphics	contig	complete_genome	plasmid	has_overlap	duplicated	biosample_id	gi	isolation_src	env_biome	env_feature	host	host_sciname	strain	all_biosample_metadata	assembly_base	gtdb_tax	ncbi_tax	same_taxonomy	domain_gtdb	phylum_gtdb	class_gtdb	order_gtdb	family_gtdb	genus_gtdb	species_gtdb


| Fields  | Description |
| ------- | -------- |
| accession   | 1  |
| assembly   | 2  |
| feature_count_nhbr   | 3  |
| title   | 3  |
| title   | 3  |
| title   | 3  |
| title   | 3  |
| title   | 3  |
| title   | 3  |
| title   | 3  |
| title   | 3  |
| title   | 3  |
| title   | 3  |
| title   | 3  |
| title   | 3  |
| title   | 3  |
| title   | 3  |
| title   | 3  |
| title   | 3  |
| title   | 3  |
| title   | 3  |
| title   | 3  |
| title   | 3  |
| title   | 3  |
| title   | 3  |
| title   | 3  |
| title   | 3  |
| title   | 3  |
| title   | 3  |
| title   | 3  |
| title   | 3  |
| title   | 3  |
| title   | 3  |
| title   | 3  |
