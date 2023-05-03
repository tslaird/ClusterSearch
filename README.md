

# sample workflow:

```
make_blastdb.make_blastdb_parallel()
cluster_blast.cluster_blast(query_file = "query_cluster.fa", database_file= "fasta_files_combined/fasta_files_combined_master")
parse_blastp_helpers.filter_blasttable("results/cluster_blast_out.txt")
parse_blastp_helpers.filter_blasttable("results/cluster_blast_out.txt", filter_params="qseqid == 'IacA' and bitscore >= 250 or qseqid == 'IacB' and bitscore >= 480 or qseqid == 'IacC' and bitscore >= 90 or qseqid == 'IacD' and bitscore >= 140 or qseqid == 'IacE' and bitscore >= 240")

accs_list = parse_blastp_helpers.fetch_accessions("results/blastout.csv")
parse_blastp_helpers.parseblastout2_parallel(accs_list,"results/blastout.csv", n_cpus=1, max_gene_dist=10000, min_cluster_number=1, gene_color_dict=None)

```
