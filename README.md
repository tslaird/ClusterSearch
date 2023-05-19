
# Installation

This package can be installed directly from github within a conda environment using the following commands:
```
conda create -n ClusterSearch
conda activate ClusterSearch
conda install git pip
pip install git+https://github.com/tslaird/ClusterSearch.git
```
ClusterSearch relies on the BLAST software.
The blast executables can be downloaded from (https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/)
or installed with conda using the following command:
```
conda install -c bioconda blast

```

Once blast is installed you can run the following workflow (modified as needed) as long as you have activated the ClusterSearch conda environment and are running a python interpreter in your shell by typing and executing ```python``` in a Terminal window.

Then you can import the different modules for each set of commands:
```
from ClusterSearch import download_genomes
from ClusterSearch import parse_gbff
from ClusterSearch import make_blastdb
from ClusterSearch import cluster_blast
from ClusterSearch import parse_blastp
```

### Sample workflow:

The following is a sample workflow that searches for the well characterized phenylacetic acid catabolic gene cluster within complete (fully assembled) **Pseudomonas putida** genomes in the RefSeq database. The commands are meant to be run in a python interpreter.

### Downloading gbff files
The first step is to download files from RefSeq.
This can be done manually via a web interface if desired however the ```download_gbff_files_parallel()``` function enables doing this from the command line. The function first downloads the RefSeq assembly summary file and then gets the ftp paths for each genome for subsequent download. If you are running the workflow on an HPC or with lots of storage you can certainly download all RefSeq genomes. However, the function takes a parameter ```filter_params``` which lets you subset the entire RefSeq dataset to only include geomes of interest (for example just from a particular genus or species). For the example here the filter will allow only complete genomes which contain the word "Pseudomonas putida" and that are the latest version of the genome entry. You can filter based on any column of the assembly summary file using the query syntax available in the pandas python package (see more here: https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.query.html).

```
download_gbff_files_parallel(filter_params="assembly_level == 'Complete Genome' and organism_name.str.contains('Pseudomonas putida')")
```

The files will be downloaded to a folder named ```gbff_files```

The files will then be unzipped to a folder named ```gbff_unzipped_files```

### Extracting protein information from gbff files

After downloading and extracting the gbff files the next step is to parse the files to construct fasta protein files for a subsequent blast search.
Each protein within the gbff file gets extracted along with relevant information which is stored in the fasta header.
A sample extracted protein is listed below. Each component of the header is separated by a delimeter of two exclamation points ```!!``` Interestingly it was only one exclamation point but the parser got confused because there was a genome name that contained an exclamation point. Some people give their genome entries strange names.
The header entries are ```>filename!!Assembly_id!!Nucleotide_accession!!Genome_name!!Locus_tag!!Old_locus_tag!!Biosample_id!!Protein_name!!DNA_coordinates!!Protein_id!!Pseudogene```

The last entry in the header will be NULL unless the gene was classified as a pseudogene by the RefSeq annotation.
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
parse_gbff.parse_gbff_parallel(n_cpus=4, input_directory="gbff_files_unzipped", output_directory="fasta_files")
```

### Constructing a blast database
First combine the genome files:
```
make_blastdb.combine_files_parallel(input_dir="fasta_files",output_dir="fasta_files_combined",n_cpus=1,chunk_size=5000)
```
Then make a blast database:
```
make_blastdb.make_blastdb_parallel()
```

### Performing a blast search of the database

To search for the Paa operon, first create a file in the working directory that contains the Paa protein sequences below. These sequences were obtained from the Microbes Online database from the *Escherichia coli* K-12 MG1655 genome. Name the file: ```Paa_operon_EcoliK-12_MG1655.fa```

The file text should be:
```
>PaaA
MTQEERFEQRIAQETAIEPQDWMPDAYRKTLIRQIGQHAHSEIVGMLPEGNWITRAPTLR
RKAILLAKVQDEAGHGLYLYSAAETLGCAREDIYQKMLDGRMKYSSIFNYPTLSWADIGV
IGWLVDGAAIVNQVALCRTSYGPYARAMVKICKEESFHQRQGFEACMALAQGSEAQKQML
QDAINRFWWPALMMFGPNDDNSPNSARSLTWKIKRFTNDELRQRFVDNTVPQVEMLGMTV
PDPDLHFDTESGHYRFGEIDWQEFNEVINGRGICNQERLDAKRKAWEEGTWVREAALAHA
QKQHARKVA
>PaaB
MSNVYWPLYEVFVRGKQGLSHRHVGSLHAADERMALENARDAYTRRSEGCSIWVVKASEI
VASQPEERGEFFDPAESKVYRHPTFYTIPDGIEHM
>PaaC
MNQLTAYTLRLGDNCLVLSQRLGEWCGHAPELEIDLALANIGLDLLGQARNFLSYAAELA
GEGDEDTLAFTRDERQFSNLLLVEQPNGNFADTIARQYFIDAWHVALFTRLMESRDPQLA
AISAKAIKEARYHLRFSRGWLERLGNGTDVSGQKMQQAINKLWRFTAELFDADEIDIALS
EEGIAVDPRTLRAAWEAEVFAGINEATLNVPQEQAYRTGGKKGLHTEHLGPMLAEMQYLQ
RVLPGQQW
>PaaD
MGMQRLATIAPPQVHEIWALLSQIPDPEIPVLTITDLGMVRNVTQMGEGWVIGFTPTYSG
CPATEHLIGAIREAMTTNGFTPVQVVLQLDPAWTTDWMTPDARERLREYGISPPAGHSCH
AHLPPEVRCPRCASVHTTLISEFGSTACKALYRCDSCREPFDYFKCI
>PaaE
MTTFHSLTVAKVESETRDAVTITFAVPQPLQEAYRFRPGQHLTLKASFDGEELRRCYSIC
RSYLPGEISVAVKAIEGGRFSRYAREHIRQGMTLEVMVPQGHFGYQPQAERQGRYLAIAA
GSGITPMLAIIATTLQTEPESQFTLIYGNRTSQSMMFRQALADLKDKYPQRLQLLCIFSQ
ETLDSDLLHGRIDGEKLQSLGASLINFRLYDEAFICGPAAMMDDAETALKALGMPDKTIH
LERFNTPGTRVKRSVNVQSDGQKVTVRQDGRDREIVLNADDESILDAALRQGADLPYACK
GGVCATCKCKVLRGKVAMETNYSLEPDELAAGYVLSCQALPLTSDVVVDFDAKGMA
>PaaF
MSELIVSRQQRVLLLTLNRPAARNALNNALLMQLVNELEAAATDTSISVCVITGNARFFA
AGADLNEMAEKDLAATLNDTRPQLWARLQAFNKPLIAAVNGYALGAGCELALLCDVVVAG
ENARFGLPEITLGIMPGAGGTQRLIRSVGKSLASKMVLSGESITAQQAQQAGLVSDVFPS
DLTLEYALQLASKMARHSPLALQAAKQALRQSQEVALQAGLAQERQLFTLLAATEDRHEG
ISAFLQKRTPDFKGR
>PaaG
MMEFILSHVEKGVMTLTLNRPERLNSFNDEMHAQLAECLKQVERDDTIRCLLLTGAGRGF
CAGQDLNDRNVDPTGPAPDLGMSVERFYNPLVRRLAKLPKPVICAVNGVAAGAGATLALG
GDIVIAARSAKFVMAFSKLGLIPDCGGTWLLPRVAGRARAMGLALLGNQLSAEQAHEWGM
IWQVVDDETLADTAQQLARHLATQPTFGLGLIKQAINSAETNTLDTQLDLERDYQRLAGR
SADYREGVSAFLAKRSPQFTGK
>PaaH
MMINVQTVAVIGSGTMGAGIAEVAASHGHQVLLYDISAEALTRAIDGIHARLNSRVTRGK
LTAETCERTLKRLIPVTDIHALAAADLVIEAASERLEVKKALFAQLAEVCPPQTLLTTNT
SSISITAIAAEIKNPERVAGLHFFNPAPVMKLVEVVSGLATAAEVVEQLCELTLSWGKQP
VRCHSTPGFIVNRVARPYYSEAWRALEEQVAAPEVIDAALRDGAGFPMGPLELTDLIGQD
VNFAVTCSVFNAFWQERRFLPSLVQQELVIGGRLGKKSGLGVYDWRAEREAVVGLEAVSD
SFSPMKVEKKSDGVTEIDDVLLIETQGETAQALAIRLARPVVVIDKMAGKVVTIAAAAVN
PDSATRKAIYYLQQQGKTVLQIADYPGMLIWRTVAMIINEALDALQKGVASEQDIDTAMR
LGVNYPYGPLAWGAQLGWQRILRLLENLQHHYGEERYRPCSLLRQRALLESGYES
>PaaI
MSHKAWQNAHAMYENDACAKALGIDIISMDEGFAVVTMTVTAQMLNGHQSCHGGQLFSLA
DTAFAYACNSQGLAAVASACTIDFLRPGFAGDTLTATAQVRHQGKQTGVYDIEIVNQQQK
TVALFRGKSHRIGGTITGEA
>PaaJ
MREAFICDGIRTPIGRYGGALSSVRADDLAAIPLRELLVRNPRLDAECIDDVILGCANQA
GEDNRNVARMATLLAGLPQSVSGTTINRLCGSGLDALGFAARAIKAGDGDLLIAGGVESM
SRAPFVMGKAASAFSRQAEMFDTTIGWRFVNPLMAQQFGTDSMPETAENVAELLKISRED
QDSFALRSQQRTAKAQSSGILAEEIVPVVLKNKKGVVTEIQHDEHLRPETTLEQLRGLKA
PFRANGVITAGNASGVNDGAAALIIASEQMAAAQGLTPRARIVAMATAGVEPRLMGLGPV
PATRRVLERAGLSIHDMDVIELNEAFAAQALGVLRELGLPDDAPHVNPNGGAIALGHPLG
MSGARLALAASHELHRRNGRYALCTMCIGVGQGIAMILERV
>PaaK
MITNTKLDPIETASVDELQALQTQRLKWTLKHAYENVPMYRRKFDAAGVHPDDFRELSDL
RKFPCTTKQDLRDNYPFDTFAVPMEQVVRIHASSGTTGKPTVVGYTQNDIDNWANIVARS
LRAAGGSPKDKIHVAYGYGLFTGGLGAHYGAERLGATVIPMSGGQTEKQAQLIRDFQPDM
IMVTPSYCLNLIEELERQLGGDASGCSLRVGVFGAEPWTQAMRKEIERRLGITALDIYGL
SEVMGPGVAMECLETTDGPTIWEDHFYPEIVNPHDGTPLADGEHGELLFTTLTKEALPVI
RYRTRDLTRLLPGTARTMRRMDRISGRSDDMLIIRGVNVFPSQLEEEIVKFEHLSPHYQL
EVNRRGHLDSLSVKVELKESSLTLTHEQRCQVCHQLRHRIKSMVGISTDVMIVNCGSIPR
SEGKACRVFDLRNIVGA
>PaaX
MSKLDTFIQHAVNAVPVSGTSLISSLYGDSLSHRGGEIWLGSLAALLEGLGFGERFVRTA
LFRLNKEGWLDVSRIGRRSFYSLSDKGLRLTRRAESKIYRAEQPAWDGKWLLLLSEGLDK
STLADVKKQLIWQGFGALAPSLMASPSQKLADVQTLLHEAGVADNVICFEAQIPLALSRA
ALRARVEECWHLTEQNAMYETFIQSFRPLVPLLKEAADELTPERAFHIQLLLIHFYRRVV
LKDPLLPEELLPAHWAGHTARQLCINIYQRVAPAALAFVSEKGETSVGELPAPGSLYFQR
FGGLNIEQEALCQFIR
>PaaY
MPIYQIDGLTPVVPEESFVHPTAVLIGDVILGKGVYVGPNASLRGDFGRIVVKDGANIQD
NCVMHGFPEQDTVVGEDGHIGHSAILHGCIIRRNALVGMNAVVMDGAVIGENSIVGASAF
VKAKAEMPANYLIVGSPAKAIRELSEQELAWKKQGTHEYQVLVTRCKQTLHQVEPLREIE
PGRKRLVFDENLRPKQ
>PaaZ
MQQLASFLSGTWQSGRGRSRLIHHAISGEALWEVTSEGLDMAAARQFAIEKGAPALRAMT
FIERAAMLKAVAKHLLSEKERFYALSAQTGATRADSWVDIEGGIGTLFTYASLGSRELPD
DTLWPEDELIPLSKEGGFAARHLLTSKSGVAVHINAFNFPCWGMLEKLAPTWLGGMPAII
KPATATAQLTQAMVKSIVDSGLVPEGAISLICGSAGDLLDHLDSQDVVTFTGSAATGQML
RVQPNIVAKSIPFTMEADSLNCCVLGEDVTPDQPEFALFIREVVREMTTKAGQKCTAIRR
IIVPQALVNAVSDALVARLQKVVVGDPAQEGVKMGALVNAEQRADVQEKVNILLAAGCEI
RLGGQADLSAAGAFFPPTLLYCPQPDETPAVHATEAFGPVATLMPAQNQRHALQLACAGG
GSLAGTLVTADPQIARQFIADAARTHGRIQILNEESAKESTGHGSPLPQLVHGGPGRAGG
GEELGGLRAVKHYMQRTAVQGSPTMLAAISKQWVRGAKVEEDRIHPFRKYFEELQPGDSL
LTPRRTMTEADIVNFACLSGDHFYAHMDKIAAAESIFGERVVHGYFVLSAAAGLFVDAGV
GPVIANYGLESLRFIEPVKPGDTIQVRLTCKRKTLKKQRSAEEKPTGVVEWAVEVFNQHQ
TPVALYSILTLVARQHGDFVD
```

Next, blast the query file against the created database using the following command
```
cluster_blast.cluster_blast(query_file = "Paa_operon_EcoliK-12_MG1655.fa", database_file= "fasta_files_combined/fasta_files_combined_master")
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
Any of the columns of the blast output file can be filtered in a similar manner.

### Parsing the filtered tabular blast output to identify gene clusters

First obtain all unique genome accessions to be examined for gene clusters
```
accs_list = parse_blastp_helpers.fetch_accessions("results/blastout.csv")
```

Next, input those accessions to the following function:
```
parse_blastp_helpers.parseblastout2_parallel(accs_list,"results/blastout.csv", n_cpus=4, max_gene_dist=10000, min_cluster_number=6, gene_color_dict=None)
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
parse_blastp_helpers.fetch_neighborhood_parallel("results/cluster_hits.csv", features_upstream=10, features_downstream=10, gene_prefix="Paa")
```

The result is a large tab separated file that contains the following fields:

| Fields  | Description |
| ------------- | -------------------------------------------------------- |
| accession | the accession number of the contig containing the cluster |
| assembly | the name of the RefSeq genome assembly the cluster is found in (e.g. "GCF_XXXXXXX.X) |
| title | the name of the contig containing the cluster |
| feature_count_nhbr | a count of genomic features in the cluster neighborhood |
| cluster_len_nhbr | the length in nucleotides of the cluster neighborhood |
| synteny_nhbr | text rendering of the cluster in the format GeneA-GeneB-x-GeneC (x represents genes not matching to query proteins)|
| synteny_alphabet_nhbr | text rendering of the cluster in the format GENEAgeneBGENEC (uppercase represents forward direction, lowercase represents reverse direction) |
| synteny_dir_dist_nhbr | text rendering of the cluster in the format \| A 〉-39-\| B 〉-77-〈 x \| -20-\| C 〉where the interspersed values are the distance between genes |
| synteny_dir_nhbr | text rendering of the cluster in the format \| A 〉\| B 〉〈 x \|\| C 〉 |
| cluster_number | a unique identifier for the cluster if there are multiple clusters on a particular genome or contig |
| adj_coord_list | the coordinate list of the cluster genes (start, stop, direction [1=fwd or -1 = reverse], query hit name) in the format of a python list of tuples adjusted so that a majority of genes are in the forward direction |
| tared_adj_coord_list | similar to the adj_coord_list but with the coordinates edited to start at 0 |
| nhbrhood_hit_list | a list of the hits matching the query proteins (hits not matching any query protein identified as "x")|
| nhbrhood_locus_tags | a list of RefSeq locus tags from proteins in the cluster neighborhood |
| nhbrhood_old_locus_tags | a list of RefSeq old locus tags from proteins in the cluster neighborhood |
| nhbrhood_prot_ids | a list of protein ids from proteins in the cluster neighborhood |
| nhbrhood_prot_name | a list of protein names from proteins in the cluster neighborhood |
| nhbrhood_prot_seq | a list of protein sequences from proteins in the cluster neighborhood (newline characters in the sequence were replaced with two question mark characters "??")|
| clusterGC | the GC content of the cluster |
| genomeGC | the GC content of the entire genome |
| diffGC | the difference in GC content between the cluster and entire genome |
| four_mer_freq_cluster | the tetranucleotide frequency of each of 256 possible tetranucleotides within the gene cluster (organized in a list according to the python function ```[''.join(p) for p in itertools.product(['a','t','g','c'], repeat=4)``` |
| four_mer_freq_genome | similar to four_mer_freq_cluster but of the entire genome excluding the cluster |
| four_mer_distance | the Manhattan distance calculated by scipy.spatial.distance.cityblock(four_mer_freq_cluster,four_mer_freq_genome) |
| cluster_seq | the DNA sequence of the cluster |
| filename | the local filename of the genome ".gbff" file |
| biosample | the biosample id of the genome |
| number_of_hits | the number of unique hits to query proteins within a cluster |
| cluster_length | the length in nucleotides from the beginning of the first gene in the cluster to the end of the last gene |
| synteny synteny_dir_dist | similar to synteny_dir_nhbr but lacking information on hits to non query proteins |
| synteny_dir_pident | text rendering of the cluster in the format \| A:82.0 〉\| B:99.0 〉\| C:50.0 〉where the numbers represent blast percent identity matches to the query protein |
| synteny_dir_bitscore | text rendering of the cluster in the format \| A:450 〉\| B:200 〉\| C:490 〉where the numbers represent blast bitscore values to the query protein |
| synteny_dir_length | text rendering of the cluster in the format \| A:230 〉\| B:120 〉\| C:320 〉where the numbers represent the length of the gene |
| synteny_dir | similar to synteny_dir_nhbr but lacking information on hits to non query proteins |
| name | duplicate of title |
| hit_list | a list of the "locus tags" for a cluster parsed from the gbff file |
| old_locus_hit_list | a list of the "old locus tags" for a cluster parsed from the gbff file |
| protein_name_list | a list of the protein names for a cluster parsed from the gbff file
| protein_id_list | a list of the protein ids for a cluster parsed from the gbff file |
| pseudogene_list | a list of all the hits that are pseudogenes in the gene neghborhood |
| query_list | the list of the hits to the query proteins |
| coord_list | the coordinate list of the cluster genes (start, stop, direction [1=fwd or -1 = reverse], query hit name) in the format of a python list of tuples |
| ncbi_graphics | a url link to the gene cluster overlay in the ncbi nucleotide graphics viewer |
| contig | whether or not the cluster resides on a contig (not on a complete genome) (based on pattern matching with the name field) |
| complete_genome | whether or not the cluster comes from a complete genome (based on pattern matching with the name field) |
| plasmid | whether or not the cluster is on a plasmid (based on pattern matching with the name field) |
| has_overlap | whether or not genes in the cluster are annotated as overlapping |
| duplicated | whether or not the cluster appears multiple times in a genome |
| biosample_id | ncbi biosample id |
| gi | ncbi gi number |
| isolation_src | the isolation source fetched from the biosample report (if available) |
| env_biome | the environmental biome fetched from the biosample report (if available) |
| env_feature | the environmental feature fetched from the biosample report (if available) |
| host | the host name fetched from the biosample report (if available) |
| host_sciname | the host scientific name fetched from the biosample report (if available) |
| strain | strain information fetched from the biosample report (if available) |
| all_biosample_metadata | all the metadata fetched from the biosampple report stored as a string resembling a python list of lists|
| assembly_base | the first component of the assembly string without the version number (example: GCF_000019125) |
| gtdb_tax | all gtdb taxonomy category classifications according to the gtdb database file (if available) |
| ncbi_tax | all ncbi taxonomy category classifications according to the gtdb database file (if available) |
| same_taxonomy | whether or not the gtdb and ncbi taxonomies are the same |
| domain_gtdb | the gtdb assigned domain (if available) |
| phylum_gtdb | the gtdb assigned phylum (if available) |
| class_gtdb | the gtdb assigned class (if available) |
| order_gtdb | the gtdb assigned order (if available) |
| family_gtdb | the gtdb assigned family (if available) |
| genus_gtdb | the gtdb assigned genus (if available) |
| species_gtdb | the gtdb assigned species (if available) |

### Downstream analysis

Often times you may want to do phylogenomic comparisons between certain genes within each gene cluster.
In order to better facilitate downstream comparisons you can use the function ```extract_cluster_prot```


For example to fetch the PaaA protein sequences from the analysis you could run:

```
extract_cluster_prot(hit_name="PaaA", neighborhood_file="results/cluster_positive_neighborhoods.tsv", outfile_name = "resulst/PaaA_hits.fa", additional_fields=['assembly','accession','title'], filter_params=None)
```

You could add additional info to each resulting fasta header by adding another one of the column names (such as species_gtdb or synteny) to the ```additional_fields``` list. Furthermore if you would like to subset the data to only include certain entries you can input a query string for the ```filter_params``` parameter which will subset the dataframe constructed from the neighborood_file.
