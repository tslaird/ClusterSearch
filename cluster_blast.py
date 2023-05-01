import psutil
import os
import subprocess

#blast the db using blastp


def blastout2blasttable(input_file):
    blastout = pd.read_csv(input_file, sep='\t',low_memory = False, names=["qseqid","qgi","qacc","sseqid","sallseqid","sgi","sallgi","sacc","sallacc","qstart","qend","sstart","send","qseq","sseq","evalue","bitscore","score","length","qlen","slen","pident","nident","mismatch","positive","gapopen","gaps","ppos","frames","qframe","sframe","btop","stitle","salltitles","sstrand","qcovs","qcovhsp","qcovus"])
    blastout.to_csv(input_file,sep='\t', index=False)


def cluster_blast(query_file,database_file, eval=0.1, n_cpus=1, results_directory='results'):
    if not os.path.exists("results/"):
        os.mkdir("results/")
    subprocess.call(['blastp','-db',database_file,'-query',query_file,'-out',results_directory+"/"+'cluster_blast_out.txt','-outfmt','6 qseqid qgi qacc sseqid sallseqid sgi sallgi sacc sallacc qstart qend sstart send qseq sseq evalue bitscore score length qlen slen pident nident mismatch positive gapopen gaps ppos frames qframe sframe btop stitle salltitles sstrand qcovs qcovhsp qcovus','-max_target_seqs','50000','-num_threads',str(n_cpus),'-evalue',str(eval)])
    blastout2blasttable(input_file=results_directory+"/"+'cluster_blast_out.txt',)
