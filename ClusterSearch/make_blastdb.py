import glob
import os
import concurrent.futures
import psutil
import subprocess
from random import shuffle
import itertools


def combine_files(file_list,list_number, output_dir= "fasta_files_combined"):
    with open(output_dir+'/combo_'+str(list_number), 'w') as outfile:
        for f in file_list:
            with open(f) as infile:
                outfile.write(infile.read())


def combine_files_parallel(input_dir="fasta_files",output_dir="fasta_files_combined",n_cpus=1,chunk_size=5000):
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
    files = glob.glob(input_dir+"/*.fa")
    print(files)
    file_chunks = [files[i * chunk_size:(i + 1) * chunk_size] for i in range((len(files) + chunk_size - 1) //chunk_size )]
    with concurrent.futures.ProcessPoolExecutor(max_workers=n_cpus) as executor:
        for _ in executor.map(combine_files, file_chunks,range(len(file_chunks)), itertools.repeat(output_dir) ):
            pass

def make_blastdb(file_list_index):
    combined_file = 'fasta_files_combined/combo_'+str(file_list_index)
    print("Making blastdb for: "+ combined_file)
    subprocess.call(['makeblastdb','-in',combined_file,'-dbtype','prot'])

def make_blastdb_parallel(input_dir="fasta_files_combined", n_cpus=1):
    combined_files_list = glob.glob(input_dir+"/combo_*")
    print(combined_files_list)
    with concurrent.futures.ProcessPoolExecutor(max_workers=n_cpus) as executor:
        for _ in executor.map(make_blastdb, range(len(combined_files_list))):
            pass
    combined_files_str = " ".join(combined_files_list)
    print(combined_files_str)
    subprocess.call(['blastdb_aliastool','-dblist',combined_files_str,'-dbtype','prot','-out',input_dir+'/fasta_files_combined_master','-title','all_proteins_combined_master'])
