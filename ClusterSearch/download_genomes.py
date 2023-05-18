import pandas as pd
import urllib.request
import concurrent.futures
import glob
import re
import sys
import psutil
import hashlib
import os
import gzip

# if os.path.exists("assembly_summary.txt"):
#     print("assembly summary already downloaded")
# else:
#     url = 'ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt'
#     urllib.request.urlretrieve(url,"assembly_summary.txt")
#
# if os.path.exists('assembly_summary_filtered.txt'):
#     assembly_summary_filtered= pd.read_csv("assembly_summary_filtered.txt", sep='\t')
# else:
#     assembly_summary = pd.read_csv('assembly_summary.txt', sep = '\t', low_memory = False, skiprows=[0])
#     #assembly_summary_filtered = assembly_summary[((assembly_summary['assembly_level'] == 'Complete Genome') & (assembly_summary['version_status']=='latest') & assembly_summary['organism_name'].str.contains('Pseudomonas putida'))]
#     assembly_summary_filtered = assembly_summary[(assembly_summary['version_status']=='latest')]
#     assembly_summary_filtered.to_csv("assembly_summary_filtered.txt", sep='\t')
# all_paths = assembly_summary_filtered['ftp_path']+ "/"+[i.split('/')[-1] for i in assembly_summary_filtered['ftp_path']] + \
#     "_genomic.gbff.gz"
# all_paths.to_csv("ftp_paths.txt", index=False)
# sample_names_df=pd.DataFrame([re.search("GCF_.+genomic(?=.gbff.gz)?",i.split('/')[-1])[0] for i in all_paths])
# sample_names_df.to_csv("sample_names.txt", sep='\t', index=False, header=False)
# sample_names=[re.search("GCF_.+genomic(?=.gbff.gz)?",i.split('/')[-1])[0] for i in all_paths]
# output_gbff_names= ['gbff_files/'+i.split('/')[-1] for i in all_paths]
# output_gbff_unzip_names= ['gbff_files_unzipped/'+i.split('/')[-1][0:-3] for i in all_paths]
# download_dict= dict(zip(output_gbff_names, all_paths))

# chunk_size=500
# fasta_names= ["fasta_files/"+i+"_proteins.fa" for i in sample_names]
# fasta_file_chunks= [fasta_names[i * chunk_size:(i + 1) * chunk_size] for i in range((len(fasta_names) + chunk_size - 1) // chunk_size )]
# combined_fasta_chunks_index=list(range(0,len(fasta_file_chunks)))

#number_of_cpus = psutil.cpu_count()
#print("using "+str(number_of_cpus)+" cpus")

def fetch_gbff_files(path, retries=10, unzip=True, redownload=True):
    name = 'gbff_files/'+path.split('/')[-1]
    output_name = 'gbff_files_unzipped/'+name.split('/')[1][:-3]
    md5path= path.rsplit('/',1)[0]+"/md5checksums.txt"
    md5name= 'gbff_files/'+path.split('/')[-1]+"_md5checksums.txt"
    if os.path.exists(name) and os.path.exists(md5name) and redownload==False:
        with open(name, 'rb') as file_to_check:
            data = file_to_check.read()
            md5_returned = str(hashlib.md5(data).hexdigest())
        with open(md5name,'r') as md5_file:
            md5_file_text=md5_file.read()
            md5_sum=str(re.findall('\w+(?=\s\s.+gbff.gz)',md5_file_text)[0])
            if md5_returned==md5_sum:
                print(name + " already downloaded")
    else:
        while(retries >= 1):
            try:
                urllib.request.urlretrieve(path,name)
                urllib.request.urlretrieve(md5path,md5name)
                with open(name, 'rb') as file_to_check:
                    data = file_to_check.read()
                    md5_returned = str(hashlib.md5(data).hexdigest())
                with open(md5name,'r') as md5_file:
                    md5_file_text=md5_file.read()
                    md5_sum=str(re.findall('\w+(?=\s\s.+gbff.gz)',md5_file_text)[0])
                    print(md5_returned+" : "+md5_sum)
                    if md5_returned==md5_sum:
                        print("Fetched " + name)
                        if unzip ==True:
                            f = gzip.open(name, 'rb')
                            file_content = f.read()
                            f.close()
                            print('Unzipping '+ name)
                            output = open(output_name, 'wb')
                            output.write(file_content)
                            output.close()
                            print("Done unzipping "+ name)
                        break
            except:
                print("Retrying download from " + path)
                retries = retries - 1
                continue
        if not os.path.exists(name) and not os.path.exists(md5name):
            print("Failed downloading "+ path)
        else:
            return

def unzip_file(name):
    output_name = 'gbff_files_unzipped/'+name.split('/')[1].split('.gz')[0]
    if os.path.exists(output_name):
        print(name + " already unzipped")
    else:
        f = gzip.open(name, 'rb')
        file_content = f.read()
        f.close()
        print('Unzipping '+ name)
        output = open(output_name, 'wb')
        output.write(file_content)
        output.close()
        print("Done unzipping "+ name)




def download_gbff_files_parallel(filter_params=None, num_cpus=4):
    if os.path.exists("assembly_summary.txt"):
        print("assembly summary already downloaded")
    else:
        url = 'ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt'
        urllib.request.urlretrieve(url,"assembly_summary.txt")
    if os.path.exists("gbff_files"):
        pass
    else:
        os.mkdir("gbff_files")
    if os.path.exists("gbff_files_unzipped"):
        pass
    else:
        os.mkdir("gbff_files_unzipped")
    if os.path.exists('assembly_summary_filtered.txt'):
        assembly_summary_filtered= pd.read_csv("assembly_summary_filtered.txt", sep='\t')
    else:
        assembly_summary = pd.read_csv('assembly_summary.txt', sep = '\t', low_memory = False, skiprows=[0])
        #assembly_summary_filtered = assembly_summary[((assembly_summary['assembly_level'] == 'Complete Genome') & (assembly_summary['version_status']=='latest') & assembly_summary['organism_name'].str.contains('Pseudomonas putida'))]
        assembly_summary_filtered = assembly_summary[(assembly_summary['version_status']=='latest')]
        if filter_params is not None:
            assembly_summary_filtered = assembly_summary.query(filter_params)
        assembly_summary_filtered.to_csv("assembly_summary_filtered.txt", sep='\t')
    all_paths = assembly_summary_filtered['ftp_path']+ "/"+[i.split('/')[-1] for i in assembly_summary_filtered['ftp_path']] + \
        "_genomic.gbff.gz"
    all_paths.to_csv("ftp_paths.txt", index=False)
    sample_names_df=pd.DataFrame([re.search("GCF_.+genomic(?=.gbff.gz)?",i.split('/')[-1])[0] for i in all_paths])
    sample_names_df.to_csv("sample_names.txt", sep='\t', index=False, header=False)
    sample_names=[re.search("GCF_.+genomic(?=.gbff.gz)?",i.split('/')[-1])[0] for i in all_paths]
    output_gbff_names= ['gbff_files/'+i.split('/')[-1] for i in all_paths]
    output_gbff_unzip_names= ['gbff_files_unzipped/'+i.split('/')[-1][0:-3] for i in all_paths]
    download_dict= dict(zip(output_gbff_names, all_paths))
    with concurrent.futures.ThreadPoolExecutor(max_workers=num_cpus) as executor:
        for _ in executor.map(fetch_gbff_files, all_paths):
            pass

#gbff_files_dir=glob.glob("gbff_files/*gbff.gz")
#print(str(len(gbff_files_dir))+ " gbff.gz files in gbff_files directory")

#gbff_files_unzip_dir=glob.glob("gbff_files_unzipped/*.gbff")
#print(str(len(gbff_files_unzip_dir))+ " gbff files in gbff_files_unzipped directory")

#files_not_fetched=list(set(output_gbff_names)- set(gbff_files_dir))
#if len(files_not_fetched) >0:
#    print(files_not_fetched)
#    print(str(len(files_not_fetched))+" files were not fetched")
#    paths_not_fetched=[]
#    for i in files_not_fetched:
#        paths_not_fetched.append(download_dict[i])
#    print("Paths not fetched")
#    print(paths_not_fetched)
#    paths_not_fetched_df=pd.DataFrame(paths_not_fetched)
#    paths_not_fetched_df.to_csv("paths_not_fetched.txt",sep='\t', index=False, header=False)
#
#else:
#    print("All files fetched")

#files_not_unzipped= ["gbff_files/"+i.split("/")[-1]+".gz" for i in list(set(output_gbff_unzip_names)- set(gbff_files_unzip_dir))]

#with concurrent.futures.ThreadPoolExecutor(max_workers=number_of_cpus) as executor:
#    for _ in executor.map(fetch_gbff_files,[download_dict[i] for i in files_not_unzipped]):
#        pass
