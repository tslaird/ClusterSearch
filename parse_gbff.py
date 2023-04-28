import re
import textwrap
import concurrent.futures
import glob
import os
import psutil
import pandas as pd

def gbff2protein(gbff_file):
    outname = [re.sub('.gbk|.gbff','_proteins.fa', gbff_file)]
    outname = ''.join(outname).split('/')[-1]
    if os.path.exists("fasta_files/"+outname):
        print(outname + " already converted")
    else:
        locus_re= re.compile('(?<=LOCUS\s{7})\w+')
        definition_re = re.compile('(?<=DEFINITION\s{2})[\s\S]+?(?=\nACC)')
        definition_sub = re.compile('\n|(\s{12})')
        locus_tag_re = re.compile('(?<=gene\W{12})[\S\s]+?locus_tag[\S\s]+?(?="\n)')
        locus_tag_sub = re.compile('[\s\S]+"')
        old_locus_tag_re = re.compile('(?<=gene\W{12})[\S\s]+?old_locus_tag[\S\s]+?(?="\n)')
        product_re = re.compile('(?<=product=")[\S\s]+?(?=")')
        coords_re = re.compile('(?<=(tRNA|rRNA|gene|ncRN)(\s|A){12})\S+')
        protein_id_re = re.compile('(?<=protein_id=")[\s\S]+?(?=")')
        protein_re = re.compile('(?<=translation=")[\S\s]+?(?=")')
        biosample_re = re.compile('(?<=BioSample:\s)\w+')
        assembly_re = re.compile('(?<=Assembly:\s)\S+')
        features_re = re.compile('(?<=\n)\W{5}(?=gene\W{12})')
        product_sub = re.compile('\n|\s{19}')
        protein_out_sub = re.compile('\n|(\s+)')
        separate_re = re.compile('//\n')
        pseudogene_re = re.compile('(?<=\/)pseudo(?=\n)')
        with open(gbff_file) as file:
            file_text=file.read()
        loci = separate_re.split(file_text)
        all_proteins = []
        for i in loci:
            if locus_re.search(i):
                locus = locus_re.search(i).group(0)
                definition = definition_re.search(i).group(0)
                definition = definition_sub.sub('',definition)
                biosample =  biosample_re.search(i)
                if biosample:
                    biosample_out = biosample.group(0)
                else:
                    biosample_out = 'NULL'
                if assembly_re.search(i):
                    assembly = assembly_re.search(i).group(0)
                else:
                    assembly = re.search('GCF_\d+?(?=\.)',gbff_file).group(0)
                features = features_re.split(i)
                for f in features[1:]:
                     locus_tag = locus_tag_re.search(f).group(0)
                     locus_tag = locus_tag_sub.sub('', locus_tag)
                     if old_locus_tag_re.search(f):
                         old_locus_tag = old_locus_tag_re.search(f).group(0)
                         old_locus_tag = locus_tag_sub.sub('', old_locus_tag)
                     else:
                         old_locus_tag = 'NULL'
                     if product_re.search(f):
                         product = product_re.search(f).group(0)
                         product = product_sub.sub('',product)
                     else:
                         product = 'NULL'
                     coords = coords_re.search(f).group(0)
                     protein_id = protein_id_re.search(f)
                     if pseudogene_re.search(f):
                         pseudogene = "PSEUDOGENE"
                     else:
                         pseudogene = "NULL"
                     if protein_id:
                         protein_id_out = protein_id.group(0)
                     else:
                         protein_id_out = 'NULL'
                     protein = protein_re.search(f)
                     if protein:
                         protein_out = protein.group(0)
                         protein_out = protein_out_sub.sub('', protein_out)
                         protein_out = textwrap.fill(protein_out, 70)
                         almost_all_items = ">%s!!%s!!%s!!%s!!%s!!%s!!%s!!%s!!%s!!%s!!%s" % (gbff_file.split("/")[-1],assembly,locus,locus_tag,old_locus_tag,definition,biosample_out,product,coords, protein_id_out, pseudogene)
                         almost_all_items = re.sub('\n|\s{4,}','',almost_all_items)
                         all_items = "%s\n%s\n" % (almost_all_items, protein_out)
                     else:
                         almost_all_items = ">%s!!%s!!%s!!%s!!%s!!%s!!%s!!%s!!%s!!%s!!%s" % (gbff_file.split("/")[-1],assembly,locus,locus_tag,old_locus_tag,definition,biosample_out,product,coords, protein_id_out, pseudogene)
                         almost_all_items = re.sub('\n|\s{4,}','',almost_all_items)
                         all_items = "%s\n" % (almost_all_items)
                     #all_proteins.add(">"+assembly+'•'+locus+'•'+locus_tag+'•'+biosample+'•'+product+'•'+coords+'•'+protein_id+'\n'+protein+'\n')
                     all_items= re.sub(' |,','_',all_items)
                     all_proteins.append(all_items)
        result= '\n'.join(all_proteins)+'\n'
        with open("fasta_files/"+outname,'w+') as output:
            output.writelines(result)
        print("extracted proteins for: " + gbff_file)


number_of_cpus = psutil.cpu_count()
print("using "+str(number_of_cpus)+" cpus")

if not os.path.exists("fasta_files"):
    os.mkdir("fasta_files")

gbff_inputs = glob.glob('gbff_files_unzipped/*.gbff')
print("Will process "+str(len(gbff_inputs))+" gbff files")

with concurrent.futures.ProcessPoolExecutor(max_workers=number_of_cpus*2) as executor:
    for _ in executor.map(gbff2protein,gbff_inputs):
        pass

fasta_files_dir=glob.glob("fasta_files/*fa")
print(str(len(fasta_files_dir))+ " fasta files in fasta_files directory")

output_names=[re.sub("gbff_files_unzipped","fasta_files",i) for i in gbff_inputs]
output_names=[re.sub('.gbk|.gbff','_proteins.fa', i) for i in output_names]

files_not_converted=list(set(output_names)-set(fasta_files_dir))

if len(files_not_converted) >0:
    print(str(len(files_not_converted))+" files were not converted")
    print("Files not converted")
    print(files_not_converted)
    files_not_converted_df=pd.DataFrame(files_not_converted)
    files_not_converted_df.to_csv("gbff_files_not_converted.txt",sep='\t', index=False, header=False)
else:
    print("All gbff files converted")
