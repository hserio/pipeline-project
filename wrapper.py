#Track 2 - Genome assembly
import os
import argparse
import sys
import logging

#function to parse command line arguments
def check_arg(args=None):
     parser = argparse.ArgumentParser(description='Write pattern count given pattern and text')
     parser.add_argument('-i', '--input',
                         help='path to input file',
                         nargs='+', #allow multiple input files
                         required=True
                         )
     parser.add_argument('-o', '--output',
                         help='output file name',
                         required=True
                         )
     return parser.parse_args(args)

#retrieve command line arguments and assign to variables
args = check_arg(sys.argv[1:])
infile = args.input
outfile = args.output

#create a directory for the whole pipeline
os.makedirs("PipelineProject_Hannah_Serio", exist_ok=True)
#go to the new directory
os.chdir("PipelineProject_Hannah_Serio")

#get logging set up
logging.basicConfig(filename='pipeline.log', level=logging.INFO, format='%(asctime)s - %(message)s')

#use Biopython to search NCBI to get the fasta file
from Bio import Entrez
from Bio import SeqIO
Entrez.email = "hserio@luc.edu"
handle = Entrez.efetch(db="nucleotide", id=["NC_006273.2"], rettype="fasta")
#get the list of SeqIO objects in fasta format
reference = SeqIO.parse(handle, "fasta")
#write reference to a fasta file
HCMV = "HCMV.fasta"
with open(HCMV, "w") as f:
    SeqIO.write(reference, f, "fasta")

#create a directory for index files
os.makedirs("index", exist_ok=True)
#go to the new directory for index
os.chdir("index")
#use bowtie2 to build index
os.system(f'bowtie2-build ../{HCMV} HCMV')
#get current working directory
index_path = os.getcwd() + '/HCMV'
#go back to main Pipeline directory
os.chdir("../")

#create a directory for fastq files
os.makedirs("fastq", exist_ok=True)

#function that finds the number of reads in a fastq file
def find_reads(file):
    with open(file) as f:
        count = sum(1 for line in f) // 4 #each read has 4 rows
    return count


#readcount = {} #create a dictionary that has file name as keys and read counts as values
fastq_pairs = [] #initialize list of fastq pairs
#loop through all input files to move them to one fastq directory
for file in infile:
    #open each input file and read contents
    with open(file) as input_file:
        #get file name
        file_name = os.path.basename(file)
        #get path to output file in fastq directory
        out_path = os.path.join("fastq", file_name)
        #open out file to write
        with open(out_path, 'w') as output_file:
            #copy original input file to new directory
            output_file.write(input_file.read())
        #separate 1 & 2 fastq files
        if file_name.endswith("_1.fastq"):
            read2 = file.replace("_1.fastq", "_2.fastq") #get the read 2 file name
            if read2 in infile:
                fastq_pairs.append((os.path.basename(file), os.path.basename(read2))) #append pair of fastq files as list
    #get the original read count of the file 
    if file.endswith("1.fastq"):
        readcount= find_reads(file) #find the initial read count and append
        logging.info(f'{os.path.basename(file)} had {readcount} read pairs before Bowtie2 filtering')

#currently still in Pipeline directory 
#go to the directory with fastq files
os.chdir("fastq")
#get path for fastq file directory
fastq_path = os.getcwd() + '/'
#go back to main Pipeline directory
os.chdir("../")

#create a directory for mapped files
os.makedirs("mapped", exist_ok=True)
#go to the new directory for mapping
os.chdir("mapped")
#map sequences to reference using bowtie2
for reads in fastq_pairs:
    fastq1, fastq2 = reads #get file names of both fastq files
    outputname = f'{os.path.basename(fastq1)[:-6]}' #isolate basenames for output file, removing the '.fastq' part
    fastq1_full = fastq_path + fastq1 #get full paths to fastq files
    fastq2_full = fastq_path + fastq2
    os.system(f'bowtie2 -x {index_path} -1 {fastq1_full} -2 {fastq2_full} -S {outputname}.sam --al-conc {outputname}.fastq')

#go back to main Pipeline directory
os.chdir("../")

#get the read count of aligned reads after mapping
for reads in fastq_pairs:
    fastq1, fastq2 = reads #get file names of both fastq files in pair
    mapped_q = f'{os.path.basename(fastq1)[:-6]}.1.fastq' #isolate basenames for output file with mapping suffix
    mapped_path = os.path.join("mapped", mapped_q)
    new_count = find_reads(mapped_path)
    logging.info(f'{mapped_q} had {new_count} read pairs after Bowtie2 filtering')

#go to the directory with mapped fastq files
os.chdir("mapped")

#get list of mapped fastq files for spades
mapped_fastq = []
for file in sorted(os.listdir()): #sort list so files are in correct order
    if file.endswith('.fastq'): #need same suffix
        mapped_fastq.append(file)

#use spades to create transcriptome
command = f"spades.py -k 77,99,127 -t 2 --only-assembler --pe-1 1 {mapped_fastq[0]} --pe-2 1 {mapped_fastq[1]} --pe-1 2 {mapped_fastq[2]} --pe-2 2 {mapped_fastq[3]} --pe-1 3 {mapped_fastq[4]} --pe-2 3 {mapped_fastq[5]} --pe-1 4 {mapped_fastq[6]} --pe-2 4 {mapped_fastq[7]} -o transcriptome"
os.system(command)
#write command to pipeline.log
logging.info(command)

#go to the directory with mapped fastq files
os.chdir("transcriptome")

#use BioPython to parse fasta and get the number of contigs with length > 1000 and total assembly length
contigs_count = 0 #initialize contig count
assembly_length = 0 #initialize assembly length
max_contig = "" #initialize max contig
with open(f"contigs.fasta") as fasta_file: #read contigs fasta file
    for record in SeqIO.parse(fasta_file, "fasta"):
        if len(record.seq) > 1000:
            contigs_count += 1 #add one to get count of contigs over 1000 bp
            assembly_length += len(record.seq) #sum basepair number for all contigs over 1000 bp
        if len(record.seq) > len(max_contig): #find the longest contig for blast
            max_contig = record.seq #

#write contigs to pipeline.log
logging.info(f"There are {contigs_count} contigs > 1000 bp in the assembly.")

#write assembly length to pipeline.log
logging.info(f"There are {assembly_length} bp in the assembly.")

#create fasta file with longest contig
longest_contig = "longest_contig.fasta"
with open(longest_contig, 'w') as f:
    f.write(f">longest_contig\n{max_contig}\n") #write the longest contig to the file in fasta format

#use NCBI's dataset function to get the dataset for Betaherpesvirinae subfamily (taxa id 10357)
os.system(f"datasets download virus genome taxon 10357 --refseq --include genome")
os.system(f"unzip ncbi_dataset.zip")

#make database with blast+ program
current_dir = os.getcwd()
input_fna = current_dir + "/ncbi_dataset/data/genomic.fna" #get full path to .fna file
os.system(f"makeblastdb -in {input_fna} -out Betaherpesvirinae -title Betaherpesvirinae -dbtype nucl")
#use blast+ to see if assembly aligns with other strains, include requested output options
os.system(f"blastn -query longest_contig.fasta -db Betaherpesvirinae -out myresults.csv -outfmt '6 qseqid sseqid sacc pident length qstart qend sstart send bitscore evalue stitle' -max_target_seqs 10")

#get correct header for log file
log_header = 'sacc\tpident\tlength\tqstart\tqend\tsstart\tsend\tbitscore\tevalue\tstitle'
logging.info(log_header)
#get output file
hits = "myresults.csv"
#add csv to log file
with open(hits, 'r') as result:
    for line in result:
        logging.info(line.strip())

#with open(outfile, 'w') as o:
#    o.write(str(current_dir))