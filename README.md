# Pipeline Project (Track 2)

## Dependencies

All of the tools below must be installed in order for the code to run. Versions to download are also listed.

Python - (3.10) Python is necessary to run the code.

BioPython - BioPython is used to parse fasta files and must be installed prior to running the script.

Bowtie2 - (2.5.2) Bowtie2 is used to create an index and map fastq files to a reference.

SPAdes genome assembler - (3.15.5) SPAdes is used to assemble a genome using mapped reads.

Nucleotide BLAST - (2.15.0+) BLAST is used to search for virus strain matches. This also includes 'makeblastdb' which is necessary for creating a database.

## Running the code

Included in this repository are several sample data files that can be used to test the functionality of the code. Multiple input fastq files can be added, as long as they are separated by a space and read 1 is before read 2.

````bash
python wrapper_final.py -i /path/to/input/_1.fastq /path/to/input/_2.fastq
````

A log file 'PipelineProject.log' will be generated with the top 10 hits.

To test the script using the test data provided in the repository, download the .fastq files and run the command below.

````bash
python wrapper_final.py -i donor1d2_1.fastq donor1d2_2.fastq donor1d6_1.fastq donor1d6_2.fastq donor3d2_1.fastq donor3d2_2.fastq donor3d6_1.fastq donor3d6_2.fastq
````
