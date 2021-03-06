# PERJ
# Pair-End Reads Joining (PERJ) for joining paired-end reads into one fragment
Function:

Join two illumina paired-end reads from two files, with options of 5 end triming options, outputing in fasta (1) or fastq (2), filling in Ns between left and right reads, trim low quality reads. The output will be a file containing the joint peseudo reads.

#input data are the name of two fastq files and options

# How to use this tool
The PERJ can be run in command line or graphic interface. Either is working.

## Usage for command terminal:
perl PERJ.pl [Options]

Options:

     -l File:  Left reads file in fastq format

     -r File:  Right reads file in fastq format

     -format integer: 1 for output sequence in fasta or 2 for output sequence in fastq

     -trimleft integer:  trim off how many of bp of nucleotides for the left read, default: 0

     -trimright integer:  trim off how many of bp of nucleotides for the  right read, default: 0

     -n integer: the number of bp of nucleotide N inserted into the gap between left and right reads, default: 0, for standard Illumina reads value 200 is recommended

## Example: 
A testing data comes with PERJ and run the following command in the terminal
`perl PERJ.pl -l TTGGATGG_1_50.fq -r TTGGATGG_2_50.fq -trimleft 1 -trimright 2 -format 1 -n 0 -o testOUT.txt`

## Usage for graph interface
To use , user can use one of floowing way:

1. `java -jar PERJ.jar`
2. Double click the PERJ.jar to initiate the graphic version, then select the read files and run
Then you will see the following interface of PERJ
![What is this](PERJ_graphic.png)

## Manual: 
Please see PERJ_manual_en_V11.pdf

## Citing this tool
Peng Qi, Davis Gimode, Dipnarayan Saha, Stephan Schröder, Debkanta Chakraborty, Xuewen Wang, Mathews M Dida, Russell L Malmberg, Katrien M Devos, (2018), UGbS-Flex, a novel bioinformatics pipeline for imputation-free SNP discovery in polyploids without a reference genome: finger millet as a case study, BMC Plant Biology, 18:117, doi: 10.1186/s12870-018-1316-3, [full text](https://bmcplantbiol.biomedcentral.com/articles/10.1186/s12870-018-1316-3)

