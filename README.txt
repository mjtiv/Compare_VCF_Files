PROGRAM: Compare_VCF_Files

BRIEF DESCRIPTION: Program takes in genotyping and RNA-Seq results (VCF files),
filters the VCF files based on user defined parameters and compares
the overall genotyping calls for matching samples and reports all
the final statistics of the analysis. 

AUTHOR: M. Joseph Tomlinson IV
DATE CREATED: 4-01-2018

Program runs using python 3.6 and requires three imported python modules.

Required Python Modules:
sys
os.path
time

Compare_VCF_Files is a program that specifically takes in two VCF files. One VCF file
produced from genotyping data from Axiom software and another VCF file produced from 
GATK's best practices pipeline for variant calling. THe program filters both files based
on user parameters and then compares the files for matching variants. Once the matching variants
are identified the program then compares genotypes of matching samples. A reports are created of
matching variants and samples to examine the overall concordance of data.


Required Input Files (4 Files)

1. Genotyping VCF File 
	Genotyping file is a VCF file produced from Axiom software and for each sample only
	genotype calls are made. 

2. "Probeset_id" txt File
	A list of all passing probe set IDs created from the Axiom software for all
	passing variants.

3. RNA-Seq VCF File
	VCF file created from GATK's Best Practices for RNA-Seq for variant calling.
	Specifically its a vcf file with multiple samples, but same tissue

4. Matching Key File
	A matching key txt file that relates sample IDs from the genotyping file to
	sample IDs in RNA-Seq vcf files. Can easily be created in excel and saved as txt file.
	Header for first line is "Key" and "Value". First column is the genotyping panel sample names
	and second column is RNA-Seq sample names.


User Defined Parameters
Note: All parameters must start with the -- and "parameter_name" followed by the parameter
(actual parameter space seperated as seen below).  

--Exclusion_Region_Length NUMB_VALUE
	Exclusion length surrounding INDELS for filtering
	possible sequencing errors caused by INDEL alignment issues
	that could affect neighboring variants.

--Minimum_Quality_Score NUMB_VALUE
	The minimum Phred quality score for a variant to considered
	in the analysis

--Minimum_Read_Count NUMB_VALUE
	Minimum number of read counts for a sample (reference and alternative)
	for the sample to considered in the analysis procuess.

Note for Testing Purposes:
--Match_Limit NUMB_VALUE
	Limits the total number of matches identified in the two files being
	compared. Built in test if program is running correctly and also if
	input files are all in the correct format for later parsing. 


Major Output Files (created in same directory as program)

1. summary_report.txt 
	Summarizes the entire analysis process and gives summary statistics of the data

2. sample_report.txt 
	Summarizes the overall sample statistics (matches, no matches, no calls etc.). 
	The far right three columns summarize statistics about the no matches (NM) and breaks down
	the no matches into three categories; NM_Homozygous_Reference, NM_Homozygous_Alternative,
	NM_Discordant. Which breakdown the type of mis-match that can occur between the genotypes verus RNA-Seq.

3. variant_report.txt
	Summarizes the overall sample statistics (matches, no matches, no calls etc.). 
	The far right three columns summarize statistics about the no matches (NM) and breaks down
	the no matches into three categories; NM_Homozygous_Reference, NM_Homozygous_Alternative,
	NM_Discordant. Which breakdown the type of mis-match that can occur between the genotypes verus RNA-Seq.

4. summary_matrix_report.txt
	Sample matrix to show relationship between all samples. Samples should show high concordance with self
	and lower concordance with other samples (relationship of random chance versus true match). 

Log Directory
	A log directory exisits that contains various files created during intermediate steps of the anlaysis
	process. File names are self-explantory and files only should be used when trouble shooting. 


RUNNING PROGRAM

The program can be run on a local computer or on the command line. However, this program requires the parameter
file is found in the same directory as the program and all approriate parameters filled in. 
If running locally, run using the standard run command in a python gui. If running in an HPC environment
the program use the following command "python Compare_VCF_Files2.1.4.py". Version of the program may change,
so modification of the version number may need to occur. All output is produced in the same directory as the program
for simplicity purposes.  



