#!/usr/bin/env python3.6


"""
PROGRAM: Compare_VCF_Files
DESCRIPTION: Program takes in genotyping and RNA-Seq results (VCF files),
filters the VCF files based on user defined parameters and compares
the overall genotyping calls for matching samples and reports all
the final statistics of the analysis.

AUTHOR: M. Joseph Tomlinson IV

DATE CREATED: 3-19-2018
"""

# Required python modules
import sys
import os.path
import time


def split_variable(line):
    """
    Splits the line
    :param line: input line from file
    :return value: returns the line split as a tab
    """
    key,value = line.split("\t")
    return value


def read_parameter_file ():
    """
    Parsed apart a user input file to get location of files being analyzed,
    key file used to match sample IDs and parameters for the program
    param: NONE
    return: dictionary of various parameters for the program to run
    """

    input_file = open('Input_Parameter_File.txt', 'r') # Open the file in python

    testing_program_match_limit = 'off'
    
    for line in input_file:
        if line.startswith("--"):
            line=line.rstrip('\n')
            parsed_parameters = line.split("--")
            for parameters in range (1, len(parsed_parameters)):
                inputs = parsed_parameters[parameters].split(" ")

                if inputs[0] == "Genotyping_Array_File":
                    genotyping_file = inputs[1]

                elif inputs[0] == "600K_SNP_List":
                    list_600K_snps = inputs[1]

                elif inputs[0] == "RNA_Seq_File":
                    rna_seq_file = inputs[1]

                elif inputs[0] == "ID_Key_File":
                    id_key_file = inputs[1]

                elif inputs[0] == "Exclusion_Region_Length":
                    exclusion_region_length = inputs[1]

                elif inputs[0] == "Minimum_Read_Count":
                    min_total_read_count = inputs[1]

                elif inputs[0] == "Minimum_Quality_Score":
                    quality_score_min = inputs[1]

                elif inputs[0] == "Match_Limit":
                    testing_program_match_limit = inputs[1]
                    
        else:
            pass

    input_file.close() 
    
    return{'genotyping_file': genotyping_file, 'rna_seq_file': rna_seq_file,
           'id_key_file': id_key_file, 'exclusion_region_length': exclusion_region_length,
           'min_total_read_count': min_total_read_count, 'quality_score_min': quality_score_min,
           'list_600K_snps': list_600K_snps, 'testing_program_match_limit': testing_program_match_limit}


def get_file_name(input_file_name):
    """
    Parses apart the pathway from the raw file name, so the file name can be passed
    to the rest of the program as a simple variable
    :param input_file_name: full pathway of file being analyzed
    :return: the file name with pathway removed
    """

    input_file_name = input_file_name.split("\\")
    input_file_name = input_file_name[-1]
    input_file_name = input_file_name.split("/")
    file_name = input_file_name[-1]
    return (file_name)


#################################################################################

def filter_600K_Chip(file_name, snp_list_file):
    """
    Filters the 600K data for SNPs that pass user defined filters
    param file_name: name of the 600K SNP chip
    param snp_list: filtered SNP list created using Axiom software
    """

    microarray_stats={'raw_count': 0, 'snp_list_count': 0, 'kept_snp_count': 0, 'filtered_snp_count': 0,
                      'filtered_no_calls_count': 0}


    microarray_file_name = get_file_name (file_name)


    # Reading the SNP list
    snp_list_file = open (snp_list_file, "r")

    #SNP list to filter 600K
    snp_list = []
    for line in snp_list_file:
        if line.startswith("probeset_id"):
            pass
        else:
            line = line.rstrip('\r\n')
            snp = line
            microarray_stats['snp_list_count'] += 1
            snp_list.append(snp)
    snp_list_file.close()

    
    # Starting to Parse Files
    raw_600K_counts = count_variants(file_name) 
    
    input_file = open (file_name, "r")

    filtered_600K_file_name = "Log_Directory/Filtered_" + microarray_file_name

    filtered_600K_file = open (filtered_600K_file_name, "w")
    removed_600K_variants_file = open ("Log_Directory/Removed_" + microarray_file_name, "w")

    for line in input_file:

        microarray_stats['raw_count'] += 1
        
        #Reading Header and Printing to File
        if line.startswith("#CHROM"):
            filtered_600K_file.write(line)
            removed_600K_variants_file.write("Filter_Applied\t" + line)
        
        elif line.startswith(("##", "#", " #", "'#")):
            continue

        #Getting to actual variant data  
        else:
            parsed_line = line.split("\t")
            probe_id = parsed_line[2]
            probe_id = probe_id.rstrip('\r\n')
            
            if probe_id in snp_list:
                sample_data = parsed_line[9:]
                no_calls_counter = 0
                sample_total = len(sample_data)
                for sample in sample_data:
                    sample_genotype = sample.split(":")
                    if sample_genotype == "./.":
                        no_calls_counter +=1
                    else:
                        pass
                #finding all variants with no data recorded aka (./.)
                if no_calls_counter == sample_total:
                    removed_600K_variants_file.write("No_Calls\t" + line)
                    microarray_stats['filtered_no_calls_count'] += 1
                    
                else:
                    filtered_600K_file.write(line)
                    microarray_stats['kept_snp_count'] += 1

            else:
                removed_600K_variants_file.write("Filtered_SNPs\t" + line)
                microarray_stats['filtered_snp_count'] += 1
                 

    # Closing Files           
    input_file.close()
    filtered_600K_file.close()
    removed_600K_variants_file.close()

    return {'filtered_600K_file_name': filtered_600K_file_name, 'microarray_stats': microarray_stats}


#################################################################################
# RNA-Seq Filtering Code


def failure_Report(failure_File_Out, parsed_line, failure):
    """
    Writes all the failure variants to a seperate file with the exact
    failure defined to allow validation of filtering
    :param failure_File_Out: failure file being written to
    :param parsed_line: variant data being passed to program
    :param failure: exact failure for the variant
    :return: NONE
    """
    line=('\t'.join(map(str,parsed_line)))
    failure_File_Out.write(str(failure)+"\t"+line)
    return()


def removeSNPsInExclusionZone(parsed_line, no_of_exclusions, indel_exclusion_regions):
    """
    Identifes variants that are found in INDEL exclusion zones for filtering
    param parsed_line: variant being examined
    param noex: number of INDEL exclusion zones
    param indel_exclusion_regions: identified INDEL exclusion regions for data
    """
    #Judge Region Exclusion (local counter to pass or fail SNPs)
    judregex = 0
    #examining the overlap between the SNP and indel regions
    for i in range(no_of_exclusions): 
    #nested "if" loops to examine if SNP chromosomes match then compares the
    #the regions to see if SNP is Indel region (moves couter)
        if parsed_line[0] == indel_exclusion_regions[i][0]:
            if indel_exclusion_regions[i][1] < int(parsed_line[1]) < indel_exclusion_regions[i][2]:
                #judge if variant is in exclusion region counter, if fail will be
                #excluded from future analysis
                judregex += 1
                continue
    return{'line_of_data':parsed_line, 'judging_score':judregex}


def identify_INDEL_Regions(input_vcf, exclusion_region_length):
    """
    Identifies all INDELs in the dataset that could complicate variant analysis using information
    from both the reference and alternative alleles recorded in a vcf file
    :param inputvcf: input vcf file that is being analyzed for INDELs
    :param exreglen: exclusion region based on sequencing length of data
    :return: a dictionary with all the exclusion regions and total number of regions excluded
    """

    # infile being analyzed for INDELS
    vcf_input_file = open(input_vcf, "r")

    file_name = get_file_name(input_vcf)

    indel_log_file = open("Log_Directory/Identified_Indel_Regions_" + file_name, "w")
    
    ###Section of Code looks at the Indels That Passed and Counts Their Numbers

    #Creating a list of exclusion regions for Indels
    indel_exclusion_regions = []
    #Start Counter for no_of_exclusions
    no_of_exclusions = 0

    for line in vcf_input_file:
        
        #Get Header from file
        if line.startswith('#CHROM'):
            parsed_line = line.split('\t')
            header_info = parsed_line[:8]
            header_info = ('\t'.join(map(str,header_info)))
            indel_log_file.write("Chromosome\tStart_INDEL_Zone\tStop_INDEL_Zone\t" + header_info + "\n")

            #indel_log_file.write("Start\tStop\t
    
        elif line.startswith(("##", "#", " #", "'#")):
                    pass
                
        else:
            parsed_line  = line.rstrip().split('\t')
            #Is looking at the column called FILTER-- if "PASS" it passed all GATK filters, if failed
            # will list the filter it failed examples DP=filtered depth issue 
            if parsed_line[6] != 'PASS':  #This filter removes TONS of samples
                continue
            else:
                #Identifying the idels in the VCF file (columns 3 and 4)
                #Retrieve the alternative allele column from VCF
                alleleR = parsed_line[3].split(',') #splitting for different types of variants
                alleleA = parsed_line[4].split(',') #splitting for different types of variants
        
                #find the INDELs in either Ref or Alt using for loops
                Indel= 0
                for i in alleleR:
                    i = i.rstrip('"')
                    i = i.lstrip('"')
                    if len(i) > 1:
                        Indel += 1
                #Identify alternative alleles where indels or a deletion (asterix symbol)
                for j in alleleA:
                    j = j.rstrip('"')
                    j = j.lstrip('"')
                    if len(j) > 1 or j=="*":
                        Indel += 1
                        
                #Adds the Indel to a list and calculates the foward and reverse distance of it
                if Indel > 0:
                    indel_exclusion_regions.append([])
                    indel_exclusion_regions[no_of_exclusions].append(parsed_line[0])
                    indel_exclusion_regions[no_of_exclusions].append(int(parsed_line[1]) - exclusion_region_length)
                    indel_exclusion_regions[no_of_exclusions].append(int(parsed_line[1]) + exclusion_region_length)

                    #Counter for number of exclusions
                    no_of_exclusions += 1

                    #Write Info to Log file 
                    start_exclusion_zone = int(parsed_line[1]) - exclusion_region_length
                    stop_exclusion_zone = int(parsed_line[1]) + exclusion_region_length

                    indel_log_file.write(parsed_line[0] + "\t" + str(start_exclusion_zone)
                                         + "\t" + str(stop_exclusion_zone) +"\t")

                    variant_data = (parsed_line[:8])
                    variant_data = ('\t'.join(map(str,variant_data)))

                    indel_log_file.write(variant_data + "\n")
                                    
    vcf_input_file.close()  
    indel_log_file.close()
                    
    return {'indel_exclusion_regions':indel_exclusion_regions, 'no_of_exclusions':no_of_exclusions}


def testing_variant_counts(parsed_line, min_total_read_count):
    """
    Tests all samples for a variant to see if the variant can be futher analyzed.
    Various filtering criterias are applied to see if the variant is "testable."
    Can only examine mono-allelic or bi-allelic samples 
    param: parsed_line
    """

    #Variants starts as 'fail' if it makes it through the list will become 'pass'
    verdict = 'fail'

    #Local Counters---so if all samples show behavior triggers larger counter
    low_read_count_local_counter = 0
    low_freq_count_local_counter = 0
    no_value = 0
    
    #Testable Variants Counter
    testable_SNPs=0

    #starts the range of i to capture all samples genotype information
    #format for each sample genotype is (GT:AD:DP:GQ:PL)
    for sample in parsed_line[9:]:
        sample = sample.rstrip('"')
        sample = sample.lstrip('"')

        # Verify genotype data contains more information other files may not
        if sample.count(':'):
            
        #get the genotype information as a list for each sample (ex: ['0/1', '135,464', '605', '99', '11857,0,2154'])
            genotyping_information = sample.split(':')

            #Skipping all no recorded genotype values
            if genotyping_information[0] == './.':
                no_value +=1
                #print ("No Record")
                continue
    
            #get the alleles for a genotype (ex: ['0', '1'])
            alleles = genotyping_information[0].split('/')

            #Returning the physical genotype values, so can be later used to INDEX the data for
            #SNPs with 2 alternative alleles
            allele_one=(alleles[0])
            allele_two=(alleles[1])

            #get the count of the genotype counts (ex: ['135', '464'])
            genotyping_counts = genotyping_information[1].split(',')

            #Gets the allele total (ex: 135 + 464 = 599)
            total_counts = int(genotyping_counts[0])+int(genotyping_counts[1])

            #if the number of alleles is less than user input (standard value = 20) stop analysis
            if total_counts < min_total_read_count:
                low_read_count_local_counter+= 1
                #print ("Number of Reads Too Low")
                continue

            #Allele count for individuals [alt, ref] (ex: 464, 135)
            count_list = []

            # Note: alleles get loaded alternate first followed by reference
            if int(genotyping_counts[0]) > 0:
                count_list.append(int(genotyping_counts[0]))
                
            if int(genotyping_counts[1]) > 0:
                count_list.append(int(genotyping_counts[1]))
                
            # Finds the minimum in the list
            lowest_allele_count = min(count_list)

            # Checks to see minimum allele is < (1% of total alleles for bi-allelic samples
            # Issue with this trigger is monoallelic samples that get one biallelic count
            if lowest_allele_count <= (0.01*total_counts) and allele_one != allele_two:
                low_freq_count_local_counter+=1
                #print ("Lowest Allele Count Too Low")
                #print (sample)
                continue
      
            # sample passed all thresholds---passes and can be tested
            verdict = 'pass'
                                  
    return {'verdict':verdict,
            'low_read_count_local_counter':low_read_count_local_counter,
            'low_freq_count_local_counter':low_freq_count_local_counter,
            'no_value': no_value} 


def printing_filtering_final_data(filtered_rna_seq_file, parsed_line, min_total_read_count):
    """
    Prints final sample data to file, but after filtering genotypes with low
    read counts for a sample
    param filtered_rna_seq_file: filtered file being written to
    param parsed_line: variant data from file
    Return: NONE
    """

    variant_info = parsed_line[:8]
    variant_info = ('\t'.join(map(str,variant_info)))
    filtered_rna_seq_file.write(variant_info+"\t")
    filtered_rna_seq_file.write("Genotype/Filter:Counts/Geno:Counts\t")

    for sample in parsed_line[9:]:
        sample = sample.rstrip('"')
        sample = sample.lstrip('"')

        # Verify genotype data contains more information other files may not
        if sample.count(':'):

            #get the genotype information as a list for each sample (ex: ['0/1', '135,464', '605', '99', '11857,0,2154'])
            genotyping_information = sample.split(':')

            #Skipping all no recorded genotype values
            if genotyping_information[0] == './.':
                # Write to file and move on
                filtered_rna_seq_file.write(str(genotyping_information[0])+"\t")
                continue

            #get the alleles for a genotype (ex: ['0', '1'])
            alleles = genotyping_information[0].split('/')

            #Returning the physical genotype values, so can be later used to INDEX the data for
            #SNPs with 2 alternative alleles
            allele_one=(alleles[0])
            allele_two=(alleles[1])

            #get the count of the genotype counts (ex: ['135', '464'])
            genotyping_counts = genotyping_information[1].split(',')
            
            #Gets the allele total (ex: 135 + 464 = 599)
            total_counts = int(genotyping_counts[0])+int(genotyping_counts[1])

            #if the number of alleles is less than user input (standard value = 20) stop analysis
            if total_counts < min_total_read_count:
                filtered_rna_seq_file.write("Filter:"+str(genotyping_information[0])+":"+str(genotyping_information[1])+"\t")
                continue

            #Allele count for individuals [alt, ref] (ex: 464, 135)
            count_list = []

            # Note: alleles get loaded alternate first followed by reference
            if int(genotyping_counts[0]) > 0:
                count_list.append(int(genotyping_counts[0]))
                
            if int(genotyping_counts[1]) > 0:
                count_list.append(int(genotyping_counts[1]))
                
            # Finds the minimum in the list
            lowest_allele_count = min(count_list)

            # Checks to see minimum allele is < (1% of total alleles for bi-allelic samples
            # Issue with this trigger is monoallelic samples that get one biallelic count
            if lowest_allele_count <= (0.01*total_counts) and allele_one != allele_two:
                filtered_rna_seq_file.write("Filter"+"\t")
                continue

            filtered_rna_seq_file.write(str(genotyping_information[0])+":"+str(genotyping_information[1])+"\t")
    filtered_rna_seq_file.write("\n")
               

def filter_RNA_Seq_Data(input_file, exclusion_region_length, quality_score_min, min_total_read_count):
    """
    Filters the raw RNA-Seq data using various parameters to create a testable filtered dataset that
    will be matched with the other VCF file to identify common variants
    Param input_file: input RNA-Seq VCF file being filtered
    Param exclusion_region_length: length of the region surrounding indels to filter neighboring SNPs
    Param quality_score_min: minimum quality score to filter a variant
    Param min_total_read_count: minimum per sample read count
    Return filtered_rna_seq_file_name: name of the filtered file to use later in the program
    Return raw_rna_seq_stats: stats about the overall filtering process  
    """

    # Get the RNA_Seq File Name With Extra Pathway Stuff
    rna_seq_file = input_file
    rna_seq_file_name = get_file_name (rna_seq_file)

    # Open file
    rna_seq_file = open (input_file, 'r')


    #Passed Variants to be Tested
    filtered_rna_seq_file_name = "Log_Directory/Filtered_" + rna_seq_file_name
    
    filtered_rna_seq_file = open (filtered_rna_seq_file_name, 'w')

    #Log file of variants that fail and could NOT be tested (quality controls)
    failure_File_Out=open("Log_Directory/Failing_RNA_Seq_variants.txt",'w')

    
    #File Stats Dictionary
    raw_rna_seq_stats = {'variants_in_file': 0,
                        'failed_GATK_SNP_filter' : 0, 'failed_Qual_filter': 0,
                        'numb_SNPs_excl_Indels' : 0,
                        'ref_allele_mutiple_forms' : 0, 'alt_allele_mutiple_forms': 0,
                        'no_genotype_values': 0, 'low_read_count': 0, 'low_freq_count': 0,
                        'combo_filter_failure': 0,
                        'no_indel_exclusion_regions': 0,
                        'passing_variants': 0}
    

    #Identifies the INDEL regions in the dataset based on a certain length
    indel_data=identify_INDEL_Regions(input_file, exclusion_region_length)
    indel_exclusion_regions=indel_data['indel_exclusion_regions']
    no_of_exclusions=indel_data['no_of_exclusions']

    raw_rna_seq_stats['no_indel_exclusion_regions'] = no_of_exclusions

    #Starts parsing through the data
    for line in rna_seq_file:

        #Get the header from the file and print to new file
        if line.startswith("#CHROM"):
            file_header_info=line
            filtered_rna_seq_file.write(file_header_info)

            failure_File_Out.write("Failure\t"+line)

            #Getting sample IDs from list
            line = line.rstrip('\n')
            line = line.split("\t")
            sample_list=line[9:]
            
        elif line.startswith(("##", "#", " #", "'#")):
                    pass
        else:
            #Counting the Number of Probes in File
            raw_rna_seq_stats['variants_in_file']+= 1

            #Splitting the line
            parsed_line = line.split("\t")

            #Removes all SNPs without a filter passing score
            if parsed_line[6] != 'PASS':
                raw_rna_seq_stats['failed_GATK_SNP_filter']+=1
                failure="GATK_Filter"
                failure_Report(failure_File_Out, parsed_line, failure)
                continue
            
            #Quality score filter of data (QUAL column)
            if float(parsed_line[5]) < quality_score_min:
                raw_rna_seq_stats['failed_Qual_filter']+=1
                failure="Qual_Score_Filter"
                failure_Report(failure_File_Out, parsed_line, failure) 
                continue

            #Filters all SNPs found in INDEL Exclusion Zone
            exclusion_filter = removeSNPsInExclusionZone(parsed_line, no_of_exclusions, indel_exclusion_regions)
            parsed_line = exclusion_filter['line_of_data']
            judregex=exclusion_filter['judging_score']
            if judregex > 0:
                raw_rna_seq_stats['numb_SNPs_excl_Indels']+=1
                failure="Indel_Region"
                failure_Report(failure_File_Out, parsed_line, failure)
                continue

            #Splitting the Reference Alleles
            alleles_Ref = parsed_line[3].split(',')

            #Filters if Reference Allele has more than one form by chance that INDEL filter didn't catch
            # should occur, but built in for safety
            no_Ref_Alleles = len(alleles_Ref)
            if no_Ref_Alleles > 1:
                raw_rna_seq_stats['ref_allele_mutiple_forms']+= 1
                failure="Multiple_Ref_Alleles"
                failure_Report(failure_File_Out, parsed_line, failure)
                
                continue

            #Splitting the Alternative Alleles
            alleleA = parsed_line[4].split(',')

            #Filters alternative allele if greater than input threshold alleles
            noAtl = len(alleleA)
            if noAtl > 1:
                raw_rna_seq_stats['alt_allele_mutiple_forms']+=1
                failure="To_Many_Alt_Alleles"
                failure_Report(failure_File_Out, parsed_line, failure)
                continue
            
            # Testing if variant can actually be analyzed
            variant_counts=testing_variant_counts(parsed_line, min_total_read_count)

            #Various filtering Variables
            low_read_count_local_counter = variant_counts['low_read_count_local_counter']
            low_freq_count_local_counter = variant_counts['low_freq_count_local_counter']
            no_value_local_counter = variant_counts['no_value']

            #True or False statement for if variant passes and should be examined further
            verdict = variant_counts['verdict']

            #Recording the type of failure of data (most SNPs fail for a variety of reasons)
            noind = len(parsed_line[9:])
            if verdict == 'fail':

                if no_value_local_counter == noind:
                    raw_rna_seq_stats['no_genotype_values'] +=1
                    failure="No_Genotype_Value"
                    failure_Report(failure_File_Out, parsed_line, failure)
                    continue
                
                if low_read_count_local_counter == noind:
                    raw_rna_seq_stats['low_read_count'] +=1
                    failure="Samples_Low_Read_Count"
                    failure_Report(failure_File_Out, parsed_line, failure)
                    continue

                if low_freq_count_local_counter == noind:
                    raw_rna_seq_stats['low_freq_count']+=1
                    failure="Samples_Low_Freq"
                    failure_Report(failure_File_Out, parsed_line, failure)
                    continue
            
                if (low_read_count_local_counter + low_freq_count_local_counter + no_value_local_counter)== noind:
                    raw_rna_seq_stats['combo_filter_failure']+=1    
                    failure="Sample_Combo_Filter"
                    failure_Report(failure_File_Out, parsed_line, failure)
                    continue


            # Final Results Print to File for further testing
            else:
                raw_rna_seq_stats['passing_variants']+=1
                printing_filtering_final_data(filtered_rna_seq_file, parsed_line, min_total_read_count)


    #Flusing Data if writing slow down
    failure_File_Out.flush()
    rna_seq_file.flush()
    
    #Closing Files            
    rna_seq_file.close()
    filtered_rna_seq_file.close()
    failure_File_Out.close

    return {'filtered_rna_seq_file_name': filtered_rna_seq_file_name, 'raw_rna_seq_stats':raw_rna_seq_stats}


def count_variants(file_name):
    """
    Counts the number of variants in the file
    Param file_name: file to have it lines counted
    Return counts: number of variants in file
    """
    input_file = open (file_name, "r")

    counts = 0
    
    for line in input_file:
        #Skipping all extra information parts of file
        if line.startswith(("##", "#", " #", "'#")):
            pass
        #Getting to actual variant data  
        else:
            #print ("Next outer loop")
            
            #Counting the Number of Probes in File
            counts+= 1
    return (counts)



def compare_and_parse(input_file_one, input_file_two, testing_program_match_limit):
    """
    Compares and parses through the two genotyping files to find out which variants match (chromosomes and alleles)
    and prints the common variants to new files and also prints matching variants but with non-matching alleles to
    another file for further investigation of the discordance
    Param input_file_one: file one being compared
    Param input_file_two: file two being compared
    Return common_snps_file_one_name: name of the new file one with common variants
    Return common_snps_file_two_name: name of the new file two with common variants
    Return compare_and_parse_stats: statistics about the matching of the variant files
    """

    print ("Input File One")
    print (input_file_one)

    print ("Input File Two")
    print (input_file_two)
    

    #Setup Tallying Dictionary
    compare_and_parse_stats = {'input_file_one': 0,
                        'input_file_two' : 0, 'match': 0,
                        'non_matching_ref' : 0,
                        'non_matching_alt' : 0, 'both_non_matching': 0}
    

    count_file_one = count_variants(input_file_one)
    compare_and_parse_stats['input_file_one'] = count_file_one
    
    count_file_two = count_variants(input_file_two)
    compare_and_parse_stats['input_file_two'] = count_file_two

    # Get File names without pathway included
    print ("Running Program to Find Common SNPs Between VCF Files")
    print ("Files being Analyzed are:")
    
    file_one_name = get_file_name (input_file_one)
    
    print (file_one_name)
    
    file_two_name = get_file_name (input_file_two)
    print (file_two_name)

    # Creating Common SNPs Output Files
    common_snps_file_one_name = "Log_Directory/Common_" + file_one_name
    common_snps_file_two_name = "Log_Directory/Common_" + file_two_name
    
    common_snps_file_one = open ("Log_Directory/Common_" + file_one_name, "w")
    common_snps_file_two = open ("Log_Directory/Common_" + file_two_name, "w")
    non_matching_alleles_file = open("Log_Directory/Common_Variants_Non_Matching_Alleles_.txt", "w")

    # Header for non-matching report
    non_matching_alleles_file.write("Issue\tFile_One_Info\tFile_Two_Info\n")

    # Creating Dictionaries to Store Data
    file_one_stats_dict={}
    file_two_stats_dict={}

    # Input Files
    input_file_one= open(input_file_one, 'r')
    input_file_two= open(input_file_two, 'r')

    #################TESTING PURPOSES########
    matching_probes_count = 0

    header_flag_one = "true"
    header_flag_two = "true"
    #Creating nested loops to sort through data

    #Getting information from input file one
    for file_one_line in input_file_one:
        if file_one_line.startswith("#CHROM"):
            file_one_header_info=(file_one_line)
            if header_flag_one == "true":
                common_snps_file_one.write(file_one_header_info)
                header_flag_one = "false"
            else:  
                pass
        #Capturing possible VCF variations
        elif file_one_line.startswith(("##", "#", " #", "'#")):
            pass
        #Getting to actual variant data  
        else:
            #print ("Next outer loop")
            
            #Counting the Number of Probes in File
            file_one_variants_data = file_one_line.split("\t")
            
            file_one_variant=str(file_one_variants_data[0]+":"+file_one_variants_data[1])
            #print (file_one_variant)

            #Second nested loop of data
            #With each loop need to reset the file to Zero
            input_file_two.seek(0)
            for file_two_line in input_file_two:
                if file_two_line.startswith("#CHROM"):
                    file_two_header_info=(file_two_line)
                    if header_flag_two == "true":
                        common_snps_file_two.write(file_two_header_info)
                        header_flag_two = "false"
                    else:  
                        pass
                elif file_one_line.startswith(("##", "#", " #", "'#")):
                    pass
                else:
                    file_two_variants_data = file_two_line.split("\t")
                    file_two_variant=str(file_two_variants_data[0]+":"+file_two_variants_data[1])

                    # If Perfect Match for the Variants
                    if (file_one_variant==file_two_variant
                        and file_one_variants_data[3]==file_two_variants_data[3]
                        and file_one_variants_data[4]==file_two_variants_data[4]):
                        
                        common_snps_file_one.write(file_one_line)
                        common_snps_file_two.write(file_two_line)
                        #after matching data found, no point to continue loop
                        compare_and_parse_stats['match']+= 1
                       #print ("Match")
                        matching_probes_count+=1
                        break

                    # If Ref Allele Does Not Match for Variants
                    elif (file_one_variant==file_two_variant
                        and file_one_variants_data[3]!= file_two_variants_data[3]
                        and file_one_variants_data[4] == file_two_variants_data[4]):

                        file_one_variant_info=str(file_one_variants_data[0] + ":"
                                                  + file_one_variants_data[1] + ":"
                                                  + file_one_variants_data[2] + ":"
                                                  + file_one_variants_data[3] + ":"
                                                  + file_one_variants_data[4])

                        file_two_variant_info=str(file_two_variants_data[0] + ":"
                                                  + file_two_variants_data[1] + ":"
                                                  + file_two_variants_data[2] + ":"
                                                  + file_two_variants_data[3] + ":"
                                                  + file_two_variants_data[4])
                                                  
                        non_matching_alleles_file.write("Non_Match_Ref_Allele\t" + file_one_variant_info + "\t" + file_two_variant_info + "\n") 
                        #after matching data found, no point to continue loop
                        compare_and_parse_stats['non_matching_ref']+= 1
                        break

                    # If Alt Allele Does Not Match for Variants
                    elif (file_one_variant==file_two_variant
                        and file_one_variants_data[3] == file_two_variants_data[3]
                        and file_one_variants_data[4] != file_two_variants_data[4]):

                        file_one_variant_info=str(file_one_variants_data[0] + ":"
                                                  + file_one_variants_data[1] + ":"
                                                  + file_one_variants_data[2] + ":"
                                                  + file_one_variants_data[3] + ":"
                                                  + file_one_variants_data[4])

                        file_two_variant_info=str(file_two_variants_data[0] + ":"
                                                  + file_two_variants_data[1] + ":"
                                                  + file_two_variants_data[2] + ":"
                                                  + file_two_variants_data[3] + ":"
                                                  + file_two_variants_data[4])
                                                  
                        non_matching_alleles_file.write("Non_Match_Alt_Allele\t" + file_one_variant_info + "\t" + file_two_variant_info + "\n") 
                        #after matching data found, no point to continue loop
                        compare_and_parse_stats['non_matching_alt']+= 1
                        break

                    # If Alt Allele Does Not Match for Variants
                    elif (file_one_variant==file_two_variant
                        and file_one_variants_data[3] != file_two_variants_data[3]
                        and file_one_variants_data[4] != file_two_variants_data[4]):

                        file_one_variant_info=str(file_one_variants_data[0] + ":"
                                                  + file_one_variants_data[1] + ":"
                                                  + file_one_variants_data[2] + ":"
                                                  + file_one_variants_data[3] + ":"
                                                  + file_one_variants_data[4])

                        file_two_variant_info=str(file_two_variants_data[0] + ":"
                                                  + file_two_variants_data[1] + ":"
                                                  + file_two_variants_data[2] + ":"
                                                  + file_two_variants_data[3] + ":"
                                                  + file_two_variants_data[4])
                                                  
                        non_matching_alleles_file.write("Both_Non-Matching\t" + file_one_variant_info + "\t" + file_two_variant_info + "\n") 
                        #after matching data found, no point to continue loop
                        compare_and_parse_stats['both_non_matching']+= 1
                        break

                    #Keep Going Through Loop Until Match Found     
                    else:
                        pass

        ##########Testing Purposes---to Limit Loops############
        # Just test if a number, should be a string otherwise 'off'
                if testing_program_match_limit.isdigit():
                    if matching_probes_count == int(testing_program_match_limit):
                        break     
        ################Testing Purposes#######################
                
    common_snps_file_one.close()
    common_snps_file_two.close()
    non_matching_alleles_file.close()
    
    return{'common_snps_file_one_name': common_snps_file_one_name , 'common_snps_file_two_name':  common_snps_file_two_name,
           'compare_and_parse_stats': compare_and_parse_stats }

    
########################################################################################

def parse_key_file(key_dict_file):
    """ Parses the sample key file created by the user
        to create a key dictionary of sample names between
        the two groups
        param: key_dict_file (file provided by the user)
        return: sample_id_dict (key dictionary of sample ids between files)
    
    """
    # setting up sample ID dictionary
    sample_id_key_dict = {}
    tallying_sample_dict = {}

    #Cycling through the key dictionary file for IDs
    key_dict_file = open(key_dict_file, "r")
    for line in key_dict_file:
        if line.startswith('Key'):
            pass
        else:
            line = line.rstrip('\r\n')
            key,value = line.split('\t')
            sample_id_key_dict.update({key:value})
            # Tallying Dictionary
            # Sample: Match, One_Allele_Expressed, No_Match
            tallying_sample_dict.update({key:{'match':0, 'no_match': 0, 'total_counts': 0, 'homo_ref': 0, 'homo_alt': 0, 'discordant': 0,
                                              'no_call_DNA': 0, 'no_call_RNA': 0, 'no_call_Both': 0, 'total_no_calls': 0}})
    
    key_dict_file.close()

    return{'sample_id_key_dict': sample_id_key_dict, 'tallying_sample_dict': tallying_sample_dict}


def samples_data (file_name):
    """
    Takes in a VCF file and returns a dictionary of samples with corresponding genotype data
    allowing for the samples genotypes to be compared and tabulated (assuming other file is in same order)
    param: input file (VCF file)
    return: sample dictionary with values corresponding to sample genotypes
    return: variant_list
    """
    file = open(file_name, 'r')

    #Variant Dictionary
    variant_dict = {}
    variant_counter = 0

    # Sample_List
    sample_list = []
    
    #Starts looping through file
    for line in file:
        if line.startswith("#CHROM"):
            line = line.rstrip("\r\n")
            headers=line.split("\t")
            samples_IDs = headers[9:]

            #Create a dictionary of samples to put genotype calls into
            file_samples_genotype_dict={}
            for sample in samples_IDs:
                sample_list.append(sample)
                sample = sample.rstrip('\r\n')
                file_samples_genotype_dict.update({sample: []})
            #print (file_samples_genotype_dict)
            
        else:
            line = line.rstrip("\r\n")
            # Need to remove extra tab put in filtering step
            line = line.rstrip("\t")
            line = line.split("\t")

            # Dealing with variant data
            # Variant data from line
            variant_info=line[0:5]
            # Putting variant information into a dictionary
            key = variant_counter
            variant_dict.update({key:{'variant_info':variant_info,
                                      'match':0, 'no_match': 0, 'total_counts': 0, 'homo_ref': 0,
                                      'homo_alt': 0, 'discordant': 0, 'no_call_DNA': 0,
                                      'no_call_RNA': 0,'no_call_Both': 0, 'total_no_calls': 0 }})
            #Moving Variant Counter
            variant_counter+=1

            # Dealing with Genotyping Data from Samples
            #Genotype data from line
            genotypes = line[9:]
            #print (genotypes)

            #print (genotypes)
            index_position = 0  
            for sample_genotype in genotypes:    
                sample_genotype = sample_genotype.rstrip("\r\n")
                sample_ID=samples_IDs[index_position]
                sample_ID = sample_ID.rstrip("\r")

                sample_genotype = sample_genotype.split(":")
                sample_genotype = sample_genotype[0]
                sample_genotype = sample_genotype.lstrip('"')

                # Append genotypes to sample dictionary
                file_samples_genotype_dict[sample_ID].append(sample_genotype)
                    
                #Move Position Marker
                index_position+= 1
                pass

    return{'file_samples_genotype_dict': file_samples_genotype_dict, 'variant_dict': variant_dict, 'sample_list': sample_list}



def matching_samples(input_file_one_genotypes, input_file_two_genotypes, tallying_sample_dict, sample_id_key_dict, variant_one_dict):
    """
    Function accepts sample file dictionaries and tallies the number of matches and non-matches
    Between samples. The final results are outputed into various summary dictionaries

    param: input_file_one_genotypes & input_file_two_genotypes: dictionaries of all the samples and corresponding genotypes
    param: tallying_sample_dict: dictionary created of overlapping samples for tallying (empty values)
    param: sample_id_key_dict: key that matches IDs between files from various inputs
    param: variant_one_dict: a dictionary of variants with corresponding information, keys are indices based on line number

    return: tallying_sample_dict: a filled in dictionary of all the sample matching results
    return: total_results_dictionary: a summary report dictionary of all the results
    returnn: variant_one_dict: summary dictionary of all the tallying for each variant
    """

    #Creating a totals dictionary
    total_results_dictionary={'match': 0, 'no_match': 0, 'total_counts': 0, 'homo_ref': 0, 'homo_alt': 0, 'discordant': 0,
                              'no_call_DNA': 0, 'no_call_RNA': 0, 'no_call_Both': 0, 'total_no_calls': 0}

    #Start with first sample ID (x values)
    for x in (input_file_one_genotypes):

        #Some values may not exist in files being matched
        try:
            matching_sample_id = sample_id_key_dict[x]
        except KeyError:
            continue

        #Now move into the input file two data sample IDs (y values)
        for y in (input_file_two_genotypes):
            #matches sample IDs
            if matching_sample_id == y :

                #gets lists of genotypes back
                list_1=(input_file_one_genotypes[x])

                list_2=(input_file_two_genotypes[y])

                #Start Cycling Through Data
                for k in range(len(list_1)):

                    #Taking Slices of String
                    if (list_1[k] == './.' or list_1[k] == 'Filter') and (list_2[k] == './.' or list_2[k] == 'Filter'):
                        tallying_sample_dict[x]['no_call_Both']+=1
                        total_results_dictionary['no_call_Both']+=1
                        variant_one_dict[k]['no_call_Both']+=1
                        #Total No_Calls Tallying (Slightly redundant, but good for verification purposes)
                        tallying_sample_dict[x]['total_no_calls']+=1
                        total_results_dictionary['total_no_calls']+=1
                        variant_one_dict[k]['total_no_calls']+=1

                    elif (list_1[k] == './.' or list_1[k] == 'Filter'):
                        tallying_sample_dict[x]['no_call_DNA']+=1
                        total_results_dictionary['no_call_DNA']+=1
                        variant_one_dict[k]['no_call_DNA']+=1
                        #Total No_Calls Tallying (Slightly redundant, but good for verification purposes)
                        tallying_sample_dict[x]['total_no_calls']+=1
                        total_results_dictionary['total_no_calls']+=1
                        variant_one_dict[k]['total_no_calls']+=1

                    if (list_2[k] == './.' or list_2[k] == 'Filter'):
                        tallying_sample_dict[x]['no_call_RNA']+=1
                        total_results_dictionary['no_call_RNA']+=1
                        variant_one_dict[k]['no_call_RNA']+=1
                        #Total No_Calls Tallying (Slightly redundant, but good for verification purposes)
                        tallying_sample_dict[x]['total_no_calls']+=1
                        total_results_dictionary['total_no_calls']+=1
                        variant_one_dict[k]['total_no_calls']+=1
                    
                    elif list_1[k]==list_2[k]:
                        tallying_sample_dict[x]['match']+=1
                        tallying_sample_dict[x]['total_counts']+=1
                        total_results_dictionary['match']+=1
                        total_results_dictionary['total_counts']+=1
                        variant_one_dict[k]['match']+=1
                        variant_one_dict[k]['total_counts']+=1

                    else:
                        tallying_sample_dict[x]['no_match']+=1
                        tallying_sample_dict[x]['total_counts']+=1
                        total_results_dictionary['no_match']+=1
                        total_results_dictionary['total_counts']+=1
                        variant_one_dict[k]['no_match']+=1
                        variant_one_dict[k]['total_counts']+=1

                        #Breakdown of non-matching data (slicing the genotype lists)
                        # position[0] /[1] position[2]
                        #Ref Allele is only expressed in RNA-Seq
                        if list_1[k][0]==list_2[k][0] and list_1[k][0]==list_2[k][2]:
                            tallying_sample_dict[x]['homo_ref']+=1
                            total_results_dictionary['homo_ref']+=1
                            variant_one_dict[k]['homo_ref']+=1

                        #Alternative Allle is Only Expressed in RNA-Seq
                        elif list_1[k][2]==list_2[k][0] and list_1[k][2]==list_2[k][2]:
                            tallying_sample_dict[x]['homo_alt']+=1
                            total_results_dictionary['homo_alt']+=1
                            variant_one_dict[k]['homo_alt']+=1

                        #Complete Discordance
                        else:
                            tallying_sample_dict[x]['discordant']+=1
                            total_results_dictionary['discordant']+=1
                            variant_one_dict[k]['discordant']+=1
        
    #print (total_results_dictionary)
                
    total_match_counter = (total_results_dictionary['match'])
    total_non_match_counter = (total_results_dictionary['no_match'])

    #print ("The total matches are:", total_match_counter)
    #print ("The total non-matches are:", total_non_match_counter)

    #print (variant_one_dict)

    return{'tallying_sample_dict': tallying_sample_dict,
           'total_results_dictionary': total_results_dictionary,
           'variant_one_dict': variant_one_dict}


def printing_tallying_sample_report(tallying_sample_dict):
    """
    Prints the tallying results on a per sample basis for all the matching of the genotypes
    Param tallying_sample_dict: overall statitics about the matching/non-matching for each sample
    Return: NONE
    """

    tallying_sample_file = open('sample_report_file.txt', 'w')
    tallying_sample_file.write("Sample_Name\tMatch\tNo_Match\tTotal_Counts\tTotal_Match_Ratio\t\tNo_Call_DNA"
                               "\tNo_Call_RNA\tNo_Call_Both\tTotal_No_Calls\t\tNM_Homozygous_Reference"
                               "\tNM_Homozygous_Alternative\tNM_Discordant\n")

    #Looping Sample Final Results Dictionary
    for key in tallying_sample_dict:
        sample_id = key.split('.')
        tallying_sample_file.write(sample_id[0]+"\t")  
        tallying_sample_file.write(str(tallying_sample_dict[key]['match'])+"\t")
        tallying_sample_file.write(str(tallying_sample_dict[key]['no_match'])+"\t")
        tallying_sample_file.write(str(tallying_sample_dict[key]['total_counts'])+"\t")

        # Total Match Ratio for help verifying matrix grid
        match = float(tallying_sample_dict[key]['match'])
        non_match  = float(tallying_sample_dict[key]['no_match'])
        total = float(tallying_sample_dict[key]['total_counts'])
        
        # Occurs in small sample tests
        if total == 0:
            total_match_ratio = 'N.A.'
        else:
            total_match_ratio = str(round((match/total), 3))
        
        tallying_sample_file.write(total_match_ratio + "\t\t")


        tallying_sample_file.write(str(tallying_sample_dict[key]['no_call_DNA'])+"\t")
        tallying_sample_file.write(str(tallying_sample_dict[key]['no_call_RNA'])+"\t")
        tallying_sample_file.write(str(tallying_sample_dict[key]['no_call_Both'])+"\t")
        tallying_sample_file.write(str(tallying_sample_dict[key]['total_no_calls'])+"\t\t")
        tallying_sample_file.write(str(tallying_sample_dict[key]['homo_ref'])+"\t")
        tallying_sample_file.write(str(tallying_sample_dict[key]['homo_alt'])+"\t")
        tallying_sample_file.write(str(tallying_sample_dict[key]['discordant'])+"\n")
        

    tallying_sample_file.close()

    return()


def printing_tallying_variant_report(variant_tallying_results, variant_two_dict):
    """
    Prints the tallying report for the variant information for matching of the genotypes
    Param variant_tallying_results: overall statistics about the matching/non-matching for each variant
    Return: NONE
    """

    tallying_variant_file = open('variant_report_file.txt', 'w')
    tallying_variant_file.write("Chromosome\tPosition\tProbe_ID\tRS_ID\tRef_Allele\tAlt_Allele\tMatch"
                                "\tNo_Match\tTotal_Counts\tTotal_Match_Ratio\t\tNo_Call_DNA\tNo_Call_RNA\tNo_Call_Both\t"
                                "Total_No_Calls\t\tNM_Homozygous_Reference\tNM_Homozygous_Alternative\tNM_Discordant\n")

    # Looping Sample Final Results Dictionary
    for index in variant_tallying_results:
        variant_info_file_one = variant_tallying_results[index]['variant_info']
        variant_info_file_two = variant_two_dict[index]['variant_info']
        tallying_variant_file.write(variant_info_file_one[0]+"\t")
        tallying_variant_file.write(variant_info_file_one[1]+"\t")
        # printing the probe ID to file for reference
        tallying_variant_file.write(variant_info_file_one[2]+"\t")
        # printing the RS ID to file for reference
        tallying_variant_file.write(variant_info_file_two[2]+"\t")
        tallying_variant_file.write(variant_info_file_one[3]+"\t")
        tallying_variant_file.write(variant_info_file_one[4]+"\t")
        # printing actual results to the file
        tallying_variant_file.write(str(variant_tallying_results[index]['match']) + "\t")
        tallying_variant_file.write(str(variant_tallying_results[index]['no_match']) + "\t")
        tallying_variant_file.write(str(variant_tallying_results[index]['total_counts']) + "\t")
        
        # Total Match Ratio for help verifying matrix grid
        match = float(variant_tallying_results[index]['match'])
        non_match  = float(variant_tallying_results[index]['no_match'])
        total = float(variant_tallying_results[index]['total_counts'])
        
        # Occurs in small sample tests
        if total == 0:
            total_match_ratio = 'N.A.'
        else:
            total_match_ratio = str(round((match/total), 3))
        
        tallying_variant_file.write(total_match_ratio + "\t\t")
        tallying_variant_file.write(str(variant_tallying_results[index]['no_call_DNA'])+"\t")
        tallying_variant_file.write(str(variant_tallying_results[index]['no_call_RNA'])+"\t")
        tallying_variant_file.write(str(variant_tallying_results[index]['no_call_Both'])+"\t")
        tallying_variant_file.write(str(variant_tallying_results[index]['total_no_calls'])+"\t\t")
        tallying_variant_file.write(str(variant_tallying_results[index]['homo_ref'])+"\t")
        tallying_variant_file.write(str(variant_tallying_results[index]['homo_alt'])+"\t")
        tallying_variant_file.write(str(variant_tallying_results[index]['discordant'])+"\n")
        

    tallying_variant_file.close()

    return()


def create_matrix_dictionary(input_sample_genotypes, filtered_sample_list_two):
    '''
    Create a new empty dictionary for inside the loop for matching and no_matching results
    Param input_sample_genotypes: list of input_file_two's samples with genotypes
    Param filtered_sample_list_two: filtered sample list to create the dictionary
    Returns matrix_dict: matrix for tallying locally
    '''
    matrix_dict= {}
    for key in input_sample_genotypes:
        if key in filtered_sample_list_two:
            matrix_dict.update({key: {'match':0, 'no_match': 0, 'no_call': 0}})
        else:
            continue

    return (matrix_dict)


def matching_matrix(input_file_one_genotypes, input_file_two_genotypes,
                    filtered_sample_list_one, filtered_sample_list_two):
    """
    Creating a matrix from the data to see the overall behavior of variants (matching versus random chance).
    All parts of matrix creation are housed under this function for simplicity
    Param input_file_one_genotypes: genotypes from file one
    Param input_file_two_genotypes: genotypes from file two
    Param filtered_sample_list_one: list of samples to loop through, keeps order consistent in matrix files
    Param filtered_sample_list_two: list of samples to loop through, keeps order consistent in matrix files
    Return: NONE
    """

    # Open Matrix file for printing
    summary_matrix_report = open('summary_matrix_report.txt', 'w')

    #printing header to file and re-setting for newline
    summary_matrix_report.write("\t")
    for sample_id in filtered_sample_list_two:
        summary_matrix_report.write(sample_id + "\t")
    summary_matrix_report.write("\n")   

    # Built loops through lists because depending on the programs enviroment dictionaries
    # can re-order themselves and be evil
    print ("Validation of Ordering")
    print ("\n")
    print ("List from of Samples from Genotyping Panel")
    print (filtered_sample_list_one)
    print ("\n")
    print ("List from of Samples from RNA-Seq Panel")
    print (filtered_sample_list_two)
    
    #Start with first sample ID (x values)
    for sample_list_one in filtered_sample_list_one:

        # Print the ID to file
        summary_matrix_report.write(sample_list_one + "\t")

        # Creating Matrix Dictionary for the loop
        matrix_dict = create_matrix_dictionary(input_file_two_genotypes, filtered_sample_list_two)
    
        #Now move into the input file two data sample IDs (y values)
        for sample_list_two in filtered_sample_list_two:
               
                #gets lists of genotypes back
                list_1=(input_file_one_genotypes[sample_list_one])

                list_2=(input_file_two_genotypes[sample_list_two])

                #Start Cycling Through Data
                for k in range(len(list_1)):

                    #Taking Slices of String
                    if list_1[k] == './.' or list_2[k] == './.' or list_1[k] == 'Filter' or list_2[k] == 'Filter':
                        matrix_dict[sample_list_two]['no_call']+=1
                        
                    elif list_1[k]==list_2[k]:
                        matrix_dict[sample_list_two]['match']+=1

                    else:
                        matrix_dict[sample_list_two]['no_match']+=1


        #print results to file for matrix analysis
        for sample_list_two in filtered_sample_list_two:
            match = float(matrix_dict[sample_list_two]['match'])
            non_match = float(matrix_dict[sample_list_two]['no_match'])
            total = match + non_match
            if total == 0 :
                ratio_of_matching = 'N.A.'
            else:
                ratio_of_matching = str(round(match/total, 3))

            #write final answer to file
            summary_matrix_report.write(ratio_of_matching + "\t")
        #Create a newline after the loop
        summary_matrix_report.write("\n")
    
    #Closing the file
    summary_matrix_report.close()
    
    return()


def find_samples(sample_id_key_dict, sample_list, key):
    """
    Parses through the sample id key dictionary to filter the sample list for
    overlapping samples
    Param sample_id_key_dict: sample ID dicitionary used to identify samples between panels
    Param sample_list: list of samples from the panels
    Param key: trigger on if it is first or second panel being parsed
    """
    filtered_keys_list = []
    filtered_values_list = []

    # Getting keys for the 600K dataset
    if key == 'true':
        for sample in sample_list:
            try:
                # Test value in dictionary (if works good)
                matching_sample_id = sample_id_key_dict[sample]
                filtered_keys_list.append(sample)
                filtered_values_list.append(matching_sample_id)
                
            except KeyError:
                continue
    
            
    # Using filtered values to get keys
    elif key == 'false':
        for sample in sample_list:
            for value in sample_id_key_dict:
                # Find Matching Value in Sample List Dictionary
                if sample == sample_id_key_dict[value]:
                    filtered_values_list.append(sample)
                    filtered_keys_list.append(value)
                else:
                    continue
    else:
        print ("Issue with Sample Dictionary")

    #print (filtered_keys_list)
    #print (filtered_values_list)
    

    return {'filtered_keys_list': filtered_keys_list, 'filtered_values_list': filtered_values_list}

################################################################################################################

def summary_report(parameters, microarray_stats, raw_rna_seq_stats, compare_and_parse_stats, total_results_dictionary, start_time):
    """
    Prints a final summary report file of all the programs stats from the various stages
    Param parameters: parameters provided to the program for running
    Param microarray_stats: statistics about the genotyping panel
    Param raw_rna_seq_stats: statistics about the rna seq data
    Param compare_and_parse_stats: statistics about the matching of the two files
    Param total_results_dictionary: statistics about the final results of the analysis
    Param start_time: start time the program begin running
    Return: NONE
    """

    # User Input Parameters
    genotyping_file = parameters['genotyping_file']
    list_600K_snps = parameters['list_600K_snps']
    rna_seq_file = parameters['rna_seq_file']
    id_key_file = parameters['id_key_file']
    exclusion_region_length = int(parameters['exclusion_region_length'])
    quality_score_min = int(parameters['quality_score_min'])
    min_total_read_count = int(parameters['min_total_read_count'])


    raw_count = microarray_stats['raw_count']
    snp_list_count = microarray_stats['snp_list_count']
    kept_snp_count = microarray_stats['kept_snp_count']
    filtered_snp_count = microarray_stats['filtered_snp_count']
    filtered_no_calls_count = microarray_stats['filtered_no_calls_count']


    # RNA_Seq Filtering Results
    variants_in_file = raw_rna_seq_stats['variants_in_file']
    failed_GATK_SNP_filter = raw_rna_seq_stats['failed_GATK_SNP_filter']
    failed_Qual_filter = raw_rna_seq_stats['failed_Qual_filter']
    no_indel_exclusion_regions = raw_rna_seq_stats['no_indel_exclusion_regions']
    numb_SNPs_excl_Indels = raw_rna_seq_stats['numb_SNPs_excl_Indels']
    ref_allele_mutiple_forms = raw_rna_seq_stats['ref_allele_mutiple_forms']
    alt_allele_mutiple_forms = raw_rna_seq_stats['alt_allele_mutiple_forms']
    no_genotype_values = raw_rna_seq_stats['no_genotype_values']
    low_read_count = raw_rna_seq_stats['low_read_count']
    low_freq_count = raw_rna_seq_stats['low_freq_count']
    combo_filter_failure = raw_rna_seq_stats['combo_filter_failure']
    passing_variants = raw_rna_seq_stats['passing_variants']

    # Comparing VCF Panels for Common SNPs
    input_file_one = compare_and_parse_stats['input_file_one']
    input_file_two = compare_and_parse_stats['input_file_two']
    match = compare_and_parse_stats['match']
    non_matching_ref = compare_and_parse_stats['non_matching_ref']
    non_matching_alt = compare_and_parse_stats['non_matching_alt']
    both_non_matching = compare_and_parse_stats['both_non_matching']
   
 
    summary_report = open('summary_report.txt', 'w')
    summary_report.write("Running VCF Comparison Tool for RNA-Seq Data to 600K")
    summary_report.write("")
    summary_report.write("")

    # User Input Parameters Printing
    summary_report.write("User Input Parameters are the following:\n")
    summary_report.write("\n")
    summary_report.write("The genotyping file is: " + str(genotyping_file) + "\n")
    summary_report.write("The SNP file to filter genotyping file is: " + str(list_600K_snps) + "\n")
    summary_report.write("The RNA-Seq file is: " + str(rna_seq_file) + "\n")
    summary_report.write("The ID key file being used is: " + str(id_key_file) + "\n")
    summary_report.write("\n")

    # Built in for testing purposes
    testing_program_match_limit = parameters['testing_program_match_limit']
    if testing_program_match_limit == 'off':
        summary_report.write("\n")
    else:
        summary_report.write("Program Testing Mode is Turned On\n")
        summary_report.write("The number of variant matches is: " + str(testing_program_match_limit))
        summary_report.write("\n")

    summary_report.write("##########################################################")
    summary_report.write("\n")
    summary_report.write("Filtering Parameters for 600K Genotyping Data\n")
    summary_report.write("\n")
    summary_report.write("Number of SNPs in 600K Genotyping Data: " + str(raw_count) + "\n")
    summary_report.write("Number of SNPs in SNP list: " + str(snp_list_count) + "\n")
    summary_report.write("Number of SNPs filtered by user cutoff : " + str(filtered_snp_count) + "\n")
    summary_report.write("Number of SNPs filtered for no calls : " + str(filtered_no_calls_count) + "\n")
    summary_report.write("Number of SNPs left after all filters : " + str(kept_snp_count) + "\n")
    summary_report.write("\n")
    summary_report.write("##########################################################")
    summary_report.write("\n")
    summary_report.write("Filtering Parameters for RNA-Seq Data\n")
    summary_report.write("\n")
    summary_report.write("Exclusion length for INDELs: " + str(exclusion_region_length) + "\n")
    summary_report.write("Minimum Quality Score (Phred) for Variants: " + str(quality_score_min) + "\n") 
    summary_report.write("Minimum Coverage for a Sample: " + str(min_total_read_count) + "\n")
    summary_report.write("\n")
    summary_report.write("Pipeline first filters RNA-Seq data then begins actual analysis of files")
    summary_report.write("\n")
    summary_report.write("SNPs in RNA-Seq File: " + str(variants_in_file) + "\n")
    summary_report.write("\n")
    summary_report.write("Variant Level Filtering\n")
    summary_report.write("RNA-Seq SNPs that failed GATK filters: " + str(failed_GATK_SNP_filter) + "\n")
    summary_report.write("RNA-Seq SNPs that failed Quality Score Filter): "+str(failed_Qual_filter) + "\n")
    summary_report.write("\n")
    summary_report.write("INDEL regions identified in VCF file: " + str(no_indel_exclusion_regions) + "\n")
    summary_report.write("RNA-Seq SNPs filtered by INDEL filtering: " + str(numb_SNPs_excl_Indels) +"\n")
    summary_report.write("\n")
    summary_report.write("RNA-Seq SNPs filtered with too many Ref Alleles: "+ str(ref_allele_mutiple_forms) + "\n")
    summary_report.write("RNA-Seq SNPs filtered with too many Alt Alleles: " + str(alt_allele_mutiple_forms) + "\n")
    summary_report.write("\n")
    summary_report.write("Sample Level Filtering Causing Variant to be Filtered\n")
    summary_report.write("RNA_Seq SNPs filtered because no genotype values: " + str(no_genotype_values) + "\n")
    summary_report.write("RNA_Seq SNPs filtered because low read count: " + str(low_read_count) + "\n")
    summary_report.write("RNA_Seq SNPs filtered because low freq count (one allele count <1% total counts): " + str(low_freq_count) + "\n")
    summary_report.write("RNA_Seq SNPs filtered because of combination of sample level filters: " + str(combo_filter_failure) + "\n")
    summary_report.write("\n")
    summary_report.write("Total number of Variants that passed: " + str(passing_variants) + "\n")
    summary_report.write("\n")
    summary_report.write("##########################################################")
    summary_report.write("\n")
    summary_report.write("Comparing VCF Files for Common Variants\n")
    summary_report.write("\n")
    summary_report.write("Total number of Variants In Genotyping Panel: " + str(input_file_one) + "\n")
    summary_report.write("Total number of Variants In Filtered RNA-Seq Panel: " + str(input_file_two) + "\n")
    summary_report.write("\n")
    summary_report.write("Total number of Matching Variants (Ref and Alt Allele Match): " + str(match) + "\n")
    summary_report.write("\n")
    summary_report.write("Discordant Matching Variants (alleles do not match)\n")
    summary_report.write("Were not further analyzed\n")
    summary_report.write("Total number of Matching Variants but Ref allele Was Discordant: " + str(non_matching_ref) + "\n")
    summary_report.write("Total number of Matching Variants but Alt allele Was Discordant: " + str(non_matching_alt) + "\n")
    summary_report.write("Total number of Matching Variants but Both Alleles Were Discordant: " + str(both_non_matching) + "\n")
    summary_report.write("\n")
    summary_report.write("\n") 
    summary_report.write("##########################################################")
    summary_report.write("\n")
    summary_report.write("Final Results of Analysis")
    summary_report.write("\n")
    summary_report.write("\n")
    summary_report.write("The total number of matches was: " + str(total_results_dictionary['match'])+"\n")
    summary_report.write("The total number of non-matches was: " + str(total_results_dictionary['no_match'])+"\n")
    summary_report.write("The total number of measurable counts was: " + str(total_results_dictionary['total_counts'])+"\n")

    #Final concordance calculation
    match_total = float(total_results_dictionary['match'])
    total_counts = float(total_results_dictionary['total_counts'])
    overall_concordance = str(round((match_total / total_counts * 100), 2 ))
    

    summary_report.write("Overall Concordance was: " + overall_concordance +"% \n")
    summary_report.write("\n")
    summary_report.write("The total number of no-calls DNA was: " + str(total_results_dictionary['no_call_DNA'])+"\n")
    summary_report.write("The total number of no-calls RNA was: " + str(total_results_dictionary['no_call_RNA'])+"\n")
    summary_report.write("The total number of no-calls Both was: " + str(total_results_dictionary['no_call_Both'])+"\n")
    summary_report.write("The total number of no-calls was: " + str(total_results_dictionary['total_no_calls'])+"\n")
    summary_report.write("\n")
    summary_report.write("Breakdown of Non-Matches")
    summary_report.write("\n")
    summary_report.write("Genotype is homozygous reference allele: " + str(total_results_dictionary['homo_ref'])+"\n")
    summary_report.write("Genotype is homozygous alternative allele: " + str(total_results_dictionary['homo_alt'])+"\n")
    summary_report.write("Genotype is discordance and represents a new allele: " + str(total_results_dictionary['discordant'])+"\n")


    #Started time in the program to see how long runs take
    total_time = time.clock() - start_time
    summary_report.write("\n")
    summary_report.write("\n")
    summary_report.write("Program ran for a total of " + str(round(total_time, 2))+" seconds")
    
    summary_report.close()
    



def main():
    
    print ("VCF Comparison Program Has Started Running")
    print ("")
    print ("")

    #Setting up program
    start_time=time.clock()
    home_directory = os.getcwd()

    # Opens the parameter file to get all the required inputs for the rest of the code
    parameters=read_parameter_file() 

    # Program parameters
    genotyping_file = parameters['genotyping_file']
    print ("The genotyping file is: ", genotyping_file)
    print ("")
    list_600K_snps = parameters['list_600K_snps']
    print ("The SNPs to filter genotyping file is: ", list_600K_snps)
    print ("")
    rna_seq_file = parameters['rna_seq_file']
    print ("The RNA-Seq file is: ", rna_seq_file)
    print ("")
    id_key_file = parameters['id_key_file']
    print ("The ID Key file is: ", id_key_file)
    print ("")
    exclusion_region_length = int(parameters['exclusion_region_length'])
    print ("The INDEL exclusion length is: ", exclusion_region_length)
    print ("")
    quality_score_min = int(parameters['quality_score_min'])
    print ("The minimum quality score for variants is: ", quality_score_min)
    print ("")
    min_total_read_count = int(parameters['min_total_read_count'])
    print ("The total minimum read count for a sample is: ", min_total_read_count)
    print ("")
    testing_program_match_limit = parameters['testing_program_match_limit']
    if testing_program_match_limit == 'off':
        print ("")
        pass
    else:
        print ("Program Testing Mode is Turned On")
        print ("The number of variant matches is: ", testing_program_match_limit)
        print ("")


    #Making a directory to put all important, but non-essential results in
    # See if directory exists otherwise make it
    verdict = os.path.exists('Log_Directory')
    if str(verdict) == 'False':
        os.makedirs('Log_Directory')
    else:
        print("Log Directory already exists")
        print("")

    print ("Filtering the 600K Chip Data -Keeping SNPs Only Found in the Provided SNP List")
    # Filter 600K Genotyping panel based on MAF and No Calls
    filter_600K_results = filter_600K_Chip(genotyping_file, list_600K_snps)
    filtered_600K_file_name = filter_600K_results['filtered_600K_file_name']
    microarray_stats = filter_600K_results['microarray_stats']

    print("")
    print ("Filtering the RNA-Seq Data Based on User Parameters")
    # Filter RNA_Seq_Data File before begining matching
    filtering_results = filter_RNA_Seq_Data(rna_seq_file, exclusion_region_length, quality_score_min, min_total_read_count)
    filtered_rna_seq_file_name = filtering_results['filtered_rna_seq_file_name']
    raw_rna_seq_stats = filtering_results['raw_rna_seq_stats']
    

    print("")
    print ("Comparing VCF files for Common Variants")
    # Compare Files to Get Common SNPs
    # Input_file_One is always 600K data (first file listed)
    # Input_file_Two is the Filtered RNA-Seq Data
    compare_and_parse_results = compare_and_parse(filtered_600K_file_name, filtered_rna_seq_file_name, testing_program_match_limit)
    common_snps_file_one_name = compare_and_parse_results['common_snps_file_one_name']
    common_snps_file_two_name = compare_and_parse_results['common_snps_file_two_name']
    compare_and_parse_stats = compare_and_parse_results['compare_and_parse_stats']
    

    #Comparing Final Common SNPs Between Files
    key_dict_file = id_key_file
    
    key_dictionaries = parse_key_file(key_dict_file)
    sample_id_key_dict = key_dictionaries['sample_id_key_dict']
    tallying_sample_dict = key_dictionaries['tallying_sample_dict']

    # Create Dictionaries of Samples (keys) and Genotypes (values)
    # Input File One
    input_file_one_data = samples_data (common_snps_file_one_name)
    input_file_one_genotypes = input_file_one_data['file_samples_genotype_dict']
    variant_one_dict = input_file_one_data['variant_dict']
    sample_list_one = input_file_one_data['sample_list']
    
    # Input File Two
    input_file_two_data = samples_data (common_snps_file_two_name)
    input_file_two_genotypes = input_file_two_data['file_samples_genotype_dict']
    variant_two_dict = input_file_two_data['variant_dict']
    sample_list_two = input_file_two_data['sample_list']

    print("")
    print ("Analyzing VCF files for Genotypes to Identify Overall Concordance")
    # Passing Only One Variant Dictionary (both files their the same)    
    matching_results = matching_samples(input_file_one_genotypes, input_file_two_genotypes,
                                        tallying_sample_dict, sample_id_key_dict, variant_one_dict)

    tallying_sample_dict = matching_results['tallying_sample_dict']
    total_results_dictionary = matching_results['total_results_dictionary']
    variant_tallying_results = matching_results['variant_one_dict']

    # Print Results from Analysis
    printing_tallying_sample_report(tallying_sample_dict)
    printing_tallying_variant_report(variant_tallying_results, variant_two_dict)

    # Get the filtered sample lists for usage with matrix output
    #Note Key can be turned to false to turn on for search for values instead of keys
    key = 'true'
    filtered_sample_list = find_samples(sample_id_key_dict, sample_list_one, key)
    #retrun lists of keys and values that are guranteed to match
    filtered_key_list = filtered_sample_list['filtered_keys_list']
    filtered_values_list = filtered_sample_list['filtered_values_list']

    print("")
    print ("Creating a Sample Matrix to See Overall Behavior Among Samples")
    # Setup a matrix of results to see overall matching behavior if radomized samples
    matching_matrix(input_file_one_genotypes, input_file_two_genotypes,
                    filtered_key_list, filtered_values_list)
    
    # Summary Stats Report File
    summary_report(parameters, microarray_stats, raw_rna_seq_stats, compare_and_parse_stats, total_results_dictionary, start_time)
       

    print ("VCF Comparison Program Has Finished Running")

if __name__ == "__main__":
    
    main()


