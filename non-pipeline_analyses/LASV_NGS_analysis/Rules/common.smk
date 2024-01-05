
# Description:
# All common functions that are used across each snakemake rule
# file are stored here. The first section controls the rules used
# based on the settings in the config.yml file. The second section
# contains basic functions to retrieve file names for different 
# rules.

# Author:
# Caleb Carr


###########################################---START TOOL SETTINGS---###########################################
def get_adapter_trimmed_fastq_files_PE():
    """ Based on the config file configuration, return the
        corresponding file names for the trimmed reads in
        a list. """
    if config['Adapter_Trimming']['BBDuk'] == 'yes_bbduk':
        read_1 = "Data/Samples/Cleaned_Reads/{accession}/{accession}_1_BBDUK_trimmed.fastq"
        read_2 = "Data/Samples/Cleaned_Reads/{accession}/{accession}_2_BBDUK_trimmed.fastq"
        return [read_1, read_2]
    elif config['Adapter_Trimming']['Fastp'] == 'yes_fastp':
        read_1 = "Data/Samples/Cleaned_Reads/{accession}/{accession}_1_FASTP_cleaned.fastq"
        read_2 = "Data/Samples/Cleaned_Reads/{accession}/{accession}_2_FASTP_cleaned.fastq"
        return [read_1, read_2]
    elif config['Adapter_Trimming']['Cutadapt'] == 'yes_cutadapt':
        read_1 = "Data/Samples/Cleaned_Reads/{accession}/{accession}_1_CUTADAPT_cleaned.fastq"
        read_2 = "Data/Samples/Cleaned_Reads/{accession}/{accession}_2_CUTADAPT_cleaned.fastq"
        return [read_1, read_2]
    else:
        print('ERROR in config.yml file: Please write only one adapter setting.')
        return 'ERROR in config.yml file: Please write only one adapter setting.'


def get_quality_filtered_fastq_files_PE():
    """ Based on the config file configuration, return the 
        corresponding file names for the quality filtered
        read file names in a list. """
    if config['Quality_Filtering']['BBDuk'] == 'yes_bbduk':
        read_1 = "Data/Samples/Cleaned_Reads/{accession}/{accession}_1_BBDUK_qfilter_cleaned.fastq"
        read_2 = "Data/Samples/Cleaned_Reads/{accession}/{accession}_2_BBDUK_qfilter_cleaned.fastq"
        return [read_1, read_2]
    else:
        print('ERROR in config.yml file: Please write only one quality filtering setting.')
        return 'ERROR in config.yml file: Please write only one quality filtering setting.'


def get_contaminant_index_log_files():
    """ Based on the config file configuration, return the
        corresponding index log files for the contaminants. """
    if config['Contaminant_Removal']['BBSplit'] == 'yes_bbsplit':
        return "Data/Contaminant_Genomes/BBSplit_Index/Merged_Index_BBSplit_contaminant_index_log.txt"
    else:
        print('ERROR in config.yml file: Please write only one contaminant removal setting.')
        return 'ERROR in config.yml file: Please write only one contaminant removal setting.'


def get_cleaned_fastq_files_PE():
    """ Based on the config file configuration, return the
        corresponding file names for the reads with the
        contaminant reads removed in a list. """
    if config['Contaminant_Removal']['BBSplit'] == 'yes_bbsplit':
        read_1 = "Data/Samples/Cleaned_Reads/{accession}/{accession}_1_BBSPLIT_cleaned.fastq"
        read_2 = "Data/Samples/Cleaned_Reads/{accession}/{accession}_2_BBSPLIT_cleaned.fastq"
        return [read_1, read_2]
    else:
        print('ERROR in config.yml file: Please write only one contaminant removal setting.')
        return 'ERROR in config.yml file: Please write only one contaminant removal setting.'


def get_ref_report():
    """ Based on the config file configuration, return the
        corresponding file name for best matched reference report. """
    if config['Reference_Checking']['BBSplit'] == 'yes_bbsplit':
        return "Data/Samples/Cleaned_Reads/{accession}/{accession}_ref_report_BBSPLIT.txt"
    else:
        print('ERROR in config.yml file: Please write only one reference checking setting.')
        return 'ERROR in config.yml file: Please write only one reference checking setting.'


def get_ref_index_log_files():
    """ Based on the config file configuration, return the
        corresponding index log files for the references. """
    if config['Alignment_Mapping']['BBMap'] == 'yes_bbmap':
        return expand("Data/Reference_Genomes/{ref_name}/BBMap_Index/{ref_name}_index_log_BBMAP.txt", ref_name=config['References'])
    elif config['Alignment_Mapping']['BWA'] == 'yes_bwa':
        return expand("Data/Reference_Genomes/{ref_name}/BWA_Index/{ref_name}_index_log_BWA.txt", ref_name=config['References'])
    elif config['Alignment_Mapping']['STAR'] == 'yes_star':
        return expand("Data/Reference_Genomes/{ref_name}/STAR_Index/{ref_name}_index_log_STAR.txt", ref_name=config['References']) 
    elif config['Alignment_Mapping']['HISAT2'] == 'yes_hisat2':
        return expand("Data/Reference_Genomes/{ref_name}/HISAT2_Index/{ref_name}_index_log_HISAT2.txt", ref_name=config['References'])
    else:
        print('ERROR in config.yml file: Please write only one alignment setting.')
        return 'ERROR in config.yml file: Please write only one alignment setting.'


def get_mapped_sam_files():
    """ Based on the config file configuration, return the
        corresponding file name for the aligned sam file. """
    if config['Alignment_Mapping']['BBMap'] == 'yes_bbmap':
        return "Data/Samples/Mapped_Reads/{accession}/{accession}_BBMAP.sam"
    elif config['Alignment_Mapping']['BWA'] == 'yes_bwa':
        return "Data/Samples/Mapped_Reads/{accession}/{accession}_BWA.sam"
    elif config['Alignment_Mapping']['STAR'] == 'yes_star':
        return "Data/Samples/Mapped_Reads/{accession}/{accession}_STAR.sam"
    elif config['Alignment_Mapping']['HISAT2'] == 'yes_hisat2':
        return "Data/Samples/Mapped_Reads/{accession}/{accession}_HISAT2.sam"
    else:
        print('ERROR in config.yml file: Please write only one alignment setting.')
        return 'ERROR in config.yml file: Please write only one alignment setting.'


def get_sorted_bam():
    """ Based on the config file configuration, return the corresponding
        file name for the sorted bam file. """
    if config['Deduplication']['No_Deduplication'] == 'no_dedup':
        return "Data/Samples/Mapped_Reads/{accession}/{accession}_sorted.bam"
    elif config['Deduplication']['Samtools'] == 'yes_samtools':
        return "Data/Samples/Mapped_Reads/{accession}/{accession}_sorted_dedup.bam"
    else:
        print('ERROR in config.yml file: Please write only one deduplication setting.')
        return 'ERROR in config.yml file: Please write only one deduplication setting.'
###########################################---END TOOL SETTINGS---###########################################


###########################################---START FILE NAME RETRIEVAL FUNCTIONS---#########################
def get_fastq_files_PE():
    """ Return a list containing the two read fastq files. """
    read_1 = "Data/Samples/Raw_Reads/{accession}/{accession}_1.fastq"
    read_2 = "Data/Samples/Raw_Reads/{accession}/{accession}_2.fastq"
    return [read_1, read_2]


def get_list_of_sample_accessions():
    """ Return a list containing all the sample accession names. """
    return pd.read_csv(config['Samples']['File'])['Run'].tolist()


def get_list_of_reference_genomes():
    """ Produces the full list of reference genomes in fasta format. """
    return expand("Data/Reference_Genomes/{name}/{name}.fasta", name=config['References'])


def get_list_of_contaminant_genomes():
    """ Produces the full list of contaminant genomes in fasta format. """
    return expand("Data/Contaminant_Genomes/{name}.fasta", name=config['Sample_Contaminants'])


def get_most_likely_ref_index(accession_number):
    """ Based on the accession number passed to the function, the path 
        to the best matched reference genome is returned. """
    if config['Reference_Checking']['BBSplit'] == 'yes_bbsplit':
        path_to_report = f"Data/Samples/Cleaned_Reads/{accession_number}/{accession_number}_ref_report_BBSPLIT.txt"
        # Check if report exists because snakemake executes params functions in the DAG phase
        # which is before the files habe been created. This 'if' statement prevents an error 
        # from being thrown.
        if os.path.exists(path_to_report):
            reference_mapped_rates = pd.read_csv(path_to_report, sep='\t')['#name'].tolist()
            # Check if reference_mapped_rates is not empty because if no reads were mapped
            # to any genome then this will be empty and cause an error
            if len(reference_mapped_rates) == 0:
                # Use default reference genome
                most_likely_ref_genome = config['Default_Reference']
            else:
                # First reference genome listed is the best matched reference genome.
                most_likely_ref_genome = reference_mapped_rates[0]
            # Print statement doesn't allow for a DAG diagram to be created
            # print("Report was created and best match genome was extracted.")
            if config['Alignment_Mapping']['BBMap'] == 'yes_bbmap':
                return f"Data/Reference_Genomes/{most_likely_ref_genome}/BBMap_Index/"
            elif config['Alignment_Mapping']['BWA'] == 'yes_bwa': 
                return f"Data/Reference_Genomes/{most_likely_ref_genome}/BWA_Index/{most_likely_ref_genome}.fasta"
            elif config['Alignment_Mapping']['STAR'] == 'yes_star':
                return f"Data/Reference_Genomes/{most_likely_ref_genome}/STAR_Index/"
            elif config['Alignment_Mapping']['HISAT2'] == 'yes_hisat2':
                return f"Data/Reference_Genomes/{most_likely_ref_genome}/HISAT2_Index/{most_likely_ref_genome}"
            else:
                print('ERROR in config.yml file: Please write only one alignment setting.')
                return 'ERROR in config.yml file: Please write only one alignment setting.'
        else:
            # Print statement doesn't allow for a DAG diagram to be created
            # print("Still in DAG phase and file has not been created.")
            return ""
    else:
        print('ERROR in config.yml file: Please write only one reference checking setting.')
        return 'ERROR in config.yml file: Please write only one reference checking setting.'


def get_most_likely_ref_genome(accession_number):
    """ Based on the accession number passed to the function, the file 
        to the best matched reference genome is returned. """
    if config['Reference_Checking']['BBSplit'] == 'yes_bbsplit':
        path_to_report = f"Data/Samples/Cleaned_Reads/{accession_number}/{accession_number}_ref_report_BBSPLIT.txt"
        # Check if report exists because snakemake executes params functions in the DAG phase
        # which is before the files habe been created. This 'if' statement prevents an error 
        # from being thrown.
        if os.path.exists(path_to_report):
            reference_mapped_rates = pd.read_csv(path_to_report, sep='\t')['#name'].tolist()
            # Check if reference_mapped_rates is not empty because if no reads were mapped
            # to any genome then this will be empty and cause an error
            if len(reference_mapped_rates) == 0:
                # Use default reference genome
                default_ref_genome = config['Default_Reference']
                return f"Data/Reference_Genomes/{default_ref_genome}/{default_ref_genome}.fasta"
            # First reference genome listed is the best matched reference genome.
            most_likely_ref_genome = reference_mapped_rates[0]
            # Print statement doesn't allow for a DAG diagram to be created
            # print("Report was created and best match genome was extracted.")
            return f"Data/Reference_Genomes/{most_likely_ref_genome}/{most_likely_ref_genome}.fasta"
        else:
            # Print statement doesn't allow for a DAG diagram to be created
            # print("Still in DAG phase and file has not been created.")
            return ""
    else:
        print('ERROR in config.yml file: Please write only one reference checking setting.')
        return 'ERROR in config.yml file: Please write only one reference checking setting.'

###########################################---END FILE NAME RETRIEVAL FUNCTIONS---############################


