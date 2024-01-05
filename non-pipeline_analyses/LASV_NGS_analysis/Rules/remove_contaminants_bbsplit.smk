
# Description:
# Remove reads that map to contaminant genomes using bbsplit. The contaminant
# reads are filtered using the contaminant genomes listed in the config.yml file.

# Author:
# Caleb Carr


rule bbsplit_remove_sample_contaminant_reads_PE:
    """ This rule runs BBsplit using the contaminant genomes as reference and 
        keeps all unmatched reads. This rule uses the previously created index 
        contaminant genomes. """
    input:
        # Index log created to signal that the contaminant
        # merged index has been created. Function in common.smk
        # retrieves the index log file for the contaminant genomes.
        get_contaminant_index_log_files(),
        # Function in common.smk to retreive the quality filtered fastq files
        read_1 = get_quality_filtered_fastq_files_PE()[0],
        read_2 = get_quality_filtered_fastq_files_PE()[1],
    output:
        out_1 = "Data/Samples/Cleaned_Reads/{accession}/{accession}_1_BBSPLIT_cleaned.fastq",
        out_2 = "Data/Samples/Cleaned_Reads/{accession}/{accession}_2_BBSPLIT_cleaned.fastq",
    params:
        # Directory where the merged reference files are stored
        ref_dir = "Data/Contaminant_Genomes/BBSplit_Index/"
    conda:
        "../Conda_Envs/read_processing_env.yml"
    shell:
        # The 'in1' and 'in2' inputs are the recently qfiltered reads while the
        # 'outu1' and 'outu2' are the two output files for only unmapped reads.
        # The 'path' flag denotes the directories where the previously created ref 
        # directories where the indexed contaminant genomes are stored. The 'minid=0.95'
        # flag refers to the approximate alignment identity to look for (higher percentages
        # means faster but less sensitive). The 'maxindel=3' flag means don't look for 
        # indels longer than 3 (smaller settings means faster performance). The 'bandwidthratio' 
        # flag set above zero restricts the alignment band to this ratio of read length 
        # (larger values are less accurate but faster). The 'bandwidth' flag sets the fraction 
        # of the read length (greater than 0 is faster but less accurate). The 'quickmatch' flag 
        # means that cigar strings are generated quicker. The 'fast' flag sets other parameters 
        # to run quicker but this is bad for RNA-seq (okay for contaminant removal). The 
        # 'minhits=2' flag is the minimum number of seed hits required for candidate sites 
        # (larger values means faster performance). The 'printunmappedcount' flag means that 
        # the number of unmapped reads is printed. The 'usemodulo' is required to reduce the 
        # amount of RAM needed because without, the human genome requires 24G compared to 
        # 12G with 'usemodulo'.
        "bbsplit.sh "
        "in1={input.read_1} in2={input.read_2} "
        "path={params.ref_dir} "
        "outu1={output.out_1} outu2={output.out_2} "
        "minid=0.9 maxindel=3 bandwidthratio=0.16 "
        "bandwidth=12 quickmatch fast minhits=2 "
        "printunmappedcount "
        "usemodulo"