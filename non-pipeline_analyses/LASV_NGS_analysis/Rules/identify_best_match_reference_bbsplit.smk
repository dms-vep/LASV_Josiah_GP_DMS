
# Description:
# Checks the best matched reference genome using bbsplit.

# Author:
# Caleb Carr



rule bbsplit_check_strain_type:
    """ This rule runs BBSplit with the reference genomes inputed as well as each 
        data sample to get the best mapping rate for each data sample to determine 
        the correct reference genome. """
    input:
        # Temporary index log created to signal that the reference
        # merged index has been created
        "Data/Reference_Genomes/Merged_Ref_BBSPLIT/merged_ref_index_BBSPLIT.txt",
        # Function in common.smk to retrieve all cleaned fastq files
        read_1 = get_cleaned_fastq_files_PE()[0],
        read_2 = get_cleaned_fastq_files_PE()[1],
    output:
        "Data/Samples/Cleaned_Reads/{accession}/{accession}_ref_report_BBSPLIT.txt"
    params:
        # Directory where the merged reference files are placed
        ref_dir = "Data/Reference_Genomes/Merged_Ref_BBSPLIT/"
    conda:
        "../Conda_Envs/genome_processing_env.yml"
    shell:
        # The 'in1' and 'in2' flags refers to the paired end read input files for each 
        # sample. The 'path' flag signals where to put the merged index for all genomes. 
        # BBSplit must have all input genomes merged into one index. The 'samplerate' flag
        # refers to the fraction of reads from each sample that are randomly chosen 
        # to align to the reference genomes. In this case, 5% of reads are aligned.
        # The 'ambiguous2=all' flag means that any ambiguous reads that map to more 
        # than one genome are mapped to all genomes. The 'refstats' flag is to collect
        # information about the mapping rates for reference genomes in an output file. 
        "bbsplit.sh "
        "in1={input.read_1} in2={input.read_2} "
        "path={params.ref_dir} "
        "samplerate=0.005 "
        "ambiguous2=all "
        "refstats={output}"