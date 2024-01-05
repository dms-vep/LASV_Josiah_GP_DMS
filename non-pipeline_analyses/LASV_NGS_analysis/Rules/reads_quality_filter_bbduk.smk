
# Description:
# Filter reads based on quality using bbduk. This also removes
# phix contaminants as well.

# Author:
# Caleb Carr


rule bbduk_remove_contaminant_phix_reads_PE:
    """ This rule runs BBDuk to remove contaminant reads associated with the phix 
        control library. The Illumina phix reference genome is listed in the config 
        file which was the same reference listed in the bbmap suite. """
    input:
        # Function in common.smk to retrieve the adapter trimmed fastq files
        read_1 = get_adapter_trimmed_fastq_files_PE()[0],
        read_2 = get_adapter_trimmed_fastq_files_PE()[1],
    output:
        out_1 = temp("Data/Samples/Cleaned_Reads/{accession}/{accession}_1_BBDUK_phix_cleaned.fastq"),
        out_2 = temp("Data/Samples/Cleaned_Reads/{accession}/{accession}_2_BBDUK_phix_cleaned.fastq"),
    params:
        # The phix reference is linked to in the config file
        phix = config['Contaminants']['Phix']
    conda:
        "../Conda_Envs/read_processing_env.yml"
    shell:
        # The input reads are the recently trimmed reads while the 
        # output reads are the phix decontaminated reads. The flags
        # 'k=31' and 'hdist=1' mean all reads with a 31-kmer match 
        # to phix reference are removed with a hamming distance of 1.
        "bbduk.sh "
        "in1={input.read_1} in2={input.read_2} ref={params.phix} "
        "out1={output.out_1} out2={output.out_2} "
        "k=31 hdist=1"


rule bbduk_quality_filter_PE:
    """ This rule runs BBDuk to filter and trim the reads based on quality 
        on either end, length of reads, and presence of a poly A tail. """
    input:
        read_1 = "Data/Samples/Cleaned_Reads/{accession}/{accession}_1_BBDUK_phix_cleaned.fastq",
        read_2 = "Data/Samples/Cleaned_Reads/{accession}/{accession}_2_BBDUK_phix_cleaned.fastq",
    output:
        out_1 = temp("Data/Samples/Cleaned_Reads/{accession}/{accession}_1_BBDUK_qfilter_cleaned.fastq"),
        out_2 = temp("Data/Samples/Cleaned_Reads/{accession}/{accession}_2_BBDUK_qfilter_cleaned.fastq"),
    conda:
        "../Conda_Envs/read_processing_env.yml"
    shell:
        # The 'qtrim=rl' flag means that both ends (rl) will be trimmed based on a quality threshold 
        # specified by the 'trimq' flag. In this case, anything below a phred score of 25 will be 
        # trimmed. The 'trimpolya=10' flag means polyA tails of at least 10 bp are cut on both sides
        # and the 'minlength=25' flag means that reads of less than 25 bp are discarded. 
        "bbduk.sh "
        "in1={input.read_1} in2={input.read_2} "
        "out1={output.out_1} out2={output.out_2} "
        "qtrim=rl trimq=25 "
        "trimpolya=10 "
        "minlength=25"