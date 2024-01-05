
# Description:
# Remove adapters using bbduk.

# Author:
# Caleb Carr


rule bbduk_adapter_trim_PE:
    """ This rule runs bbduk which is a Decontamination Using Kmers preprocessing 
        tool for raw reads. This allows for adapters to be trimmed based on common 
        adapter sequences present in the adapter reference file. """
    input:
        # Function in common.smk to retrieve fastq files
        read_1 = get_fastq_files_PE()[0],
        read_2 = get_fastq_files_PE()[1],
    output:
        out_1 = temp("Data/Samples/Cleaned_Reads/{accession}/{accession}_1_BBDUK_trimmed.fastq"),
        out_2 = temp("Data/Samples/Cleaned_Reads/{accession}/{accession}_2_BBDUK_trimmed.fastq"),
    params:
        # A list of potential adapter sequences is referenced
        # in the config file
        adapters = config['Default_Adaptors']
    conda:
        "../Conda_Envs/read_processing_env.yml"
    shell:
        # in1 and in2 are the paired read input files while out1 and out2 are the 
        # processed output files. The 'ref' flag refers to a file of potential adapter 
        # sequences to use to trim. The'ktrim=r' flag is for 3' adapters which is the 
        # most common for illumina sequencing. The 'k=23' flag means that the adapter
        # sequence is matched for 23 bp while the 'mink=11' flag means the smallest match can be
        # 11 bp. The 'hdist=1' means a hamming distance of 1 so at most one mismatch is allowed
        # for a matched adapter sequence. The 'tbo' flag detects adapter sequences based on overlap
        # while the 'tpe' flag cuts both reads to same length if only one pair was cut originally.
        "bbduk.sh "
        "in1={input.read_1} in2={input.read_2} ref={params.adapters} "
        "out1={output.out_1} out2={output.out_2} "
        "ktrim=r k=23 mink=11 hdist=1 tpe tbo"