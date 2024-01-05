
# Description:
# Remove adapters using cutadapt. Cutadapt takes a very long time, 
# especially when a long list of potential adapters is listed.

# Author:
# Caleb Carr


rule run_cutadapt_PE:
    """ This rule runs cutadapt which is a program used to remove adapter sequences as well as
        trim and filter the data for better quality. HOWEVER, cutadapt is SUBSTANTIALLY
        SLOWER than fastp and bbduk so be prepared to wait... """
    input:
        # Function in common.smk to retrieve fastq files
        read_1 = get_fastq_files_PE()[0],
        read_2 = get_fastq_files_PE()[1],
    output:
        out_1 = temp("Data/Samples/Cleaned_Reads/{accession}/{accession}_1_CUTADAPT_cleaned.fastq"),
        out_2 = temp("Data/Samples/Cleaned_Reads/{accession}/{accession}_2_CUTADAPT_cleaned.fastq"),
    params:
        # A list of potential adapter sequences is referenced
        # in the config file
        adapters = config['Default_Adaptors'],
        # Cutadapt does not have a built in polyA tail trim
        # so a polyA sequence must be provided
        polyA = "AAAAAAAAAA"
    conda:
        "../Conda_Envs/read_processing_env.yml"
    shell:
        # The '-a' flag is for the first read adapter list while the '-A' is for the second
        # read adapter list. The 'file:' notation allows a file of adapters to be passed to 
        # the tool and the params.polyA tail refers to a polyA tail sequence to trim. Both are 
        # passed to the tool as regular adapter sequences. The '-o' and '-p' flags signify 
        # the output files for the cleaned reads. The '--quality-cutoff 25,25' flag means 
        # that both the 5' and 3' ends are trimmed based on the quality thresholds specified. 
        # In this case, ends with phred scores less than 25 are trimmed. The '--minimum-length 25'
        # flag means that all reads less than 25 bp are discarded. 
        "cutadapt "
        "-a file:{params.adapters} -a {params.polyA} -A file:{params.adapters} -A {params.polyA} "
        "-o {output.out_1} -p {output.out_2} "
        "{input.read_1} {input.read_2} "
        "--quality-cutoff 25,25 "
        "--minimum-length 25"
