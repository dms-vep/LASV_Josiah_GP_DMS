
# Description:
# Remove adapters using fastp. 

# Author:
# Caleb Carr


rule run_fastp_PE:
    """ This rule runs fastp which is an all in one preprocessing tool, so it will 
        trim adapters, filter based on quality thresholds, trim polyA tails, etc.
        This rule is setup for paired end reads only because it can automatically detect
        adapters using the paired ends overlap. """
    input: 
        # Function in common.smk to retrieve fastq files
        read_1 = get_fastq_files_PE()[0],
        read_2 = get_fastq_files_PE()[1],
    output: 
        out_1 = temp("Data/Samples/Cleaned_Reads/{accession}/{accession}_1_FASTP_cleaned.fastq"),
        out_2 = temp("Data/Samples/Cleaned_Reads/{accession}/{accession}_2_FASTP_cleaned.fastq"),
        html_report = "Results/Quality_Checking/{accession}/Fastp/{accession}_fastp.html",
        json_report = "Results/Quality_Checking/{accession}/Fastp/{accession}_fastp.json",
    conda:
        "../Conda_Envs/read_processing_env.yml"
    shell: 
        # in1 and in2 are the paired end read files while out1 and out2 are the paired end 
        # output files. The html and json reports are specified by the files given. The '--qualified_quality_phred' 
        # flag is a quality threshold for each read which is set at a phred score of 25. By default 40% of the read
        # needs to be a phred score of 25 or greater to not be filter out. The '--length_required' flag is a 
        # filter based on size of each read which is set at 25. This means all reads less than 25 bp are filtered out.
        # The '--trim_poly_x' flag is a trimmer for polyA tails. The default value is set at 10 which means 10 A bp must 
        # be detected to be trimmed. This is important because viruses can have polyA tails. The '--trim_poly_g' flag is 
        # a trimmer for polyG tails. For Illumina NextSeq/NovaSeq data, polyG can happen in read tails since G means no 
        # signal in the Illumina two-color systems. fastp can detect the polyG in read tails and trim them. This feature 
        # is enabled for NextSeq/NovaSeq data by default, and you can specify -g or --trim_poly_g. This value is 10 by default.
        "fastp "
        "--in1 {input.read_1} --in2 {input.read_2} "
        "--out1 {output.out_1} --out2 {output.out_2} "
        "--html {output.html_report} --json {output.json_report} "
        "--qualified_quality_phred 25 "
        "--length_required 25 "
        "--detect_adapter_for_pe "
        "--trim_poly_x "
        "--trim_poly_g"