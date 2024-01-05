
# Description:
# Runs quality checking tools such as fastqc and then collects 
# all reports using multiqc. 

# Author:
# Caleb Carr


rule run_fastqc_PE:
    """ This rule runs fastqc on the downloaded sample fastqc files and outputs the 
        results for each sample in a folder in results under quality checking as 
        well as running fastqc on the cleaned reads after read processing. """
    input:
        # Function in common.smk to retrieve raw fastq files
        get_fastq_files_PE()[0],
        get_fastq_files_PE()[1],
        # Function in common.smk to retrieve cleaned fastq files
        get_cleaned_fastq_files_PE()[0],
        get_cleaned_fastq_files_PE()[1],
    output:
        # fastqc needs an output directory created first
        directory("Results/Quality_Checking/{accession}/Fastqc/")
    conda:
        "../Conda_Envs/qc_env.yml"
    shell:
        # Make the output directory and the '-p' flag allows for some of the parent
        # directories to already exist.
        "mkdir -p {output} "
        "&& "
        "fastqc {input} --outdir {output}"


rule run_multiqc:
    """ This rule runs multiqc which will search for quality checking documents, 
        such as htmls, and produce a consolidated version of all reports. """
    input:
        # Function in common.smk retrieves a list of all accession names which is
        # used in the snakemake expand function to create al input directories.
        expand("Results/Quality_Checking/{accession}/Fastqc/", accession=get_list_of_sample_accessions()),
    output:
        directory("Results/Quality_Checking/MultiQC")
    params:
        start_dir = directory("Results/Quality_Checking/")
    conda:
        "../Conda_Envs/qc_env.yml"
    shell:
        "multiqc {params.start_dir} -o {output}"
    