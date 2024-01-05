
# Description:
# Downloads specified samples from csv file listed in the config.yml file.
# This most likely will be the run selector table from a BioProject. 

# Author:
# Caleb Carr


rule prefetch_download_SRA:
    """ This rule downloads the given accession using prefetch.  """
    output:
        # SRA file
        temp("Data/Samples/Raw_Reads/{accession}/{accession}.sra"),
    conda:
        "../Conda_Envs/sra_env.yml"
    shell:
        # Downloads accession number to output. The flag 
        # '-p' shows progression of download while '--output-file'
        # specifies the output file.
        "prefetch "
        "-p {wildcards.accession} --output-file {output}"


rule fasterq_dump_sra_conversion_PE:
    """ This rule converts the downloaded sra data file into
        the individual read fastq files. """
    input:
        "Data/Samples/Raw_Reads/{accession}/{accession}.sra"
    output:
        # Function in common.smk to retrieve fastq files
        read_1 = get_fastq_files_PE()[0],
        read_2 = get_fastq_files_PE()[1],
    params:
        input_dir = "Data/Samples/Raw_Reads/{accession}",
        output_dir = "Data/Samples/Raw_Reads/{accession}/",
    conda:
        "../Conda_Envs/sra_env.yml"
    shell:
        # The '-p' flag shows progress while the '--split-3' 
        # flag is for paired end reads. The ouput files are
        # directed to the directory specified by the '--outdir'
        # flag. The '--temp' flag specifies where the intermediate
        # files are stored before being deleted after the file
        # has been successfully converted. 
        "fasterq-dump "
        "-p --split-3 {params.input_dir} --outdir {params.output_dir} --temp {params.output_dir}"
