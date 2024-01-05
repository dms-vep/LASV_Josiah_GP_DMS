
# Description:
# Downloads specified contaminant genomes as listed in the config.yml file as well 
# as creating the required mask and index.

# Author:
# Caleb Carr


rule download_contaminant_genomes:
    """ This rule uses urls to download the contaminant genomes 
        and store it in the contaminants directory. The urls to 
        download the contaminant genomes are stored in the config.yml 
        file. """
    output:
        # Only keep the compressed fasta file until the masked 
        # fasta file is created.
        temp("Data/Contaminant_Genomes/{contaminant_name}.fasta.gz")
    params:
        # Sample contaminant genomes are listed in the config file
        ncbi_urls = lambda wildcards: config['Sample_Contaminants'][wildcards.contaminant_name]['Download_url']
    shell:
        # The input are the urls listed in the config.yml file
        # while the output is the downloaded compressed fasta file.
        "wget {params.ncbi_urls} -O {output}"


rule bbmask_mask_contaminant_genomes:
    """ This rule runs BBMask on the contaminant genomes to prepare to 
        remove matching reads. BBMask masks regions of low quality and 
        highly repeptitive regions in addition to having the option for 
        inputing sam files for specific organisms which reduces the 
        chances of false positives during contaminant removal. """
    input:
        "Data/Contaminant_Genomes/{contaminant_name}.fasta.gz"
    output:
        "Data/Contaminant_Genomes/{contaminant_name}.fasta"
    conda:
        "../Conda_Envs/genome_processing_env.yml"
    shell:
        # The input is the newly downloaded genome while
        # output is the masked fasta version of the genome.
        "bbmask.sh in={input} out={output}"


rule bbsplit_index_contaminant_genomes:
    """ This rule runs bbsplit to created a merged index 
        for all contaminant genomes. """
    input:
        # Function in common.smk to retrieve all contaminant genome fasta files
        contaminant_genomes = get_list_of_contaminant_genomes()
    output:
        # Directory where the merged reference files are stored
        ref_dir = directory("Data/Contaminant_Genomes/BBSplit_Index/ref"),
        # Output log to use as input to the remove 
        # contaminants rule because directories cannot be used
        # as input to rules
        output_log = "Data/Contaminant_Genomes/BBSplit_Index/Merged_Index_BBSplit_contaminant_index_log.txt",
    params:
        # The list of genomes is converted to a string and spaces are removed 
        # between list items because bbsplit 'ref=' cannot process lists with 
        # commas followed by spaces.
        contaminant_genomes_string = ','.join(get_list_of_contaminant_genomes()),
        # Directory where the merged reference files are stored
        ref_dir = "Data/Contaminant_Genomes/BBSplit_Index/"
    conda:
        "../Conda_Envs/genome_processing_env.yml"
    shell:
        # The 'ref' flag refers to the contaminant genomes to be indexed
        # while the 'path' flag refers to where the merged index is 
        # placed. The 'usemodulo' is required to reduce the amount of 
        # RAM needed because without, the human genome requires 24G 
        # compared to 12G with 'usemodulo'. This needs to be set for 
        # both indexing and mapping. The standard output is captured
        # to use as input for the mapping rule with the '&>' flag.
        "bbsplit.sh "
        "ref={params.contaminant_genomes_string} path={params.ref_dir} "
        "usemodulo "
        "&> {output.output_log}"

