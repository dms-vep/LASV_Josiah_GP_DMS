
# Description:
# Parse the mpileup file to find variants and then classify
# the variants based on mutation type. In addition, parse
# the samtools stats file to get mapping rates for each 
# sample. 

# Author:
# Caleb Carr


rule python_get_mapping_rates:
    """ This rule parse the text file generated from samtools
        stats command to get mapping rates for each sample. """
    input:
        expand("Results/Alignment_Stats/{accession}/{accession}_samtools_stats.txt", accession=get_list_of_sample_accessions()),
    output:
        "Results/Alignment_Stats/sample_mapping_rates.csv",
    params:
        # List of accession numbers to use for sample 
        # names in the final mapping rate file
        samples = get_list_of_sample_accessions(),
    conda:
        "../Conda_Envs/data_file_parsing_env.yml"
    script:
        # Python script that parse the samtools stats
        # text files and calculates mapping rates for
        # each sample.
        "../Scripts/process_samtools_stats_files.py"


rule python_process_mpileup_file:
    """ This rule runs the python script to parse the mpileup file to get
        variant information for every position in the genome. """
    input:
        "Data/Samples/Mapped_Reads/{accession}/{accession}_mpileup.txt",
    output:
        "Results/Sample_Data/{accession}/{accession}_mpileup_summary.csv",
    conda:
        "../Conda_Envs/data_file_parsing_env.yml"
    script:
        # Python script that parses the mpileup text file to retrieve
        # the base counts for every position of the genome, indel info,
        # and read coverage.
        "../Scripts/process_mpileup_file.py"


rule python_get_minor_base_frequencies:
    """ This rule parses the mpileup summary file to retrieve base
        frequencies for each base as well as strand bias
        ratios and indel frequencies. """
    input:
        "Results/Sample_Data/{accession}/{accession}_mpileup_summary.csv",
    output:
        temp("Results/Sample_Data/{accession}/{accession}_minor_base_frequencies.csv"),
    conda:
        "../Conda_Envs/data_file_parsing_env.yml"
    script:
        # Python script that parses the mpileup summary csv and returns a
        # new file with base count frequencies, strand bias ratios, and
        # indel frequencies.
        "../Scripts/get_minor_base_frequencies.py"


rule python_get_variants:
    """ This rule extracts the variant bases for each position
        and indels based on frequency, bias, and depth cutoffs 
        specified in the config file. """
    input:
        "Results/Sample_Data/{accession}/{accession}_minor_base_frequencies.csv",
    output:
        temp("Results/Sample_Data/{accession}/{accession}_variants.csv"),
    params:
        # Variant parameters are passed through the config file
        minor_base_cutoff = config['Variant_Parameters']['Minor_Base_Frequency_cutoff'],
        insert_cutoff = config['Variant_Parameters']['Insertion_Frequency_cutoff'],
        del_cutoff = config['Variant_Parameters']['Deletion_Frequency_cutoff'],
        min_ratio_bias = config['Variant_Parameters']['Strand_Bias_Ratio_Min'],
        max_ratio_bias = config['Variant_Parameters']['Strand_Bias_Ratio_Max'],
        min_depth = config['Variant_Parameters']['Min_Pos_Depth'],
    conda:
        "../Conda_Envs/data_file_parsing_env.yml"
    script:
        # Python script that parses the frequencies summary file and
        # returns a new file with just the positions that meet the 
        # variant cutoffs as specified in the config file
        "../Scripts/get_variants.py"


rule python_classify_variants:
    """ This rule classifies each variant found as either synonymous,
        nonsynonymous, nonsense, UTR, or indel. """
    input:
        alignment_data = "Results/Sample_Data/{accession}/{accession}_minor_base_frequencies.csv",
        variant_data = "Results/Sample_Data/{accession}/{accession}_variants.csv",
        coding_regions_indices = "Results/Sample_Data/{accession}/{accession}_protein_indices.csv",
    output:
        "Results/Sample_Data/{accession}/{accession}_classified_variants.csv",
    conda:
        "../Conda_Envs/data_file_parsing_env.yml"
    script:
        # Python script that classifies each variant
        "../Scripts/classify_variants.py"


rule python_get_common_variants:
    """ This rule extracts counts for all variants across
        all samples based on position in genome. """
    input:
        variants = expand("Results/Sample_Data/{accession}/{accession}_classified_variants.csv", accession=get_list_of_sample_accessions()),
        alignemnt = expand("Results/Sample_Data/{accession}/{accession}_mpileup_summary.csv", accession=get_list_of_sample_accessions()),
    output:
        "Results/Sample_Data/Common_Variants/common_variants.csv",
    params:
        # Optional gene names settings and names list as 
        # specified in the config file
        gene_name_setting = config['Optional_Gene_Names']['Provide_Gene_Names'],
        gene_names = config['Optional_Gene_Names']['Gene_Names'],
        gene_list = config['Gene_list'],
        all_nonsyn_variants = "Results/Sample_Data/Common_Variants/all_nonsyn_variants.csv",
    conda:
        "../Conda_Envs/data_file_parsing_env.yml"
    script:
        "../Scripts/get_common_variants.py"