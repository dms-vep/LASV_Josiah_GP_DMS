
# Description:
# Aligns to best matched reference genome for each sample using bbmap.

# Author:
# Caleb Carr


rule bbmap_map_reads_to_ref_genome:
    """ This rule runs BBMap to align the sample reads to the best matched reference genome
        as listed in the report that checked alignment to all reference genomes. """
    input:
        # Have all indexed reference genomes listed as input because it is not known 
        # which reference will be needed to align the sample. The *_index_log.txt files
        # are dummy files to use as input for this rule. Function in common.smk retrieves
        # all index logs for the reference files.
        get_ref_index_log_files(),
        # Common function to return report for best matching reference
        get_ref_report(),
        # Function in common.smk to retrieve all cleaned fastq files
        read_1 = get_cleaned_fastq_files_PE()[0],
        read_2 = get_cleaned_fastq_files_PE()[1],
    output:
        output_sam = "Data/Samples/Mapped_Reads/{accession}/{accession}_BBMAP.sam",
    params:
        # Function in common.smk that retrieves the best matched reference 
        # for a given sample. This function only returns the directory of the
        # reference genome because the ref directory must be accessed from 
        # that directory for bbmap to work.
        most_likely_ref = lambda wildcards: get_most_likely_ref_index(wildcards.accession),
        # Soft limit on maximum indel sizes specified in the config file
        maxindel = config['Alignment_Mapping']['BBMap_maxindel'],
    conda:
        "../Conda_Envs/genome_processing_env.yml"
    shell:
        # The 'in1' and 'in2' inputs are the cleaned sample reads while the 
        # 'path' flag has the directory of the reference genome as determined 
        # by the 'get_most_likely_ref_genomes' function in the common.smk file.
        # The 'out' flag specifies where the final sam final is placed and the 
        # 'statsfile' flag specifies that the standard output with stats about 
        # the alignment is stored in a text file. The 'scafstats' flag specifies 
        # that the stats about which reads mapped to what parts of the reference 
        # genome are also stored in a text file. The 'maxindel' flag sets a soft
        # limit on the maximum size indels that will be found. In this case, the 
        # the limit is set to just under the size of the largest gene in the 
        # reference. 
        "bbmap.sh "
        "in1={input.read_1} in2={input.read_2} "
        "path={params.most_likely_ref} "
        "out={output.output_sam} "
        "maxindel={params.maxindel}"

    