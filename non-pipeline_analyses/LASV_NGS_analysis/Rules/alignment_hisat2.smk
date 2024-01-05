
# Description:
# Aligns to best matched reference genome for each sample using hisat2.

# Author:
# Caleb Carr

rule hisat2_map_reads_to_ref_genome:
    """ This rule runs HISAT2 to align the sample reads to the best matched reference 
        genome as listed in the report that checked alignment to all reference genomes. """
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
        output_sam = "Data/Samples/Mapped_Reads/{accession}/{accession}_HISAT2.sam"
    params:
        # Function in common.smk that retrieves the best matched reference 
        # for a given sample. This function only returns the directory of the
        # reference genome because the ref directory must be accessed from 
        # that directory for bbmap to work.
        most_likely_ref = lambda wildcards: get_most_likely_ref_index(wildcards.accession),
    conda:
        "../Conda_Envs/genome_processing_env.yml"
    shell:
        # The '-x' flag specifies the directory for the index 
        # of the best matched reference. The '-1' and '-2' flags
        # specify the input read files for paired end data,
        # respectively. The '-S' flag specifies the output sam
        # file.
        "hisat2 "
        "-x {params.most_likely_ref} "
        "-1 {input.read_1} -2 {input.read_2} "
        "-S {output.output_sam}"