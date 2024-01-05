
# Description:
# Aligns to best matched reference genome for each sample using bwa-mem.
# Note on the BWA algorithms:
# BWA is a software package for mapping DNA sequences against a large 
# reference genome, such as the human genome. It consists of three 
# algorithms: BWA-backtrack, BWA-SW and BWA-MEM. The first algorithm is 
# designed for Illumina sequence reads up to 100bp, while the rest two 
# for longer sequences ranged from 70bp to a few megabases. BWA-MEM and 
# BWA-SW share similar features such as the support of long reads and 
# chimeric alignment, but BWA-MEM, which is the latest, is generally 
# recommended as it is faster and more accurate. BWA-MEM also has better 
# performance than BWA-backtrack for 70-100bp Illumina reads.

# Author:
# Caleb Carr


rule bwa_mem_map_reads_to_ref_genome:
    """ This rule runs BWA Mem to align the sample reads to the best matched reference 
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
        output_sam = "Data/Samples/Mapped_Reads/{accession}/{accession}_BWA.sam"
    params:
        # Function in common.smk that retrieves the best matched reference 
        # for a given sample. This function only returns the directory of the
        # reference genome because the ref directory must be accessed from 
        # that directory for bbmap to work.
        most_likely_ref = lambda wildcards: get_most_likely_ref_index(wildcards.accession)
    conda:
        "../Conda_Envs/genome_processing_env.yml"
    shell:
        "bwa mem "
        "{params.most_likely_ref} {input.read_1} {input.read_2} > {output.output_sam}"