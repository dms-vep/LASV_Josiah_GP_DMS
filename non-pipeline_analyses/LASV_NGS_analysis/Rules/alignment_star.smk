
# Description:
# Aligns to best matched reference genome for each sample using star.

# Author:
# Caleb Carr


rule star_map_reads_to_ref_genome:
    """ This rule runs STAR to align the sample reads to the best matched reference 
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
        output_sam = "Data/Samples/Mapped_Reads/{accession}/{accession}_STAR.sam"
    params:
        # Function in common.smk that retrieves the best matched reference 
        # for a given sample. This function only returns the directory of the
        # reference genome because the ref directory must be accessed from 
        # that directory for bbmap to work.
        most_likely_ref = lambda wildcards: get_most_likely_ref_index(wildcards.accession),
        # Output file from STAR that must be renamed
        output_sam = "Data/Samples/Mapped_Reads/{accession}/Aligned.out.sam",
        # Output directory for STAR output files
        output_star = "Data/Samples/Mapped_Reads/{accession}/",
    conda:
        "../Conda_Envs/genome_processing_env.yml"
    shell:
        # The '--runMode' flag specifies what mode is being run
        # in STAR, which in this case is aligning reads. The 
        # '--genomeDir' flag specifies the indexed reference
        # genome which is the most likely reference found.
        # The '--readFilesIn' flag specifes the input read
        # files to be aligned. The '--twopassMode Basic' 
        # flag means that STAR does two passes aligning
        # the reads which helps identify splice junctions
        # because the first pass identifies possible junctions
        # and the second pass continues to align the reads 
        # with the found junctions. The '--scoreGapNoncan'
        # flag specifies how much new splice junction sites
        # are penalized. Setting this to 0 should make it 
        # just as likely to find de novo junctions. The
        # '--outFileNamePrefix' flag specifies the output
        # directory for the aligned files and the change
        # name command renames the aligned sam file to
        # the output name specified in the rule. The 
        # '--outSAMunmapped Within KeepPairs' will record
        # unmapped reads in the output sam file.
        "STAR "
        "--runMode alignReads "
        "--genomeDir {params.most_likely_ref} "
        "--readFilesIn {input.read_1} {input.read_2} "
        "--twopassMode Basic "
        "--scoreGapNoncan 0 "
        "--outFileNamePrefix {params.output_star} "
        "--outSAMunmapped Within KeepPairs "
        "&& "
        "mv {params.output_sam} {output.output_sam}"
