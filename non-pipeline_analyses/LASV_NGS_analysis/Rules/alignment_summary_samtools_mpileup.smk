
# Description:
# Create a summary file from the aligned sam file which
# can be parsed to get variant information for every position
# using samtools mpileup.

# Author:
# Caleb Carr


rule samtools_convert_sam_file_to_bam_file:
    """ This rule converts the aligned sam file to a bam file using samtools view. """
    input:
        # Function in common.smk to retrieve the sam files
        get_mapped_sam_files(),
    output:
        temp("Data/Samples/Mapped_Reads/{accession}/{accession}.bam"),
    conda:
        "../Conda_Envs/post_alignment_processing.yml"
    shell:
        # The 'b' flag means that the output is a bam file while
        # the 'S' flag means that the input is auto-detected as sam.
        "samtools view -S -b {input} -o {output}"


rule samtools_alignment_mapping_statistics:
    """ This rule runs samtools stats on the sorted bam to 
        get mapping statistics which will be used to get
        mapping rates. """
    input:
        "Data/Samples/Mapped_Reads/{accession}/{accession}.bam",
    output:
        temp("Results/Alignment_Stats/{accession}/{accession}_samtools_stats.txt"),
    conda:
        "../Conda_Envs/post_alignment_processing.yml"
    shell:
        "samtools stats {input} | grep ^SN | cut -f 2- &> {output}"


rule samtools_sort_bam_file:
    """ This rule sorts the bam file using samtools sort """
    input:
        "Data/Samples/Mapped_Reads/{accession}/{accession}.bam",
    output:
        "Data/Samples/Mapped_Reads/{accession}/{accession}_sorted.bam",
    conda:
        "../Conda_Envs/post_alignment_processing.yml"
    shell:
        # The input is the bam file of the aligned data
        # and the output is the sorted bam file.
        "samtools sort {input} -o {output}"


rule samtools_create_mpileup:
    """ This rule creates a mpileup summary of the aligned data which can
        be parsed to produce variant information. """
    input:
        # Need all reference genomes because it is not known
        # which reference goes with each sample
        expand("Data/Reference_Genomes/{ref_name}/{ref_name}.fasta.fai", ref_name=config['References']),
        # Function in common.smk to retrieve sorted bam files
        sorted_bam = get_sorted_bam(),
    output:
        "Data/Samples/Mapped_Reads/{accession}/{accession}_mpileup.txt",
    params:
        most_likely_ref = lambda wildcards: get_most_likely_ref_genome(wildcards.accession)
    conda:
        "../Conda_Envs/post_alignment_processing.yml"
    shell:
        # The input is the sorted bam file while the output is 
        # the mpileup text file. The '-d' flag signals the max 
        # depth is created for each position in the genome. Setting 
        # it to 0 means that the largest int value is used which is
        # needed for high coverage data. The '-q' flag signals for 
        # minimum map quality score (phred score) to be considered.
        # The '-Q' flag signals for a minimum base quality (phred score) 
        # to be to be considered. The '-B' flag disables the per-Base 
        # Alignment Quality scores that mpileup calculates when a reference
        # is used. If this is turned off, then the scores stored in the sam
        # file are used. When BAQ is turned on, too many reads can be filtered 
        # out. The '-f' flag means a reference is passed. The reference must 
        # be indexed and have the .fai file in the same directory as the 
        # regular fasta file. 
        # Changed to match andersen et al paper to verify N89D mutation in isolate
        "samtools mpileup -d 0 -q 1 -Q 30 -B -f {params.most_likely_ref} {input.sorted_bam} -o {output}"




