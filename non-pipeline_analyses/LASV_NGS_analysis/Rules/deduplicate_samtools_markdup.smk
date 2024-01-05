
# Description:
# Remove duplicates using samtools markdup command.

# Author:
# Caleb Carr


rule convert_sam_file_to_bam_file_dedup:
    """ This rule converts the aligned sam file to a bam file using samtools view. """
    input:
        # Function in common.smk to retrieve the sam files
        get_mapped_sam_files(),
    output:
        temp("Data/Samples/Mapped_Reads/{accession}/{accession}_dedup.bam"),
    conda:
        "../Conda_Envs/post_alignment_processing.yml"
    shell:
        # The 'b' flag means that the output is a bam file while
        # the 'S' flag means that the input is auto-detected as sam.
        "samtools view -S -b {input} -o {output}"


rule collate_bam_file_dedup:
    """ This rule converts groups reads together by name using samtools collate. """
    input:
        "Data/Samples/Mapped_Reads/{accession}/{accession}_dedup.bam",
    output:
        temp("Data/Samples/Mapped_Reads/{accession}/{accession}_collate_dedup.bam"),
    conda:
        "../Conda_Envs/post_alignment_processing.yml"
    shell:
        "samtools collate {input} -o {output}"


rule fixmate_bam_file_dedup:
    """ This rule adds mate scores using the fixmate command. This is required 
        to run samtools markdup. """
    input:
        "Data/Samples/Mapped_Reads/{accession}/{accession}_collate_dedup.bam",
    output:
        temp("Data/Samples/Mapped_Reads/{accession}/{accession}_fixmate_dedup.bam"),
    conda:
        "../Conda_Envs/post_alignment_processing.yml"
    shell:
        # The '-m' flag signals to include a mate score.
        "samtools fixmate -m {input} {output}"


rule sort_bam_after_fixmate_dedup:
    """ This rule sorts the bam file after the mate scores have been added. """
    input:
        "Data/Samples/Mapped_Reads/{accession}/{accession}_fixmate_dedup.bam",
    output:
        temp("Data/Samples/Mapped_Reads/{accession}/{accession}_sorted_fixmate_dedup.bam"),
    conda:
        "../Conda_Envs/post_alignment_processing.yml"
    shell:
        "samtools sort {input} -o {output}"


rule remove_duplicates_samtools_dedup:
    """ This rule removes duplicates from the bam file using samtools markdup. """
    input:
        "Data/Samples/Mapped_Reads/{accession}/{accession}_sorted_fixmate_dedup.bam",
    output:
        "Data/Samples/Mapped_Reads/{accession}/{accession}_sorted_dedup.bam",
    conda:
        "../Conda_Envs/post_alignment_processing.yml"
    shell:
        # The '-r' flag signals to remove duplicates while 
        # the '-s' flag signals to print some basic stats. 
        "samtools markdup -r -s {input} {output}"