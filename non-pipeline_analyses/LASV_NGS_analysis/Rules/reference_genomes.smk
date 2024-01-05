
# Description:
# Downloads reference genomes from text files listed in the config.yml file
# and then indexes them with bbmap for alignment with bbmap as well as
# indexing them with samtools for samtools mpileup.

# Author:
# Caleb Carr

rule move_ref_genomes:
    """ This rule moves the reference file from the configure directory
        to the data directory.
    """
    input:
        "Configure/Reference_Files/{ref_name}.fasta"
    output:
        "Data/Reference_Genomes/{ref_name}/{ref_name}.fasta"
    shell:
        "cp {input} {output}"


rule samtools_index_ref_genomes:
    """ Thia rule indexes the reference fasta files with samtools.
        This is needed to run samtools mpileup which requires a samtools
        indexed reference genome. """
    input:
        "Data/Reference_Genomes/{ref_name}/{ref_name}.fasta"
    output:
        "Data/Reference_Genomes/{ref_name}/{ref_name}.fasta.fai"
    conda:
        "../Conda_Envs/genome_processing_env.yml"
    shell:
        # The input is the regular fasta file while the
        # output is the indexed .fai file that must be 
        # present in the same directory for samtools.
        "samtools faidx {input} -o {output}"


rule bbsplit_index_merged_ref_genomes:
    """ This rule runs bbsplit to create the merged
        index for all reference genomes which is used
        when determining which reference each sample 
        maps to. """
    input:
        # Function in common.smk to retrieve all reference genome fastas
        reference_genomes = get_list_of_reference_genomes(),
    output:
        # Temporary index log created to signal that the reference
        # merged index has been created
        output_log = "Data/Reference_Genomes/Merged_Ref_BBSPLIT/merged_ref_index_BBSPLIT.txt",
         # Directory where the merged reference files are placed
        ref_dir = directory("Data/Reference_Genomes/Merged_Ref_BBSPLIT/ref"),
    params:
        # The list of genomes is converted to a string and spaces are removed 
        # between list items because bbsplit 'ref=' cannot process lists with 
        # commas followed by spaces.
        reference_genomes_string = ','.join(get_list_of_reference_genomes()),
         # Directory where the merged reference files are placed
        ref_dir = "Data/Reference_Genomes/Merged_Ref_BBSPLIT/"
    conda:
        "../Conda_Envs/genome_processing_env.yml"
    shell:
        # The 'ref' flag specifies the reference genomes to 
        # be indexed while the 'path' flag specifies the
        # output directory for the indexed genome. The '&>' 
        # just captures the standard output and error for 
        # the output log which is used as an input file 
        # for the checking reference rule.
        "bbsplit.sh "
        "ref={params.reference_genomes_string} path={params.ref_dir} "
        "&> {output.output_log}"


rule bbmap_index_ref_genomes:
    """ This rule indexes the downloaded reference genomes using bbmap 
        which is needed for alignment with bbmap. """
    input:
        "Data/Reference_Genomes/{ref_name}/{ref_name}.fasta"
    output:
        # BBMap creates a new directory called 'ref' that contains all
        # of the indexed files for the genome.
        output_dir = directory("Data/Reference_Genomes/{ref_name}/BBMap_Index/ref/"),
        # Temporarily create an output txt file to use as input to 
        # map_reads_to_ref_genome rule in alignmnet.smk. This is needed
        # because directories could not be used as input for rules. 
        output_log = "Data/Reference_Genomes/{ref_name}/BBMap_Index/{ref_name}_index_log_BBMAP.txt",
    params:
        out_dir = "Data/Reference_Genomes/{ref_name}/BBMap_Index/"
    conda:
        "../Conda_Envs/genome_processing_env.yml"
    shell:
        # The 'path' flag is where the indexed reference genome is
        # stored for future use. The '&>' just captures the standard 
        # output and error for the output log which is used as an 
        # input file for the alignment rule. 
        "bbmap.sh ref={input} path={params.out_dir} &> {output.output_log}"


rule bwa_index_ref_genomes:
    """ This rule runs bwa index to index all the reference
        genomes before aligning with bwa mem. """
    input:
        "Data/Reference_Genomes/{ref_name}/{ref_name}.fasta"
    output:
        # The index files created from BWA index are stored in the
        # BWA_index directory as well as a copy of the reference
        output_dir = "Data/Reference_Genomes/{ref_name}/BWA_Index/{ref_name}.fasta",
        # Temporarily create an output txt file to use as input to 
        # map_reads_to_ref_genome rule in alignmnet.smk. This is needed
        # because directories could not be used as input for rules. 
        output_log = "Data/Reference_Genomes/{ref_name}/BWA_Index/{ref_name}_index_log_BWA.txt",
    conda:
        "../Conda_Envs/genome_processing_env.yml"
    shell:
        # The fasta file must first be copied into the specified 
        # directory and then it can be indexed with bwa. The '-a'
        # flag specifies the algorithm used to index the reference
        # while the 'bwtsw' is a more complex algorithm than the
        # default which can handle larger references, such as 
        # human genomes. The '&>' just captures the standard 
        # output and error for the output log which is used as an 
        # input file for the alignment rule.
        "cp {input} {output.output_dir} "
        "&& "
        "bwa index " 
        "-a bwtsw {output.output_dir} &> {output.output_log}"


rule star_index_ref_genomes:
    """ This rule runs star to index the reference genomes
        which are then used for alignment. """
    input:
        "Data/Reference_Genomes/{ref_name}/{ref_name}.fasta"
    output:
        # The index files created from STAR index are stored in the
        # STAR_index directory as well as a copy of the reference
        output_dir = "Data/Reference_Genomes/{ref_name}/STAR_Index/{ref_name}.fasta",
        # Temporarily create an output txt file to use as input to 
        # map_reads_to_ref_genome rule in alignmnet.smk. This is needed
        # because directories could not be used as input for rules. 
        output_log = "Data/Reference_Genomes/{ref_name}/STAR_Index/{ref_name}_index_log_STAR.txt",
    params:
        output_dir = "Data/Reference_Genomes/{ref_name}/STAR_Index/",
        gSAindexNbases = config['Alignment_Mapping']['STAR_genomeSAindexNbases']
    conda:
        "../Conda_Envs/genome_processing_env.yml"
    shell:
        # The fasta file must first be copied into the specified 
        # directory and then it can be indexed with STAR. The 
        # '--runMode' flag specifies the mode used by STAR while
        # the '--genomeDir' flag specifies the output directory
        # for the indexed genome. The '--genomeFastaFiles' flag
        # specifies the input fasta files for the reference to
        # be indexed. The '--genomeSAindexNbases' flag specifies
        # the adjustment for smaller genomes. The calculations to
        # adjust the parameter is found in the config file. The '&>' 
        # just captures the standard output and error for the 
        # output log which is used as an input file for the alignment 
        # rule.
        "cp {input} {output.output_dir} "
        "&& "
        "STAR "
        "--runMode genomeGenerate --genomeDir {params.output_dir} "
        "--genomeFastaFiles {output.output_dir} --genomeSAindexNbases {params.gSAindexNbases} "
        "&> {output.output_log}"


rule hisat2_index_ref_genomes:
    """ This rule runs hisat2 to index the reference genomes
        which are then used for alignment. """
    input:
        "Data/Reference_Genomes/{ref_name}/{ref_name}.fasta"
    output:
        # The index files created from HISAT2 index are stored in the
        # HISAT2_index directory as well as a copy of the reference
        output_dir = "Data/Reference_Genomes/{ref_name}/HISAT2_Index/{ref_name}.fasta",
        # Temporarily create an output txt file to use as input to 
        # map_reads_to_ref_genome rule in alignmnet.smk. This is needed
        # because directories could not be used as input for rules. 
        output_log = "Data/Reference_Genomes/{ref_name}/HISAT2_Index/{ref_name}_index_log_HISAT2.txt",
    params:
        # Output directory
        output_dir = "Data/Reference_Genomes/{ref_name}/HISAT2_Index/{ref_name}"
    conda:
        "../Conda_Envs/genome_processing_env.yml"
    shell:
        # The fasta file must first be copied into the specified 
        # directory and then it can be indexed with HISAT2. The 
        # first argument is the fasta file of the reference to be
        # indexed while the second argument is the directory for
        # the index files. The '&>' flag just captures the standard 
        # output and error for the output log which is used as an 
        # input file for the alignment rule.
        "cp {input} {output.output_dir} "
        "&& "
        "hisat2-build "
        "{output.output_dir} "
        "{params.output_dir} "
        "&> {output.output_log}"
        




