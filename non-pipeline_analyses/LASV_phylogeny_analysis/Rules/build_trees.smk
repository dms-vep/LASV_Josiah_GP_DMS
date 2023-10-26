# Description:
# Builds phylogeny trees for S segment, protein, and 
# codon sequences using IQtree

# Author:
# Caleb Carr


rule build_S_segment_tree:
    """
    This rule builds tree using IQtree for 
    all S segment sequences
    """
    input:
        alignment = config["S_segment_alignment_deduplicated"],
        log = config["Duplicate_removal_log"],
    params:
        prefix = config["S_segment_tree_prefix"],
        outgroup = config["Outgroup"],
    output:
        config["S_segment_tree_output"],
    shell:
        #-s     is the option to specify the name of the alignment file
        #-pre   Specify a prefix for all output files. DEFAULT: either alignment file
        #       name (-s) or partition file name (-q, -spp or -sp)
        #-o     Specify an outgroup taxon name to root the tree. The output tree
        #       in .treefile will be rooted accordingly. DEFAULT: first taxon in
        #       alignment
        #-nt    Specify the number of CPU cores for the multicore version. A
        #       special option -nt AUTO will tell IQ-TREE to automatically
        #       determine the best number of cores given the current data and
        #       computer.
        #-m     is the option to specify the model name to use during the analysis.
        #       'MFP' is the default model finder which will try different ones and 
        #       choose the best fitting one
        #-quiet Silent mode, suppress printing to the screen. Note that .log file is
        #       still written.
        "iqtree -s {input.alignment} -pre {params.prefix} -nt 10 -m MFP -o {params.outgroup} -quiet"


rule build_GPC_protein_tree:
    """
    This rule builds tree using IQtree for 
    all GPC protein sequences
    """
    input:
        alignment = config["GPC_protein_alignment_deduplicated"],
        log = config["Duplicate_removal_log"],
    params:
        prefix = config["GPC_protein_tree_prefix"],
        outgroup = config["Outgroup"],
    output:
        config["GPC_protein_tree_output"],
    shell:
        #-s     is the option to specify the name of the alignment file
        #-pre   Specify a prefix for all output files. DEFAULT: either alignment file
        #       name (-s) or partition file name (-q, -spp or -sp)
        #-o     Specify an outgroup taxon name to root the tree. The output tree
        #       in .treefile will be rooted accordingly. DEFAULT: first taxon in
        #       alignment
        #-nt    Specify the number of CPU cores for the multicore version. A
        #       special option -nt AUTO will tell IQ-TREE to automatically
        #       determine the best number of cores given the current data and
        #       computer.
        #-m     is the option to specify the model name to use during the analysis.
        #       'LG' is an emperical amino-acid model (https://pubmed.ncbi.nlm.nih.gov/18367465/)
        #       '+F' is the empirical codon frequencies counted from the data. In AliSim, if users
        #       neither specify base frequencies nor supply an input alignment, AliSim
        #       will generate base frequencies from empirical distributions.
        #       '+G' Discrete Gamma model (Yang, 1994) with default 4 rate categories.
        #       The number of categories can be changed with e.g. +G8.
        #-quiet Silent mode, suppress printing to the screen. Note that .log file is
        #       still written.
        "iqtree -s {input.alignment} -pre {params.prefix} -nt 10 -m LG+F+G -o {params.outgroup} -quiet"

rule build_reduced_GPC_protein_tree:
    """
    This rule builds tree using IQtree for 
    a subet of GPC protein sequences
    """
    input:
        alignment = config["GPC_protein_reduced_alignment"],
    params:
        prefix = config["GPC_protein_reduced_tree_prefix"],
        outgroup = config["Outgroup"],
    output:
        config["GPC_protein_reduced_tree_output"],
    shell:
        #-s     is the option to specify the name of the alignment file
        #-pre   Specify a prefix for all output files. DEFAULT: either alignment file
        #       name (-s) or partition file name (-q, -spp or -sp)
        #-o     Specify an outgroup taxon name to root the tree. The output tree
        #       in .treefile will be rooted accordingly. DEFAULT: first taxon in
        #       alignment
        #-nt    Specify the number of CPU cores for the multicore version. A
        #       special option -nt AUTO will tell IQ-TREE to automatically
        #       determine the best number of cores given the current data and
        #       computer.
        #-m     is the option to specify the model name to use during the analysis.
        #       'LG' is an emperical amino-acid model (https://pubmed.ncbi.nlm.nih.gov/18367465/)
        #       '+F' is the empirical codon frequencies counted from the data. In AliSim, if users
        #       neither specify base frequencies nor supply an input alignment, AliSim
        #       will generate base frequencies from empirical distributions.
        #       '+G' Discrete Gamma model (Yang, 1994) with default 4 rate categories.
        #       The number of categories can be changed with e.g. +G8.
        #-quiet Silent mode, suppress printing to the screen. Note that .log file is
        #       still written.
        # -czb    Collapse near zero branches, so that the final tree may be
        #         multifurcating. This is useful for bootstrapping in the presence of
        #         polytomy to reduce bootstrap supports of short branches.
        # -bb   Specify number of bootstrap replicates (>=1000).
        "iqtree -s {input.alignment} -pre {params.prefix} -nt 6 -m LG+F+G -o {params.outgroup} -quiet -czb -bb 1000"


rule build_reduced_GPC_codon_tree:
    """
    This rule builds tree using IQtree for 
    all GPC codon sequences
    """
    input:
        alignment = config["GPC_codon_reduced_alignment"],
    params:
        prefix = config["GPC_codon_reduced_tree_prefix"],
        outgroup = config["Outgroup"],
    output:
        config["GPC_codon_reduced_tree_output"],
    shell:
        #-s     is the option to specify the name of the alignment file
        #-pre   Specify a prefix for all output files. DEFAULT: either alignment file
        #       name (-s) or partition file name (-q, -spp or -sp)
        #-o     Specify an outgroup taxon name to root the tree. The output tree
        #       in .treefile will be rooted accordingly. DEFAULT: first taxon in
        #       alignment
        #-nt    Specify the number of CPU cores for the multicore version. A
        #       special option -nt AUTO will tell IQ-TREE to automatically
        #       determine the best number of cores given the current data and
        #       computer.
        #-m     is the option to specify the model name to use during the analysis
        #       Codon models:
        #       '+MG' Nonsynonymous/synonymous (dn/ds) rate ratio (Muse and Gaut, 1994).
        #       '+MGK' Like MG with additional transition/transversion (ts/tv) rate ratio.
        #       Codon frequencies:
        #       '+F3X4' Unequal nucleotide frequencies and unequal nt frequencies over three
        #       codon positions. In AliSim, if users don’t supply an input alignment,
        #       the base frequencies are randomly generated based on empirical
        #       distributions, or users could specify the frequencies via
        #       '+F3X4'{<freq_0>,...,<freq_11>}
        #       '+G' Discrete Gamma model (Yang, 1994) with default 4 rate categories.
        #       The number of categories can be changed with e.g. +G8.
        #-st CODON codon model given a protein-coding DNA alignment
        #-quiet Silent mode, suppress printing to the screen. Note that .log file is
        #       still written.
        # -czb    Collapse near zero branches, so that the final tree may be
        #         multifurcating. This is useful for bootstrapping in the presence of
        #         polytomy to reduce bootstrap supports of short branches.
        # -bb   Specify number of bootstrap replicates (>=1000).
        "iqtree -s {input.alignment} -pre {params.prefix} -nt 10 -m MGK+G+F3X4 -st CODON -o {params.outgroup} -quiet -czb -bb 1000"


rule build_GPC_codon_tree:
    """
    This rule builds tree using IQtree for 
    all GPC codon sequences
    """
    input:
        alignment = config["GPC_codon_alignment_deduplicated"],
        log = config["Duplicate_removal_log"],
    params:
        prefix = config["GPC_codon_tree_prefix"],
        outgroup = config["Outgroup"],
    output:
        config["GPC_codon_tree_output"],
    shell:
        #-s     is the option to specify the name of the alignment file
        #-pre   Specify a prefix for all output files. DEFAULT: either alignment file
        #       name (-s) or partition file name (-q, -spp or -sp)
        #-o     Specify an outgroup taxon name to root the tree. The output tree
        #       in .treefile will be rooted accordingly. DEFAULT: first taxon in
        #       alignment
        #-nt    Specify the number of CPU cores for the multicore version. A
        #       special option -nt AUTO will tell IQ-TREE to automatically
        #       determine the best number of cores given the current data and
        #       computer.
        #-m     is the option to specify the model name to use during the analysis
        #       Codon models:
        #       '+MG' Nonsynonymous/synonymous (dn/ds) rate ratio (Muse and Gaut, 1994).
        #       '+MGK' Like MG with additional transition/transversion (ts/tv) rate ratio.
        #       Codon frequencies:
        #       '+F3X4' Unequal nucleotide frequencies and unequal nt frequencies over three
        #       codon positions. In AliSim, if users don’t supply an input alignment,
        #       the base frequencies are randomly generated based on empirical
        #       distributions, or users could specify the frequencies via
        #       '+F3X4'{<freq_0>,...,<freq_11>}
        #       '+G' Discrete Gamma model (Yang, 1994) with default 4 rate categories.
        #       The number of categories can be changed with e.g. +G8.
        #-st CODON codon model given a protein-coding DNA alignment
        #-quiet Silent mode, suppress printing to the screen. Note that .log file is
        #       still written.
        "iqtree -s {input.alignment} -pre {params.prefix} -nt 10 -m MGK+G+F3X4 -st CODON -o {params.outgroup} -quiet"