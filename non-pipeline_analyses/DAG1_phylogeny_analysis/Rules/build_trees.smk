# Description:
# Builds phylogeny trees for DAG1 protein sequences
# using IQtree

# Author:
# Caleb Carr


rule build_protein_tree:
    """
    This rule builds tree using IQtree for 
    DAG1 protein sequences
    """
    input:
        alignment = config["Protein_alignment"],
    params:
        prefix = config["Protein_tree_prefix"],
        outgroup = config["Outgroup"],
    output:
        config["Protein_tree_output"],
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
        "iqtree -s {input.alignment} -pre {params.prefix} -nt 10 -m LG+F+G -o {params.outgroup} -quiet -czb -bb 1000"