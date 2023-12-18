# Description:
# Rules to calculate diversity
# within natural Lassa sequences

# Author:
# Caleb Carr


rule calculate_variation:
    """
    This rule calculates site level entropy and n effective 
    amino acids based on natural protein sequences.
    """
    input:
        protein_alignment = config["GPC_protein_alignment"],
    output:
        config["GPC_protein_variation"],
    script:
        "../Scripts/calculate_variation.py"


rule count_tree_mutations:
    """
    This rule calculates mutations along the
    reconstructed ancestors of the phylogenetic tree.
    """
    input:
        tree_sequences = config["GPC_codon_tree_with_outgroup_ancestors"],
        tree = config["GPC_codon_tree_with_outgroup_prefix"] + ".treefile",
        codon_alignment = config["GPC_codon_alignment_deduplicated"],
    output:
        mutation_counts = config["GPC_tree_mutations"],
    script:
        "../Scripts/count_tree_mutations.py"


rule run_FUBAR:
    """
    This rule runs the hyphy program FUBAR
    """
    input:
        tree = config["GPC_codon_tree_prefix"] + ".treefile",
        alignment = config["GPC_codon_alignment_deduplicated"],
    output:
        config["GPC_FUBAR_results"],
    shell:
        "hyphy fubar --alignment {input.alignment}  --tree {input.tree} --output {output}"


rule run_FEL:
    """
    This rule runs the hyphy program FEL
    """
    input:
        tree = config["GPC_codon_tree_prefix"] + ".treefile",
        alignment = config["GPC_codon_alignment_deduplicated"],
    output:
        config["GPC_FEL_results"],
    shell:
        "hyphy fel --alignment {input.alignment}  --tree {input.tree} --output {output}"