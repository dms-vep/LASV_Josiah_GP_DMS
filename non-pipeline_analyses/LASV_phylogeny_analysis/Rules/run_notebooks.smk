# Description:
# Runs jupyter notebooks

# Author:
# Caleb Carr

rule visualize_GPC_diversity:
    """
    This rule runs a notebook to create
    a GPC diversity figure.
    """
    input:
        config["GPC_protein_variation"]
    output:
        config["GPC_diversity_figure"]
    notebook:
        "../Notebooks/visualize_GPC_diversity.ipynb"

rule create_reduced_tree_figure:
    """
    This rule runs a notebook to create
    the final reduced codon tree figure.
    """
    input:
        protein=config["GPC_protein_reduced_tree_output"],
        codon=config["GPC_codon_reduced_tree_output"],
    output:
        config["Protein_validation_sequences"],
        config["GPC_protein_validation_alignment_figure"],
        config["GPC_codon_reduced_tree_figure"],
        config["Amino_acid_identity_to_josiah_figure"],
        config["GPC_protein_reduced_tree_figure"],
    notebook:
        "../Notebooks/visualize_reduced_tree.ipynb"