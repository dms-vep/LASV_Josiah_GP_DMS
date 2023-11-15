# Description:
# Runs jupyter notebooks

# Author:
# Caleb Carr

rule visualize_DAG1_MSA:
    """
    This rule runs a notebook to visualize
    DAG1 MSAs.
    """
    input:
        config["Protein_alignment"],
        config["Reduced_protein_alignment"],
    output:
        config["Consensus_glycosylation_site"],
        config["Human_and_mastomys_msa"],
    notebook:
        "../Notebooks/visualize_DAG1_MSA.ipynb"