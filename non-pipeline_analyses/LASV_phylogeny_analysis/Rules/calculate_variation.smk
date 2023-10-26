# Description:
# Calculates entropy and n effective
# states for each residue position in
# protein alignment

# Author:
# Caleb Carr


rule calculate_variation:
    """
    This calculates site level entropy and n effective 
    amino acids based on natural protein sequences.
    """
    input:
        protein_alignment = config["GPC_protein_alignment"],
    output:
        config["GPC_protein_variation"],
    script:
        "../Scripts/calculate_variation.py"