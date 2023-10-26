# Description:
# Removes duplicates from alignments

# Author:
# Caleb Carr

rule remove_duplicate_sequences:
    """
    This rule removes duplicate sequences 
    from the segment, protein, and codon 
    alignments.
    """
    input:
        # Dummy file to make sure removal is done after
        # variation calculations
        protein_variation = config["GPC_protein_variation"],
        segment_alignment = config["S_segment_alignment"],
        codon_alignment = config["GPC_codon_alignment"],
        protein_alignment = config["GPC_protein_alignment"],
    params:
        outgroup = config["Outgroup"],
    output:
        log = config["Duplicate_removal_log"],
        segment_alignment = config["S_segment_alignment_deduplicated"],
        codon_alignment = config["GPC_codon_alignment_deduplicated"],
        protein_alignment = config["GPC_protein_alignment_deduplicated"]
    script:
        "../Scripts/remove_duplicates.py"


