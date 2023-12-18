# Description:
# Python script that calculates all
# mutations along a phylogenetic tree

# Author:
# Caleb Carr

# Imports
import ete3
import pandas as pd
import numpy as np
from Bio import AlignIO

# Functions
def parse_tree_sequences(tree_sequences, tree_file, codon_alignment):

    # Dictionary of codon values and corresponding amino acids
    codon_chart = {
        'TTT': 'Phe', 'TTC': 'Phe', 'TTA': 'Leu', 'TTG': 'Leu',
        'CTT': 'Leu', 'CTC': 'Leu', 'CTA': 'Leu', 'CTG': 'Leu',
        'ATT': 'Ile', 'ATC': 'Ile', 'ATA': 'Ile', 'ATG': 'Met',
        'GTT': 'Val', 'GTC': 'Val', 'GTA': 'Val', 'GTG': 'Val',
        'TCT': 'Ser', 'TCC': 'Ser', 'TCA': 'Ser', 'TCG': 'Ser',
        'CCT': 'Pro', 'CCC': 'Pro', 'CCA': 'Pro', 'CCG': 'Pro',
        'ACT': 'Thr', 'ACC': 'Thr', 'ACA': 'Thr', 'ACG': 'Thr',
        'GCT': 'Ala', 'GCC': 'Ala', 'GCA': 'Ala', 'GCG': 'Ala',
        'TAT': 'Tyr', 'TAC': 'Tyr', 'TAA': 'Stp', 'TAG': 'Stp',
        'CAT': 'His', 'CAC': 'His', 'CAA': 'Gln', 'CAG': 'Gln',
        'AAT': 'Asn', 'AAC': 'Asn', 'AAA': 'Lys', 'AAG': 'Lys',
        'GAT': 'Asp', 'GAC': 'Asp', 'GAA': 'Glu', 'GAG': 'Glu',
        'TGT': 'Cys', 'TGC': 'Cys', 'TGA': 'Stp', 'TGG': 'Trp',
        'CGT': 'Arg', 'CGC': 'Arg', 'CGA': 'Arg', 'CGG': 'Arg',
        'AGT': 'Ser', 'AGC': 'Ser', 'AGA': 'Arg', 'AGG': 'Arg',
        'GGT': 'Gly', 'GGC': 'Gly', 'GGA': 'Gly', 'GGG': 'Gly',
    }

    # Initialize results df
    mutation_counts = pd.DataFrame(columns=["site", "mutation", "mut_type"])

    # Read tree seq data
    tree_seqs = pd.read_csv(
        tree_sequences,  
        sep="\t",
        comment="#",
        usecols=["Node", "Site", "State"]
    )

    # Read codon alignment file
    alignment = AlignIO.read(codon_alignment, "fasta")
    


    # Create new column with node number
    # tree_seqs["Node Number"] = tree_seqs["Node"].apply(lambda x: int(x[4:]))

    # Create ete3 tree object
    tree = ete3.Tree(tree_file, format=1)

    # Iterate through tree
    for node in tree.traverse("preorder"):
        
        # Ignore the first node its the outgroup
        if not node.is_root():

            # Get parent node name and seq
            parent = node.up.name
            parent_seq = tree_seqs.query("Node == @parent")["State"].tolist()

            # Get current node name
            curr_node = node.name
            curr_seq = None

            if not node.is_leaf():
                # Get current node seq
                curr_seq = tree_seqs.query("Node == @curr_node")["State"].tolist()
            else:
                # Get leaf node sequences from alignment
                for fasta in alignment:
                    if fasta.id == curr_node:
                        fasta_seq = str(fasta.seq)
                        curr_seq = [fasta_seq[start:start+3] for start in range(0, len(fasta_seq), 3)]
                # Skip if sequence is not in alignment (i.e., outgroup)
                if curr_seq == None:
                    continue

            # Check to make sure all sequences are the same length
            assert len(parent_seq) == len(curr_seq), "Sequences not same length!"

            # Iterate through sequences 
            for i in range(len(parent_seq)):
                if parent_seq[i] != curr_seq[i]:

                    # Get current mutation
                    mutation = parent_seq[i] + str(i+1) + curr_seq[i]

                    # Classify mutation type
                    mut_type = None
                    if codon_chart.get(parent_seq[i], parent_seq[i]) == codon_chart.get(curr_seq[i], curr_seq[i]):
                        mut_type = "synonymous"
                    elif codon_chart.get(curr_seq[i], curr_seq[i]) == "Stp":
                        mut_type = "nonsense"
                    else:
                        mut_type = "nonsynonymous"

                    # Add new mutation to df
                    mutation_counts.loc[len(mutation_counts.index)] = [i+1, mutation, mut_type] 


    return mutation_counts


def main():
    """
    Main method
    """

    # Input files
    tree_sequences = str(snakemake.input.tree_sequences)
    tree_file = str(snakemake.input.tree)
    codon_alignment = str(snakemake.input.codon_alignment)

    # Output files
    mutation_counts = str(snakemake.output.mutation_counts)

    result_df = parse_tree_sequences(tree_sequences, tree_file, codon_alignment)

    # Write output file
    result_df.to_csv(mutation_counts, index=False)



if __name__ == "__main__":
    main()