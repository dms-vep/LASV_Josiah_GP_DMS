# Description:
# Python script that calculates site
# level entropy and n effective amino acids
# from protein alignment

# Author:
# Caleb Carr

# Imports
import pandas as pd
import numpy as np
from Bio import SeqIO, AlignIO

# Functions
def get_counts_from_multiple_sequence_alignment(alignment):
    """
    Calculate counts of each character in alignment.
    """
    
    # Intialize results
    matrix_0s = [[0 for _ in range(25)] for _ in range(491)]
    result_df = pd.DataFrame(matrix_0s, columns=[
        "site",
        "A",
        "R",
        "N",
        "D",
        "C",
        "E",
        "Q",
        "G",
        "H",
        "I",
        "L",
        "K",
        "M",
        "F",
        "P",
        "S",
        "T",
        "W",
        "Y",
        "V",
        "-",
        "X",
        "total_w_X_w_gaps",
        "total_no_X_no_gaps",
    ])

    # Add site numbering
    result_df["site"] = list(range(1, 492))

    # Get counts
    alignment = AlignIO.read(alignment, "fasta")
    alignment = alignment[:, :491]
    for seq in alignment:
        for i in range(len(seq.seq)):
            curr_char = seq.seq[i]
            # Update character count and total
            result_df.at[i, curr_char] += 1
            result_df.at[i, "total_w_X_w_gaps"] += 1
            # Update count excluding gaps and ambigous characters
            if curr_char != "-" and curr_char != "X":
                result_df.at[i, "total_no_X_no_gaps"] += 1  
            
    # Return results
    return result_df


def calculate_entropy_and_neffective_per_site(counts_df, with_no_X_no_gaps=True):
    """
    Calculate entropy and neffective amino acids per site as
    follows:

    entropy_at_site_x = −∑ πr,x*ln(πr,x) 
    where πr,x is the fraction of character r at site X

    and 

    neffective_at_site_x = e^(entropy_at_site_x)
    """

    sites = counts_df["site"]
    
    if with_no_X_no_gaps:
        # Remove ambiguous and gaps columns
        result_df = (
            counts_df.drop(
                columns=[
                    "total_w_X_w_gaps",
                    "-",
                    "X",
                    "site",
                ]
            )
        )
        # Calculate fraction per site
        result_df = (
            result_df
            .divide(result_df["total_no_X_no_gaps"], axis="rows")
            .drop(columns=["total_no_X_no_gaps"])
        )
        # Calculate frac*ln(frac) and sum across all chars
        for column in result_df.columns:
            result_df[column] = result_df[column] * np.log(result_df[column])
        result_df["entropy"] = result_df.sum(axis=1, numeric_only=True)
        # Multiply by -1
        result_df["entropy"] = result_df["entropy"]*-1
        # Number of effective states
        result_df["n_effective"] = np.exp(result_df["entropy"])
        # Reset site values
        result_df["site"] = sites
    else:
        # Remove total excluding ambiguous chars and gaps
        result_df = (
            counts_df.drop(
                columns=[
                    "total_no_X_no_gaps",
                ]
            )
        )
        # Calculate fraction per site
        result_df = (
            result_df
            .divide(result_df["total_w_X_w_gaps"], axis="rows")
            .drop(columns=["total_w_X_w_gaps"])
        )
        # Calculate frac*ln(frac) and sum across all chars
        for column in result_df.columns:
            result_df[column] = result_df[column] * np.log(result_df[column])
        result_df["entropy"] = result_df.sum(axis=1, numeric_only=True)
        # Multiply by -1
        result_df["entropy"] = result_df["entropy"]*-1
        # Number of effective states
        result_df["n_effective"] = np.exp(result_df["entropy"])
        # Reset site values
        result_df["site"] = sites

    # Return result
    return result_df

def main():
    """
    Main method
    """

    # Input files
    protein_alignment_file = str(snakemake.input.protein_alignment)

    # Output files
    output_file = str(snakemake.output)

    alignment_counts_df = get_counts_from_multiple_sequence_alignment(protein_alignment_file)
    sequence_variability_df = calculate_entropy_and_neffective_per_site(alignment_counts_df)
    # Write output file
    sequence_variability_df.to_csv(output_file, index=False)



if __name__ == "__main__":
    main()