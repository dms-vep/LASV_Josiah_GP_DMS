# Description:
# Extracts all accessions with matching mRNA and protein sequences

# Author:
# Caleb Carr

# Imports
import pandas as pd
import numpy as np


def extract_accessions(accession_csv, mRNA_accession_output, protein_accession_output):
    """
    Function to check and extract all accessions that have matching
    mRNA and protein accessions
    """

    # Read csv as dataframe and remove and rows that do not have a protein accession
    accession_df = pd.read_csv(accession_csv)
    accession_df["RefSeq Protein accessions"] = accession_df["RefSeq Protein accessions"].replace("", np.nan)
    accession_df = accession_df.dropna(subset=["RefSeq Protein accessions"])

    # Write accession columns to output files
    accession_df["RefSeq Transcript accessions"].to_csv(mRNA_accession_output, header=None, index=None)
    accession_df["RefSeq Protein accessions"].to_csv(protein_accession_output, header=None, index=None)



def main():
    """
    Main method
    """

    # Input files
    accession_csv = str(snakemake.input)

    # Output files
    mRNA_accession_output = str(snakemake.output.mRNA_accession_list)
    protein_accession_output = str(snakemake.output.protein_accession_list)

    # Run function to remove accession without a matching protein accession
    extract_accessions(accession_csv, mRNA_accession_output, protein_accession_output)

if __name__ == "__main__":
    main()