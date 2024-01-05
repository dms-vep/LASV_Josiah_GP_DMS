
# Description:
# Create a fasta sequence from the consensus sequence
# in the mpileup file.

# Author:
# Caleb Carr

# Imports
import pandas as pd 

# This function is used to prepare a pandas 
# column to write to the output file
def pandas_column_to_string(column):
    """ Convert pandas column to string """
    result = ""
    for index in range(len(column)):
        result += column.iloc[index]
    return result

# Read input csv file that has alignment data
alignment_data = pd.read_csv(str(snakemake.input))

# Open/create the output file accessed through the snakemake 
# object and allow it to be written to
output_file = open(str(snakemake.output), "w")

# Extract all aligned gene reference names
gene_list = snakemake.params.gene_list

# Iterate through all genes
for gene_name in gene_list:

    # Extract data associated with gene
    gene_data = alignment_data.loc[alignment_data['Reference'] == gene_name]

    # Write gene name with '>' symbol for fasta format
    output_file.write(">"+gene_name+"\n")

    # Extract consensus bases for given gene
    consensus_sample_bases = gene_data['Consensus Base']

    # Check consensus sequence is not empty
    if len(consensus_sample_bases) == 0:
        # Write single N as place holder if no consensus was created
        output_file.write("N"+"\n")
    else:
        # Iterate in chuncks of 70 to write the bases to the output file.
        # Chunks of 70 were chosen because that is the width of each fasta file.
        for index in range(0, len(consensus_sample_bases), 70):
            # Convert 70 base chunk of column to string
            new_line = pandas_column_to_string(consensus_sample_bases.iloc[index:index+70])
            # Write string to ouput file
            output_file.write(new_line+"\n")


# Close output files
output_file.close()