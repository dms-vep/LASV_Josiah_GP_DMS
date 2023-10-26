# Description:
# Python script to loop through the names of sequences
# and corresponding ORFs and verify only one ORF was 
# detected for each sequence. 

# Author:
# Caleb Carr

# Imports
from itertools import zip_longest

# Functions
def assert_num_orfs_match_sequences(
    segment_fasta_file, 
    protein_fasta_file, 
    codon_fasta_file, 
    log_output_file
    ):

    # Initialize lists for sequence names and flag
    segment_list = []
    protein_list = []
    codon_list = []
    verification_flag = True
    sequence_with_multiple_ORFs = ""

    with (open(segment_fasta_file, "r") as segment_fasta, 
        open(protein_fasta_file, "r") as protein_fasta, 
        open(codon_fasta_file, "r") as codon_fasta):

        for (segment_line, protein_line, codon_line) in zip_longest(segment_fasta, protein_fasta, codon_fasta):

            # Extract sequence names
            if segment_line != None and segment_line[0] == ">":
                segment_list.append(segment_line[1:-1])
            if protein_line != None and protein_line[0] == ">":
                protein_list.append(protein_line[1:].split(" ")[0][:-2])
            if codon_line != None and codon_line[0] == ">":
                codon_list.append(codon_line[1:].split(" ")[0][:-2])

    # Check to make sure equal number of codon and protein sequences detected
    assert len(protein_list) == len(codon_list), "Uneven protein and nucleotide ORFs detected!"

    # Loop through lists of names to check there is a one-to-one mapping
    segment_index = 0
    orf_index = 0
    checked_sequences = []
    for i in range(0, len(segment_list)):
        segment_name = segment_list[segment_index]
        protein_name = protein_list[orf_index]
        codon_name = codon_list[orf_index]
        # Check if all have same name and set flag to false if not
        if segment_name == protein_name == codon_name:
            checked_sequences.append(segment_name)
            segment_index += 1
            orf_index += 1
            continue
        # Check if curr ORF has not been seen yet and update segment index
        elif protein_name not in checked_sequences:
            segment_index += 1
            continue
        # Check if curr ORF has been seen and flag multiple ORFs detected
        elif protein_name in checked_sequences:
            verification_flag = False
            sequence_with_multiple_ORFs = protein_name
            break
        # Catch all issue
        else:
            verification_flag = False
            sequence_with_multiple_ORFs = segment_name
            print("Unknown error occured!")
            break

    # Check if flag was triggered
    if verification_flag:
        # Open output files
        output_log_file = open(log_output_file, "w")
        output_log_file.write("All sequences confirmed to have one-to-one ORF mappings for GPC segment!\n")
        output_log_file.write(f"From {len(segment_list)} S segment sequences, {len(protein_list)} GPC ORFs found!")
        # Close files
        output_log_file.close()
    else:
        print(f"{sequence_with_multiple_ORFs} has more than one ORF in GPC!")
            

    # Close files
    segment_fasta.close()
    protein_fasta.close()
    codon_fasta.close()




def main():
    """
    Main method
    """

    # Input files
    segment_fasta_file = str(snakemake.input.nucleotide)
    protein_fasta_file = str(snakemake.input.protein)
    codon_fasta_file = str(snakemake.input.codon)

    # Output files
    log_output_file = str(snakemake.output)

    assert_num_orfs_match_sequences(
        segment_fasta_file,
        protein_fasta_file,
        codon_fasta_file,
        log_output_file
    )


if __name__ == "__main__":
    main()