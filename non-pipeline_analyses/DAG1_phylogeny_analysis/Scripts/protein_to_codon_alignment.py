# Description:
# Python script to convert protein alignment 
# to codon alignment

# Author:
# Caleb Carr

# Imports
from Bio import Entrez, SeqIO
Entrez.email = "example@example.com" 

def convert_protein_to_codon(input_protein_alignment_file, input_CDS_fastas_file, output_codon_alignment_file):
    """
    Function to process protein alignment and
    convert to codon alignment given CDS fastas
    """

    # Open input and output files
    protein_alignment_file = open(input_protein_alignment_file, "r") 
    CDS_fasta_file = open(input_CDS_fastas_file, "r")
    codon_alignment_file = open(output_codon_alignment_file, "w")

    # Parse alignment and fasta files
    protein_alignment = SeqIO.parse(protein_alignment_file, "fasta")
    CDS_fastas = SeqIO.parse(CDS_fasta_file, "fasta")

    # Iterate through each alignment and corresponding CDS fasta
    for alignment, CDS in zip(protein_alignment, CDS_fastas):
        
        # Check to make sure alignment and CDS names match
        assert alignment.id == CDS.id, "Alignment and CDS names do not match!"

        # Initialize new codon alignment
        curr_name = alignment.id
        curr_codon_alignment = ""
        index = 0 # Index to keep track of codon position
        
        # Iterate through protein alignment and add corresponding codon
        for AA in alignment.seq:

            codon = str(CDS.seq[index:index+3])
            # Pad codon with gaps if not length of 3
            if len(codon) != 3:
                codon = codon.ljust(3, "-")

            # Add codon if not amino acid is not a gap otherwise add three gaps
            if AA != "-":
                curr_codon_alignment += codon
                index += 3
            else:
                curr_codon_alignment += "---"


        # Check to make sure codon alignment equals protein alignment*3
        assert len(curr_codon_alignment) == len(alignment.seq)*3, "Alignment and codon lengths do not match!"

        # Write current fasta sequence to output file
        codon_alignment_file.write(f">{curr_name}\n")
        codon_alignment_file.write(f"{curr_codon_alignment}\n")

    # Close files
    protein_alignment_file.close()
    CDS_fasta_file.close()
    codon_alignment_file.close()


def main():
    """
    Main method
    """

    # Input files
    protein_alignment = str(snakemake.input.protein_alignment)
    CDS_fastas = str(snakemake.input.CDS_fastas)

    # Output files
    codon_alignment = str(snakemake.output.codon_alignment)

    convert_protein_to_codon(protein_alignment, CDS_fastas, codon_alignment)


if __name__ == "__main__":
    main()