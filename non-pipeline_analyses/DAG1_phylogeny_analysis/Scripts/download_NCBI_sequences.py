# Description:
# Python script to download all fasta sequences
# specified by a list of accession numbers as well
# as metadata

# Author:
# Caleb Carr

# Imports
import datetime
import pandas as pd
from Bio import Entrez, SeqIO
Entrez.email = "example@example.com"     # Always tell NCBI who you are

# Functions
def read_and_process_accession_list(
        list_of_mRNA_accessions, 
        list_of_protein_accessions, 
        CDS_fasta_output_file_name, 
        protein_fasta_output_file_name, 
        output_metadata_file_name, 
        length_threshold
    ):
    """
    Function to read in list of accessions, download genbank files,
    parse genbank files, and extract sequences/metadata
    """

    # Open output files
    output_CDS_fasta_file = open(CDS_fasta_output_file_name, "w")
    output_protein_fasta_file = open(protein_fasta_output_file_name, "w")
    output_metadata_file = open(output_metadata_file_name, "w")

    # Write header for metadata file
    header = [
        "name",
        "molecule_type",
        "chromosome",
        "organism",
        "mRNA_accession",
        "number_mRNA_seqs",
        "protein_accession",
        "date",
        "database",
        "authors",
        "mRNA_url",
        "title",
        "journal",
        "paper_url",
    ]
    header = "\t".join(header) + "\n"
    output_metadata_file.write(header)
    
    # Open file with list of accession numbers
    with open(list_of_mRNA_accessions, "r") as mRNA_accessions, open(list_of_protein_accessions, "r") as protein_accessions:
        for mRNA_line, protein_line in zip(mRNA_accessions, protein_accessions):

            if mRNA_line == "" or protein_line == "":
                continue

            # Initialize results for metadata
            name = "MISSING"
            molecule_type = "MISSING"
            chromosome = "MISSING"
            organism = "MISSING"
            mRNA_accession = "MISSING"
            protein_accession = "MISSING"
            date = "MISSING"
            database = "MISSING"
            authors = "MISSING"
            url = "MISSING"
            title = "MISSING"
            journal = "MISSING"
            paper_url = "MISSING"
            features_flag = False
            reference_count = 0
            length = 0
            backup_date = "MISSING"
            CDS = [0,0] # assumed that the mRNA has only one CDS and no introns

            # Extract protein accession ID
            protein_accession_from_list = protein_line.split()[0]
            protein_accession = protein_accession_from_list.split(".")[0]

            # Extract mRNA accession ID
            accession_from_list = mRNA_line.split()[0]

            print(f"{accession_from_list} {protein_accession_from_list}")

            # Retrieve genbank file for accession ID
            entrez_genbank = Entrez.efetch(db="nucleotide", id=accession_from_list, rettype="genbank", retmode="text")
            # Parse genbank file line by line to retrieve all metadata
            for line in entrez_genbank:

                # Process line by removing spaces
                line = " ".join([ele for ele in line.split(" ") if ele != ""])
                split_line = line.split(" ")

                # Extract feature information from genbank file
                if features_flag == True:
                    if "/organism" in line:
                        organism = line.replace("\n", "").replace("\"", "").split("=")[1].replace(" ", "-")
                    if "/mol_type" in line:
                        molecule_type = line.replace("\n", "").replace("\"", "").split("=")[1]
                    if "/chromosome" in line:
                        chromosome = line.replace("\n", "").replace("\"", "").split("=")[1]
                    if "CDS " in line:
                        CDS_line = line.replace("\n", "").replace("CDS", "").replace("<", "").replace(">", "").split("..")
                        CDS[0] = int(CDS_line[0]) - 1 # Genbank is 1 indexed and python is 0 indexed
                        CDS[1] = int(CDS_line[1])

                
                if split_line[0] == "LOCUS":
                    unformatted_backup_date = split_line[-1].replace("\n", "")
                    if len(unformatted_backup_date.split("-")) == 1:
                        backup_date = datetime.datetime.strptime(unformatted_backup_date, "%Y").strftime("%Y") + "-XX-XX"
                    elif len(unformatted_backup_date.split("-")) == 2:
                        backup_date = datetime.datetime.strptime(unformatted_backup_date, "%b-%Y").strftime("%Y-%m") + "-XX"
                    elif len(unformatted_backup_date.split("-")) == 3:
                        test_res = True
                        try:
                            test_res = bool(datetime.datetime.strptime(unformatted_backup_date, "%Y-%m-%d"))
                        except ValueError:
                            test_res = False
                        if test_res:
                            backup_date = unformatted_backup_date
                        else:
                            backup_date = datetime.datetime.strptime(unformatted_backup_date, "%d-%b-%Y").strftime("%Y-%m-%d")
                    else:
                        print("Datetime not between 1 and 3")
                    length = int(split_line[2])
                    continue
                if split_line[0] == "ACCESSION":
                    mRNA_accession = split_line[1].replace("\n", "")
                    genbank_base_url = "https://www.ncbi.nlm.nih.gov/nuccore/NC_004297"
                    url = genbank_base_url + mRNA_accession
                    continue
                if split_line[0] == "DBLINK":
                    database = " ".join(split_line[1:]).replace("\n", "")
                    continue
                if split_line[0] == "REFERENCE":
                    reference_count += 1
                    continue
                if split_line[0] == "AUTHORS" and reference_count == 1:
                    authors = split_line[1].split(",")[0] + " et al"
                    continue
                if split_line[0] == "TITLE" and reference_count == 1:
                    title = " ".join(split_line[1:]).replace("\n", "")
                    continue
                if split_line[0] == "JOURNAL" and reference_count == 1:
                    journal = split_line[1].replace("\n", "")
                    continue
                if split_line[0] == "PUBMED" and reference_count == 1:
                    pubmed_base_url = "https://pubmed.ncbi.nlm.nih.gov/"
                    paper_url = pubmed_base_url + split_line[1].replace("\n", "")
                    continue
                if split_line[0] == "FEATURES":
                    features_flag = True # Set features flag to then extract feature information
                    continue   

            # Check length of sequence and do not add if below threshold
            if length < length_threshold[0] or length > length_threshold[1]:
                continue

            # If not collection date is found, add date from top of genbank file
            if date == "MISSING":
                date = backup_date

            # Retrieve genbank file for accession ID
            entrez_efetch = Entrez.efetch(db="nucleotide", id=accession_from_list, rettype="fasta", retmode="text")
            # Convert current genbank accession to fasta SeqIO format
            fasta = SeqIO.read(entrez_efetch, "fasta")
            CDS_fasta = fasta.seq[CDS[0] : CDS[1]]
            
            # Join name with accession and date to make sure it is unique
            name = organism + "_" + mRNA_accession + "_" + protein_accession + "_" + date

            # Make sure every sequence has a fasta strain name
            assert name != "MISSING", "Name is missing"

            # Create new metadata line
            new_metadata_line = "\t".join([
                name,
                molecule_type,
                chromosome,
                organism,
                mRNA_accession,
                protein_accession,
                date,
                database,
                authors,
                url,
                title,
                journal,
                paper_url,
            ])
            new_metadata_line += "\n"

            # Write new metadata line
            output_metadata_file.write(new_metadata_line)

            # Write current fasta sequence to output file
            output_CDS_fasta_file.write(f">{name}\n")
            output_CDS_fasta_file.write(f"{CDS_fasta}\n")

            # Retrieve protein fasta and convert to SeqIO
            protein_entrez_efetch = Entrez.efetch(db="protein", id=protein_accession_from_list, rettype="fasta", retmode="text")
            protein_fasta = SeqIO.read(protein_entrez_efetch, "fasta")

            # Write current fasta sequence to output file
            output_protein_fasta_file.write(f">{name}\n")
            output_protein_fasta_file.write(f"{protein_fasta.seq}\n")
                
    # Close files
    mRNA_accessions.close()
    protein_accessions.close()
    output_CDS_fasta_file.close()
    output_protein_fasta_file.close()
    output_metadata_file.close()


def main():
    """
    Main method
    """

    # Input files
    list_of_mRNA_accessions = str(snakemake.input.mRNA_accession_list)
    list_of_protein_accessions = str(snakemake.input.protein_accession_list)
    length_threshold = (
        int(str(snakemake.params.genome_size_threshold_lower)), 
        int(str(snakemake.params.genome_size_threshold_upper))
    )

    # Output files
    CDS_fasta_output = str(snakemake.output.CDS_fasta_sequences)
    protein_fasta_output = str(snakemake.output.protein_fasta_sequences)
    metadata_output = str(snakemake.output.metadata)

    # Run function to download and process list of accessions
    read_and_process_accession_list(
        list_of_mRNA_accessions, 
        list_of_protein_accessions, 
        CDS_fasta_output, 
        protein_fasta_output, 
        metadata_output, 
        length_threshold,
    )


if __name__ == "__main__":
    main()