# Description:
# Python script to download all fasta sequences
# specified by a list of accession numbers as well
# as metadata

# Author:
# Caleb Carr

# Imports
import datetime
import os
import pandas as pd
from Bio import Entrez, SeqIO, AlignIO
Entrez.email = "example@example.com"     # Always tell NCBI who you are

# Functions
def read_and_process_accession_list(
    input_file_name, 
    output_fasta_file_name, 
    output_metadata_file_name, 
    length_threshold, 
    desired_segment
    ):
    """
    Function to read in list of accessions, download genbank files,
    parse genbank files, and extract sequences/metadata
    """

    def country_extraction(location):
        """
        Function to extract country
        """
        # Upper case location
        location = location.upper()

        # Extract country
        if "NIGERIA" in location:
            return "Nigeria"
        elif "IVOIRE" in location:
            return "Cote d'Ivoire"
        elif "LIBERIA" in location:
            return "Liberia"
        elif "SIERRA" in location:
            return "Sierra Leone"
        elif "GUINEA" in location:
            return "Guinea"
        elif "MALI" in location:
            return "Mali"
        elif "BENIN" in location:
            return "Benin"
        elif "TOGO" in location:
            return "Togo"
        else:
            return location

    # Open output files
    output_fasta_file = open(output_fasta_file_name, "w")
    output_metadata_file = open(output_metadata_file_name, "w")
    total_accessions_count = 0
    removed_accessions_count = 0

    # Write header for metadata file
    header = [
        "strain",
        "virus",
        "segment",
        "host",
        "accession",
        "date",
        # For now, just using a single location field
        # "region",
        # "country",
        # "division",
        # "city",
        "location",
        "country",
        "database",
        "authors",
        "url",
        "title",
        "journal",
        "paper_url",
    ]
    header = "\t".join(header) + "\n"
    output_metadata_file.write(header)
    
    # Open file with list of accession numbers
    with open(input_file_name, "r") as input_file:
        for line in input_file:

            # Initialize results for metadata
            strain = "MISSING"
            virus = "MISSING"
            segment = "MISSING"
            host = "MISSING"
            accession = "MISSING"
            date = "MISSING"
            # For now, just using a single location field
            # region = ""
            # country = ""
            # division = ""
            # city = ""
            location = "MISSING"
            country = "MISSING"
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

            # Extract current accession ID
            accession_from_list = line.split()[0]
            total_accessions_count += 1

            # Check if accession is to be excluded
            if accession_from_list[:-2] in snakemake.params.accesstions_to_exclude:
                print(f"{accession_from_list} excluded based on config file!\n")
                removed_accessions_count += 1
                continue

            # Retrieve genbank file for accession ID
            entrez_genbank = Entrez.efetch(
                db="nucleotide", 
                id=accession_from_list, 
                rettype="genbank", 
                retmode="text"
                )
            # Parse genbank file line by line to retrieve all metadata
            for line in entrez_genbank:

                # Process line by removing spaces
                line = " ".join([ele for ele in line.split(" ") if ele != ""])
                split_line = line.split(" ")

                # Extract feature information from genbank file
                if features_flag == True:
                    if "/isolate" in line or "/strain" in line:
                        strain = line.replace("\n", "").replace("\"", "").split("=")[1]
                    if "/organism" in line:
                        virus = line.replace("\n", "").replace("\"", "").split("=")[1]
                    if "/host" in line:
                        host = line.replace("\n", "").replace("\"", "").split("=")[1]
                    if "/segment" in line:
                        segment = line.replace("\n", "").replace("\"", "").split("=")[1]
                    if "/country" in line:
                        # For now, just using a single location field
                        location = line.replace("\n", "").replace("\"", "").split("=")[1]
                        country = country_extraction(location)
                    if "/collection_date" in line:
                        unformatted_date = line.replace("\n", "").replace("\"", "").split("=")[1]
                        if len(unformatted_date.split("-")) == 1:
                            date = datetime.datetime.strptime(unformatted_date, "%Y").strftime("%Y") + "-XX-XX"
                        elif len(unformatted_date.split("-")) == 2:
                            date = datetime.datetime.strptime(unformatted_date, "%b-%Y").strftime("%Y-%m") + "-XX"
                        elif len(unformatted_date.split("-")) == 3:
                            test_res = True
                            try:
                                test_res = bool(datetime.datetime.strptime(unformatted_date, "%Y-%m-%d"))
                            except ValueError:
                                test_res = False
                            if test_res:
                                date = unformatted_date
                            else:
                                date = datetime.datetime.strptime(unformatted_date, "%d-%b-%Y").strftime("%Y-%m-%d")
                        else:
                            print("Datetime not between 1 and 3")
                
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
                    accession = split_line[1].replace("\n", "")
                    genbank_base_url = "https://www.ncbi.nlm.nih.gov/nuccore/"
                    url = genbank_base_url + accession
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
                    journal = " ".join(split_line[1:]).replace("\n", "")
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

            # Check if only specific segments should be kept
            if desired_segment != None and segment != desired_segment:
                continue

            # If not collection date is found, add date from top of genbank file
            if date == "MISSING":
                date = backup_date

            # Retrieve genbank file for accession ID
            entrez_efetch = Entrez.efetch(db="nucleotide", id=accession_from_list, rettype="fasta", retmode="text")
            # Convert current genbank accession to fasta SeqIO format
            fasta = SeqIO.read(entrez_efetch, "fasta")

            # Check if sequence is below ambiguous base threshold
            if str(fasta.seq).upper().count("N")/length > snakemake.params.max_frac_N:
                print(f"{accession_from_list} excluded because of high N count!\n")
                removed_accessions_count += 1
                continue
            
            # Calculate the best matching orientation of sequence based on alignement to reference
            best_match_name, best_match_sequence = align_sequences_and_remove_low_matches(fasta.seq)
            
            # Join strain/isolate name with accession and date to make sure it is unique
            if best_match_name == "not_adjusted":
                strain = strain + "_" + accession + "_" + date
            else:
                strain = strain + "_" + accession + "_" + best_match_name + "_" + date

            # Replace slashes, periods, and spaces in name with underscores
            strain = strain.replace("/", "_")
            strain = strain.replace(". ", "-")
            strain = strain.replace(" ", "_")
            strain = strain.replace(".", "-")

            # Make sure every sequence has a fasta strain name
            assert strain != "MISSING", "Virus strain name is missing"

            # Create new metadata line
            new_metadata_line = "\t".join([
                strain,
                virus,
                segment,
                host,
                accession,
                date,
                # For now, just using a single location field
                # region,
                # country,
                # division,
                # city,
                location,
                country,
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
            output_fasta_file.write(f">{strain}\n")
            output_fasta_file.write(f"{best_match_sequence}\n")

    print(f"A total of {total_accessions_count} were processed and ")
    print(f"{total_accessions_count-removed_accessions_count} were retained!\n")   
    # Close files
    input_file.close()
    output_fasta_file.close()
    output_metadata_file.close()


def align_sequences_and_remove_low_matches(sequence):
    """
    Function to align a sequence, reverse of sequence, and 
    reverse complement of a sequence and calculate the best match
    after aligning to reference.  
    """

    # Intialize temp df for comparisons of sequence, reverse, and reverse complement
    temp_df = pd.DataFrame(columns=["Name", "Seq", "PID"])
    S_segment_reference = str(snakemake.input.reference_genome) 
    # Temp files to use for alignment of the three orientations to reference
    temp_fasta_file = str(snakemake.params.temp_fasta_file) 
    temp_alignment_file = str(snakemake.params.temp_alignment_file) 

    # Add sequence not adjusted, reverse, and reverse complement
    # Not adjusted
    temp_df.loc[len(temp_df.index)] = ["not_adjusted", sequence, 0] 
    # Reverse
    sequence = sequence[::-1] 
    temp_df.loc[len(temp_df.index)] = ["reverse", sequence, 0] 
    # Reverse complement
    sequence = sequence.upper()
    sequence = (
        sequence
        .replace("A", "t")
        .replace("C", "g")
        .replace("T", "a")
        .replace("G", "c")
        .replace("N", "n")
    )
    sequence = sequence.upper()
    temp_df.loc[len(temp_df.index)] = ["reverse_complement", sequence, 0] 

    # Open output file
    temp_output_file = open(temp_fasta_file, "w")

    # Write the three sequences to a temp file
    for i in range(len(temp_df)):
        temp_output_file.write(f">{temp_df.at[i, 'Name']}\n")
        temp_output_file.write(f"{temp_df.at[i, 'Seq']}\n")

    # Close files
    temp_output_file.close()

    # Align not adjusted, reverse, and reverse complement
    os.system(f"mafft --6merpair --keeplength --quiet --addfragments {temp_fasta_file} {S_segment_reference} > {temp_alignment_file}")

    # Calculate percent identities for each aligned sequence
    temp_alignment = AlignIO.read(temp_alignment_file, "fasta")
    pid_list = []
    for seq in temp_alignment:
        pid_list.append(percent_ids(temp_alignment[0], seq.seq))
    temp_df["PID"] = pid_list[1:]

    best_match_name = temp_df.at[temp_df["PID"].idxmax(), "Name"]
    best_match_sequence = temp_df.at[temp_df["PID"].idxmax(), "Seq"]

    # Return name of sequence orientation and sequence
    return best_match_name, best_match_sequence


def percent_ids(seq1, seq2):
    """
    Function to calculate percent similarity between two sequences
    """
    # Make sure aligned sequences have same length
    assert len(seq1) == len(seq2), "Aligned sequences do not have same length!"

    # Count base similarities 
    length = len(seq1)
    num_gaps = seq1.count('-')
    num_bases = length - num_gaps
    matching_bases = 0
    for i in range(length):
        if seq1[i] == seq2[i] and seq1[i] != "-":
            matching_bases += 1

    # Return fraction similarity between the two sequences
    return (matching_bases)/num_bases


def main():
    """
    Main method
    """

    # Input files
    list_of_accessions = str(snakemake.input.accession_list) # "../Data/all_LASV_accessions_041423.txt"
    # Params
    length_threshold = (
        int(str(snakemake.params.genome_size_threshold_lower)), 
        int(str(snakemake.params.genome_size_threshold_upper))
    )
    # Empty String to None Conversion
    None_conversion = lambda i : i or None
    desired_segment = None_conversion(str(snakemake.params.desired_segment))

    # Output files
    fasta_output = str(snakemake.output.fasta_sequences) # "../Data/all_LASV_fasta_041423.fasta"
    metadata_output = str(snakemake.output.metadata) # "../Data/all_LASV_metadata_041423.tsv"

    read_and_process_accession_list(
        list_of_accessions, 
        fasta_output, 
        metadata_output, 
        length_threshold, 
        desired_segment
        )


if __name__ == "__main__":
    main()