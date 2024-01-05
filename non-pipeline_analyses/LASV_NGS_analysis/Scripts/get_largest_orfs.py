
# Description:
# Parses the output of running the EMBOSS command
# 'getorf' to find the largest (most likely) protein
# sequences for a given fasta file. Returns a file
# with the complete fasta sequences as well as a 
# short file with just the indices for each coding
# region for each gene.

# Author:
# Caleb Carr


# Open the input file accessed through the snakemake object
# and only read lines from it
first_pass_orf_file = open(str(snakemake.input), "r")

# Open/create the output files accessed through the snakemake 
# object and allow it to be written to
output_file = open(str(snakemake.output.sequences), "w")
indicies_output_file = open(str(snakemake.output.positions), "w")

# Create header for indices output file
indicies_output_file.write("Reference,Start CDS,End CDS\n")

# Create a dictionary to hold the indices and sizes of the 
# largest protein sequences found in the file for each gene.
# The index is which line in the file the protein sequence
# is found which is used to extract the sequence later.
largest_orfs = {}

# Iterate through the input file
for index, line in enumerate(first_pass_orf_file):

    # Split each line based on white space
    sectioned = line.split()

    # Lines with gene names start with '>' symbols
    if sectioned[0][0] == '>':
        
        # Extract the gene name and indices designating
        # size of the protein sequence and calculate 
        # protein sequence.
        gene_name = sectioned[0].split('_')[0] 
        first_index = int(sectioned[1][1:])
        second_index = int(sectioned[3][:-1])
        curr_orf_size = second_index - first_index

        # Based on the indices of the size of protein sequence,
        # positive orf sizes correspond to forward strand sequences
        if curr_orf_size > 0:

            # Check if current gene is already in the dictionary
            if gene_name in largest_orfs:
                # If the gene is in the dictionary and size is 
                # larger, then update the size associated with 
                # the gene and index
                if curr_orf_size > largest_orfs[gene_name][1]:
                    # Update dictionary with new size and index
                    largest_orfs[gene_name] = [index, curr_orf_size, first_index, second_index]
            else:
                # Add gene and current size to dictionary if it
                # hasn't been encountered before
                largest_orfs[gene_name] = [index, curr_orf_size, first_index, second_index]


# Iterate through dictionary and put contents in output file
for key in largest_orfs:
    # Add gene name
    reference = key[1:] # Slice off the '>' symbol from the fasta format
    # Add start and end indices
    start = largest_orfs[key][2]
    end = largest_orfs[key][3]
    # Write new line to indices output file
    indicies_output_file.write(f"{reference},{start},{end}\n")


# Open the input file accessed through the snakemake object
# and only read lines from it for a second pass through file
second_pass_orf_file = open(str(snakemake.input), "r")

# Set a boolean to signal when the line correspond to 
# protein sequence
coding_sequence = False
# Iterate through input file
for index, line in enumerate(second_pass_orf_file):

    # Split each line based on white space
    sectioned = line.split() 
    # Extract gene name
    gene_name = sectioned[0].split('_')[0]

    # If first character is '>' and coding_sequence is true
    # then the current protein sequence was the above lines
    # before and the boolean needs to be reset to false
    if sectioned[0][0] == '>' and coding_sequence == True:
        # End gene segement with endline
        output_file.write("\n")
        coding_sequence = False
    
    # If coding_sequence is true then the current lines are 
    # the protein sequence so write these lines to the output 
    # file
    if coding_sequence == True:
        output_file.write(line)

    # If first character is '>' and index is in the dictionary of 
    # largest orfs, then start writing the protein sequence to the
    # output file
    if sectioned[0][0] == '>' and index == largest_orfs[gene_name][0]:
        output_file.write(f"{gene_name} {largest_orfs[gene_name][2]}-{largest_orfs[gene_name][3]}\n")
        coding_sequence = True
    

# Close input and output files
first_pass_orf_file.close()
second_pass_orf_file.close()
indicies_output_file.close()
output_file.close()