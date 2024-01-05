
# Description:
# Checks the best matched reference genome and then aligns to that 
# reference genome for each sample.

# Author:
# Caleb Carr


# Open the input file accessed through the snakemake object
# and only read lines from it
mpileup_file = open(str(snakemake.input), "r")

# Open/create the output file accessed through the snakemake 
# object and allow it to be written to
output_file = open(str(snakemake.output), "w")

# Create header for output file
output_file.write(("Reference,Depth,Position,Reference Base,Consensus Base,"
                    "A For,A Rev,C For,C Rev,G For,G Rev,T For,T Rev,"
                    "Insertions,Deletions,"
                    "InsertsLTE10,InsertsLTE50,InsertsLTE100,InsertsLTE500,InsertsLTE1000,InsertsLTE1500,InsertsGT1500,"
                    "DelsLTE10,DelsLTE50,DelsLTE100,DelsLTE500,DelsLTE1000,DelsLTE1500,DelsGT1500\n"))

# Iterate through each line of the mpileup file
for line in mpileup_file:
    # Sections of the mpileup file:
    # -First section is the name of mapped region
    # -Second section is the position in the reference
    # -Third section is the reference genome base
    # -Fourth section is the read bases

    # Split each line of the mpileup file based on white space
    sectioned = line.split()
    # Reference gene name
    reference = sectioned[0]
    # Position in gene
    position = sectioned[1]
    # Reference base
    reference_base = sectioned[2].upper()
    
    # total indel counts
    insertions = 0 # counts for insertions
    deletions = 0 # counts for deletions
    # base counts
    A_count_for = 0 # counts for adenine on forward strand
    T_count_for = 0 # counts for thymine on forward strand
    C_count_for = 0 # counts for cytosine on forward strand
    G_count_for = 0 # counts for guanine on forward strand
    A_count_rev = 0 # counts for adenine on reverse strand
    T_count_rev = 0 # counts for thymine on reverse strand
    C_count_rev = 0 # counts for cytosine on reverse strand
    G_count_rev = 0 # counts for guanine on reverse strand
    # counts for the specified indel size bin
    insert_10 = 0 # counts of insert size 10 or less
    insert_50 = 0 # counts of insert size 50 to 11
    insert_100 = 0 # counts of insert size 100 to 51
    insert_500 = 0 # counts of insert size 500 to 101
    insert_1000 = 0 # counts of insert size 1000 to 501
    insert_1500 = 0 # countsof insert size 1500 to 1001
    insert_GT_1500 = 0 # counts for indels > 1500
    dels_10 = 0 # counts of dels size 10 or less
    dels_50 = 0 # counts of dels size 50 to 11
    dels_100 = 0 # counts of dels size 100 to 51
    dels_500 = 0 # counts of dels size 500 to 101
    dels_1000 = 0 # counts of dels size 1000 to 501
    dels_1500 = 0 # counts of dels size 1500 to 1001
    dels_GT_1500 = 0 # counts for indels > 1500
    # total depth for each position
    total = 0 # total count of bases and indels (read depth)

    # Read base data for current position
    base_characters = sectioned[4]
    index = 0 # Index for iterating through all base data positions
    # Iterate through all base data positions
    while index < len(base_characters):

        if base_characters[index] == "^":
            # Skip next character because '^' signals that the position
            # is the start of a read and the following character is the 
            # read quality score.
            index +=2
            continue

        if base_characters[index] == "+":
            insertions +=1
            # Skip next number of characters because '+' or '-' signals that
            # the position has an insertion or deletion and the following 
            # digits are the size of the indel.
            num_index = index + 1
            indel_size = ""
            while base_characters[num_index].isdigit():
                indel_size = indel_size + base_characters[num_index]
                num_index +=1
            # Increment indel bin size count
            if int(indel_size) <= 10:
                insert_10 += 1
            elif int(indel_size) <= 50:
                insert_50 += 1
            elif int(indel_size) <= 100:
                insert_100 += 1
            elif int(indel_size) <= 500:
                insert_500 += 1
            elif int(indel_size) <= 1000:
                insert_1000 += 1
            elif int(indel_size) <= 1500:
                insert_1500 += 1
            else:
                insert_GT_1500 += 1
            # Increment the index to skip the indel
            index = index + 2 + int(indel_size)
            continue
        
        if base_characters[index] == "-":
            deletions += 1
            # Skip next number of characters because '+' or '-' signals that
            # the position has an insertion or deletion and the following 
            # character is the size of the indel.
            num_index = index + 1
            indel_size = ""
            while base_characters[num_index].isdigit():
                indel_size = indel_size + base_characters[num_index]
                num_index +=1
            # Increment indel bin size count
            if int(indel_size) <= 10:
                dels_10 += 1
            elif int(indel_size) <= 50:
                dels_50 += 1
            elif int(indel_size) <= 100:
                dels_100 += 1
            elif int(indel_size) <= 500:
                dels_500 += 1
            elif int(indel_size) <= 1000:
                dels_1000 += 1
            elif int(indel_size) <= 1500:
                dels_1500 += 1
            else:
                dels_GT_1500 += 1
            # Increment the index to skip the indel
            index = index + 2 + int(indel_size)
            continue

        if base_characters[index] == "*":
            # Increment index if a '*' is found because that signals a 
            # deleted base.
            index +=1
            continue

        if base_characters[index] == ".":
            # Increment the corresponding base count if a '.' is 
            # found because that corresponds to a correct match with the
            # reference genome on the forward strand.
            if reference_base == "A":
                A_count_for +=1
            if reference_base == "C":
                C_count_for +=1
            if reference_base == "G":
                G_count_for +=1
            if reference_base == "T":
                T_count_for +=1
            index +=1
            continue

        if base_characters[index] == ",":
            # Increment the corresponding base count if a ',' is 
            # found because that corresponds to a correct match with the
            # reference genome on the reverse strand.
            if reference_base == "A":
                A_count_rev +=1
            if reference_base == "C":
                C_count_rev +=1
            if reference_base == "G":
                G_count_rev +=1
            if reference_base == "T":
                T_count_rev +=1
            index +=1
            continue
        
        #################---FORWARD STRAND---#################
        # Increment adenine count forward if 'A' is encountered
        if base_characters[index] == "A":
            A_count_for +=1
            index +=1
            continue
        # Increment cytosine count forward if 'C' is encountered
        if base_characters[index] == "C":
            C_count_for +=1
            index +=1
            continue
        # Increment guanine count forward if 'G' is encountered
        if base_characters[index] == "G":
            G_count_for +=1
            index +=1
            continue
        # Increment thymine count forward if 'T' is encountered
        if base_characters[index] == "T":
            T_count_for +=1
            index +=1
            continue
        #################---FORWARD STRAND---#################

        #################---REVERSE STRAND---#################
        # Increment adenine count reverse if 'a' is encountered
        if base_characters[index] == "a":
            A_count_rev +=1
            index +=1
            continue
        # Increment cytosine count reverse if 'c' is encountered
        if base_characters[index] == "c":
            C_count_rev +=1
            index +=1
            continue
        # Increment guanine count reverse if 'g' is encountered
        if base_characters[index] == "g":
            G_count_rev +=1
            index +=1
            continue
        # Increment thymine count reverse if 't' is encountered
        if base_characters[index] == "t":
            T_count_rev +=1
            index +=1
            continue
        #################---REVERSE STRAND---#################

        # Increment index for next position in the base characters string
        index+=1

    # Create base dictionary to make it easier to get the most prevalent nucleotide
    base_counts = {'A':A_count_for+A_count_rev, 'C':C_count_for+C_count_rev, 'G':G_count_for+G_count_rev, 'T':T_count_for+T_count_rev}
    
    most_prevalent_base = max(base_counts, key=base_counts.get)

    # Total the base counts
    total = (A_count_for+A_count_rev 
            + C_count_for+C_count_rev 
            + G_count_for+G_count_rev 
            + T_count_for+T_count_rev
            + insertions+deletions)
    
    # Create new line and write it to the file
    new_line = (f"{reference},{total},{position},{reference_base},{most_prevalent_base},"
                f"{A_count_for},{A_count_rev},{C_count_for},{C_count_rev},{G_count_for},"
                f"{G_count_rev},{T_count_for},{T_count_rev},{insertions},{deletions},"
                f"{insert_10},{insert_50},{insert_100},{insert_500},{insert_1000},{insert_1500},{insert_GT_1500},"
                f"{dels_10},{dels_50},{dels_100},{dels_500},{dels_1000},{dels_1500},{dels_GT_1500}\n")  
    output_file.write(new_line)

# Close input and output files
mpileup_file.close()
output_file.close()