
# Description:
# Parse the samtools stats files generated
# for each sample sam file

# Author:
# Caleb Carr

# List of samtools stats files
stat_file_list = snakemake.input
# List of sample accesion numbers to 
# us as names in the output file
sample_names = snakemake.params.samples

# Open/create the output file accessed through the snakemake 
# object and allow it to be written to
output_file = open(str(snakemake.output), "w")

# Create header for output file
output_file.write(("Sample,Mapping Rate\n"))

# Index to access sample names
sample_index = 0
# Iterate through list of input files
for stat_file in stat_file_list:

    # Open current stat file
    curr_stat_file = open(str(stat_file), "r")
    
    # Count of total reads for current sample
    total_reads = 0
    # Count of mapped reads for current sample
    mapped_reads = 0
    # Flag for curr sample that has reads mapped
    curr_sample_written = False
    # Iterate through each line of the stats file
    for line in curr_stat_file:

        # Split the line based on colon
        sectioned = line.split(':')

        # Update total reads when title 
        # specifying total reads is found
        if sectioned[0] == 'raw total sequences':
            # Update total read count
            total_reads += int(sectioned[1])
        
        # Update mapped reads when title
        # specifying mapped reads is found
        if sectioned[0] == 'reads mapped':
            # Update mapped read count
            mapped_reads += int(sectioned[1])
        
        # Write new line once total and 
        # mapped reads were found
        if total_reads > 0 and mapped_reads > 0:
            # Write new line with sample accession number and mapping rate
            output_file.write(f"{sample_names[sample_index]},{mapped_reads/total_reads}\n")
            curr_sample_written = True
            break
            
    # Write new line with sample accession number and 0 if no mapping rate was found
    if curr_sample_written == False:    
        output_file.write(f"{sample_names[sample_index]},{0}\n")

    # Close current stat file
    curr_stat_file.close()
    # Update sample name index
    sample_index += 1
    
# Close output files
output_file.close()