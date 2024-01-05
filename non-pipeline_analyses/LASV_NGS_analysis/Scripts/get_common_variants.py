
# Description:
# This script extracts all variants that
# are common among more than one sample.

# Author:
# Caleb Carr

# Imports
import pandas as pd 

# Read input csv file that has alignment and variant data
samples_list_variants = snakemake.input.variants
samples_list_alignemnt = snakemake.input.alignemnt

# Params options for optional gene naming provided
# through config file
provide_gene_names = str(snakemake.params.gene_name_setting)
optional_gene_names = snakemake.params.gene_names
expected_gene_list = snakemake.params.gene_list

# All sample nonsynonymous variants
positions = [str(x) for x in list(range(1,2000))]
columns = ['Sample', 'Gene'] + positions
all_nonsyn_variants = pd.DataFrame(columns=columns)


# Create a list of dataframes for each sample from 
# the variant file
sample_variant_data_frame_list = []
for sample in samples_list_variants:
    sample_variant_data_frame_list.append(pd.read_csv(str(sample)))

# Create a list of dataframes for each sample, get
# the number of genes for each sample (should be the
# same otherwise there will be an error), and a 
# dictionary for each sample and its corresponding 
# genes
sample_data_frame_list = []
num_genes_per_sample = []
sample_gene_list_dictionary = {}
for sample in samples_list_alignemnt:

    # Create a dataframe for sample and add it to the list
    sample_data_frame_list.append(pd.read_csv(str(sample)))

    # Get number of gene per sample
    curr_data_frame = pd.read_csv(str(sample))
    gene_list = curr_data_frame['Reference'].unique()
    num_genes_per_sample.append(len(gene_list))
    
    # Create dictionary of samples with corresponding 
    # lists of genes
    sample_gene_list_dictionary[str(sample)] = gene_list

# Make a list of generic names for each gene (gene_1, etc.)
generic_gene_names = []
for i in range(len(expected_gene_list)):
    generic_gene_names.append(f"Gene_{i+1}")

# Create dictionary for accession gene names and
# provided optional gene names if setting is 
# configured in the config file
optional_name_dictionary = {}
if provide_gene_names == 'yes_provide':
    # Create dictionary from the two sets of lists
    optional_name_dictionary = dict(zip(generic_gene_names, optional_gene_names)) 

# Dictionary for common variants
# The format is:
# Gene,Position : [syn, nonsyn, nonsense, utr, inserts, dels]
common_variant_dict = {}
# Iterate through list of genes 
for i, gene_name in enumerate(generic_gene_names):

    # Iterate through all samples
    for j, (sample_data_frame, variant_data_frame) in enumerate(zip(sample_data_frame_list, sample_variant_data_frame_list)):
        # Extract data for each corresponding gene if current sample has the expected number of genes
        if len(sample_gene_list_dictionary[samples_list_alignemnt[j]]) == len(expected_gene_list):
            gene_data = variant_data_frame.loc[variant_data_frame['Reference'] == sample_gene_list_dictionary[samples_list_alignemnt[j]][i]]


            # Create empty row for current sample nonsyn mutations
            new_nonsyn_row = [0] * 2001
            new_nonsyn_row[0] = str(samples_list_variants[j]).split('/')[-1].split('_')[0]
            if gene_name == 'Gene_1':
                new_nonsyn_row[1] = 'GPC'
            else:
                new_nonsyn_row[1] = 'NP'

            

            # Iterate through current variant dataframe for 
            # specified gene data
            for k in range(0, len(gene_data['Position'])):
                curr_row = gene_data.iloc[k]
                # Create key for current variant in the form of:
                # gene,position
                curr_key = f"{curr_row['Reference']},{curr_row['Position']}"
                # Create blank list for counts of each variant
                curr_item = [0, 0, 0, 0, 0, 0]

                # Update variants based on count found
                if curr_row['Synonymous'] > 0:
                    curr_item[0] += 1
                if curr_row['Nonsynonymous'] > 0:
                    curr_item[1] += 1
                    # Add new nonsyn data
                    nonsyn = curr_row['Nonsynonymous']
                    position = int(curr_row['Position']) + 1 # Adjust for offset
                    new_nonsyn_row[position] = nonsyn
                if curr_row['Nonsense'] > 0:
                    curr_item[2] += 1
                if curr_row['UTR'] > 0:
                    curr_item[3] += 1
                if curr_row['Insert Freq'] > 0:
                    curr_item[4] += 1
                if curr_row['Del Freq'] > 0:
                    curr_item[5] += 1

                # Check if variant is in dictionary
                if curr_key in common_variant_dict:
                    # Update exisiting item
                    prev_item = common_variant_dict[curr_key]
                    new_item = [x + y for x, y in zip(prev_item, curr_item)]
                    common_variant_dict[curr_key] = new_item
                else:
                    # Add new item
                    common_variant_dict[curr_key] = curr_item

            # Add current sample to nonsyn dataframe
            all_nonsyn_variants.loc[len(all_nonsyn_variants.index)] = new_nonsyn_row

# Open/create the output file accessed through the snakemake 
# object and allow it to be written to
output_file = open(str(snakemake.output), "w")

# Create header for output file
output_file.write(("Gene,Position,Synonymous,Nonsynonymous,Nonsense,UTR,Insertions,Deletions\n"))
        
# Iterate through dictionary of common variants
for key, item in common_variant_dict.items():

    # Split key to extract gene name and position
    sectioned_key = key.split(',')
    # Create and write new line to csv output file
    new_line = (f"{sectioned_key[0]},{sectioned_key[1]},{item[0]},{item[1]},{item[2]},{item[3]},{item[4]},{item[5]}\n") 
    output_file.write(new_line)

# Close output files
output_file.close()

# Write nonsyn df to csv
all_nonsyn_variants.to_csv(str(snakemake.params.all_nonsyn_variants), index=False)

    








