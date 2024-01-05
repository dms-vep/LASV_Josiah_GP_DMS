
# Description:
# Parse the csv files with the listed variants 
# and classify the variants as either synonymous, 
# nonsynonymous, or in the UTR region.

# Author:
# Caleb Carr


# Imports
import pandas as pd 

# Dictionary of codon values and corresponding amino acids
codon_chart = {
    'TTT': 'Phe', 'TTC': 'Phe', 'TTA': 'Leu', 'TTG': 'Leu',
    'CTT': 'Leu', 'CTC': 'Leu', 'CTA': 'Leu', 'CTG': 'Leu',
    'ATT': 'Ile', 'ATC': 'Ile', 'ATA': 'Ile', 'ATG': 'Met',
    'GTT': 'Val', 'GTC': 'Val', 'GTA': 'Val', 'GTG': 'Val',
    'TCT': 'Ser', 'TCC': 'Ser', 'TCA': 'Ser', 'TCG': 'Ser',
    'CCT': 'Pro', 'CCC': 'Pro', 'CCA': 'Pro', 'CCG': 'Pro',
    'ACT': 'Thr', 'ACC': 'Thr', 'ACA': 'Thr', 'ACG': 'Thr',
    'GCT': 'Ala', 'GCC': 'Ala', 'GCA': 'Ala', 'GCG': 'Ala',
    'TAT': 'Tyr', 'TAC': 'Tyr', 'TAA': 'Stp', 'TAG': 'Stp',
    'CAT': 'His', 'CAC': 'His', 'CAA': 'Gln', 'CAG': 'Gln',
    'AAT': 'Asn', 'AAC': 'Asn', 'AAA': 'Lys', 'AAG': 'Lys',
    'GAT': 'Asp', 'GAC': 'Asp', 'GAA': 'Glu', 'GAG': 'Glu',
    'TGT': 'Cys', 'TGC': 'Cys', 'TGA': 'Stp', 'TGG': 'Trp',
    'CGT': 'Arg', 'CGC': 'Arg', 'CGA': 'Arg', 'CGG': 'Arg',
    'AGT': 'Ser', 'AGC': 'Ser', 'AGA': 'Arg', 'AGG': 'Arg',
    'GGT': 'Gly', 'GGC': 'Gly', 'GGA': 'Gly', 'GGG': 'Gly',
}


def get_codon(codon_position, position, alignment_data, gene_name):
    """ Given a variant codon position (0, 1, or 2), gene position,
        alignment summary datafram, and gene name, return the 
        major codon. """
    
    # Extract data for specific gene
    gene_data = alignment_data.loc[alignment_data['Reference'] == gene_name]

    # Initialize return value for codon
    codon = ""
    if codon_position == 0:
        # Position is the genomic position (starts at 1) 
        # while the index is 0 indexed. Therefore, position 
        # needs to be adjusted accordingly to access the 
        # create bases.
        codon += gene_data['Major Base'].iloc[position-1]
        codon += gene_data['Major Base'].iloc[position]
        codon += gene_data['Major Base'].iloc[position+1]
        
    elif codon_position == 1:
        # Position is the genomic position (starts at 1) 
        # while the index is 0 indexed. Therefore, position 
        # needs to be adjusted accordingly to access the 
        # create bases.
        codon += gene_data['Major Base'].iloc[position-2]
        codon += gene_data['Major Base'].iloc[position-1]
        codon += gene_data['Major Base'].iloc[position]

    else:
        # Position is the genomic position (starts at 1) 
        # while the index is 0 indexed. Therefore, position 
        # needs to be adjusted accordingly to access the 
        # create bases.
        codon += gene_data['Major Base'].iloc[position-3]
        codon += gene_data['Major Base'].iloc[position-2]
        codon += gene_data['Major Base'].iloc[position-1]

    # return the extracted codon
    return codon


def get_altered_codon(consensus_codon, codon_position, new_base):
    """ Given a consensus codon, codon position (0, 1, or 2), and 
        a the variant base, return the adjusted codon with the 
        new base. """
    # Slice the consensus codon and place new base in desired position
    return consensus_codon[:codon_position] + new_base + consensus_codon[codon_position+1:]


def mutation_type(consensus_codon, altered_codon):
    """ Given the consensus codon and the mutated codon, return 
        the corresponding mutation type depending on if the 
        coded amino acid changed. """
    # Return synonymous mutation if codons code for same amino acid
    if codon_chart[consensus_codon] == codon_chart[altered_codon]:
        return 'syn'
    # Return nonsense mutation if new codon codes for stop
    elif codon_chart[altered_codon] == 'Stp':
        return 'nonsense'
    # Return nonsynonymous mutation if codons code for different 
    # amino acids
    else:
        return 'nonsyn'


# Read input csv file that has variant data, all alignment data
# and coding region indices
variant_data = pd.read_csv(str(snakemake.input.variant_data))
alignment_data = pd.read_csv(str(snakemake.input.alignment_data))
coding_regions_indices = pd.read_csv(str(snakemake.input.coding_regions_indices))

# Open/create the output file accessed through the snakemake 
# object and allow it to be written to
output_file = open(str(snakemake.output), "w")

# Create header for output file
output_file.write(("Reference,Position,Synonymous,Syn Counts,Nonsynonymous,Nonsyn Counts,Nonsense,Nonsense Counts,UTR,"
                    "Insert Freq,Del Freq,"
                    "InsertsLTE10,InsertsLTE50,InsertsLTE100,InsertsLTE500,InsertsLTE1000,InsertsLTE1500,InsertsGT1500,"
                    "DelsLTE10,DelsLTE50,DelsLTE100,DelsLTE500,DelsLTE1000,DelsLTE1500,DelsGT1500\n"))

# Genomic positions of all variants
positions = variant_data['Position']

# Iterate through all positions
for i in range(len(positions)):

    # Extract gene name for current variant
    gene_name = variant_data['Reference'].iloc[i]

    # Extract data associated with gene
    gene_indices = coding_regions_indices.loc[coding_regions_indices['Reference'] == gene_name]

    # Extract starting and ending indices for protein coding regions
    start_CDS = gene_indices['Start CDS'].iloc[0]
    end_CDS = gene_indices['End CDS'].iloc[0]

    # Initialize values for writing new line
    # Gene name and position
    reference = gene_name
    position = positions.iloc[i]
    # Mutation frequencies and counts
    synonymous = 0
    syn_counts = 0
    nonsynonymous = 0
    nonsyn_counts = 0
    nonsense = 0
    nonsense_counts = 0
    utr = 0
    # Total insert frequencies and bin sized frequencies
    insert_freq = variant_data['Insert Freq'].iloc[i]
    insert_10 = variant_data['InsertsLTE10'].iloc[i]
    insert_50 = variant_data['InsertsLTE50'].iloc[i]
    insert_100 = variant_data['InsertsLTE100'].iloc[i]
    insert_500 = variant_data['InsertsLTE500'].iloc[i]
    insert_1000 = variant_data['InsertsLTE1000'].iloc[i]
    insert_1500 = variant_data['InsertsLTE1500'].iloc[i]
    insert_GT_1500 = variant_data['InsertsGT1500'].iloc[i]
    # Total del frquencies and bin sized frequencies
    del_freq = variant_data['Del Freq'].iloc[i]
    dels_10 = variant_data['DelsLTE10'].iloc[i]
    dels_50 = variant_data['DelsLTE50'].iloc[i]
    dels_100 = variant_data['DelsLTE100'].iloc[i]
    dels_500 = variant_data['DelsLTE500'].iloc[i]
    dels_1000 = variant_data['DelsLTE1000'].iloc[i]
    dels_1500 = variant_data['DelsLTE1500'].iloc[i]
    dels_GT_1500 = variant_data['DelsGT1500'].iloc[i]

    # If variant occurs outside of protein coding region, label it as a UTR variant
    if position < start_CDS or position > end_CDS:

        # Add all minor base frequencies for UTR variants because they cannot be classified
        # as synonymous or nonsynonymous
        utr = (variant_data['A Freq'].iloc[i] + variant_data['C Freq'].iloc[i] 
            + variant_data['G Freq'].iloc[i] + variant_data['T Freq'].iloc[i])
    
    # Variant occurs in coding region
    else:

        # Calculate the codon position of variant (0, 1, or 2)
        codon_position = (position - start_CDS) % 3
        # Get consensus codon for that position
        consensus_codon = get_codon(codon_position, position, alignment_data, gene_name)

        # If the minor base freq is greater than 0, then that base 
        # is considered a minor variant
        if variant_data['A Freq'].iloc[i] > 0:

            # Get altered codon for given base
            minor_codon = get_altered_codon(consensus_codon, codon_position, "A")

            # Check mutation type and add variant frequency to 
            # corresponding mutation type
            if mutation_type(consensus_codon, minor_codon) == 'syn':
                synonymous += variant_data['A Freq'].iloc[i]
                syn_counts += 1
            elif mutation_type(consensus_codon, minor_codon) == 'nonsyn':
                nonsynonymous += variant_data['A Freq'].iloc[i]
                nonsyn_counts += 1
            else:
                nonsense += variant_data['A Freq'].iloc[i]
                nonsense_counts += 1

        # If the minor base freq is greater than 0, then that base 
        # is considered a minor variant
        if variant_data['C Freq'].iloc[i] > 0:

            # Get altered codon for given base
            minor_codon = get_altered_codon(consensus_codon, codon_position, "C")

            # Check mutation type and add variant frequency to 
            # corresponding mutation type
            if mutation_type(consensus_codon, minor_codon) == 'syn':
                synonymous += variant_data['C Freq'].iloc[i]
                syn_counts += 1
            elif mutation_type(consensus_codon, minor_codon) == 'nonsyn':
                nonsynonymous += variant_data['C Freq'].iloc[i]
                nonsyn_counts += 1
            else:
                nonsense += variant_data['C Freq'].iloc[i]
                nonsense_counts += 1

        # If the minor base freq is greater than 0, then that base 
        # is considered a minor variant
        if variant_data['G Freq'].iloc[i] > 0:

            # Get altered codon for given base
            minor_codon = get_altered_codon(consensus_codon, codon_position, "G")

            # Check mutation type and add variant frequency to 
            # corresponding mutation type
            if mutation_type(consensus_codon, minor_codon) == 'syn':
                synonymous += variant_data['G Freq'].iloc[i]
                syn_counts += 1
            elif mutation_type(consensus_codon, minor_codon) == 'nonsyn':
                nonsynonymous += variant_data['G Freq'].iloc[i]
                nonsyn_counts += 1
            else:
                nonsense += variant_data['G Freq'].iloc[i]
                nonsense_counts += 1

        # If the minor base freq is greater than 0, then that base 
        # is considered a minor variant
        if variant_data['T Freq'].iloc[i] > 0:

            # Get altered codon for given base
            minor_codon = get_altered_codon(consensus_codon, codon_position, "T")

            # Check mutation type and add variant frequency to 
            # corresponding mutation type
            if mutation_type(consensus_codon, minor_codon) == 'syn':
                synonymous += variant_data['T Freq'].iloc[i]
                syn_counts += 1
            elif mutation_type(consensus_codon, minor_codon) == 'nonsyn':
                nonsynonymous += variant_data['T Freq'].iloc[i]
                nonsyn_counts += 1
            else:
                nonsense += variant_data['T Freq'].iloc[i]
                nonsense_counts += 1


    # Write mutation types to output file
    new_line = (f"{reference},{position},{synonymous},{syn_counts},{nonsynonymous},{nonsyn_counts},{nonsense},{nonsense_counts},{utr},"
                f"{insert_freq},{del_freq},"
                f"{insert_10},{insert_50},{insert_100},{insert_500},{insert_1000},{insert_1500},{insert_GT_1500},"
                f"{dels_10},{dels_50},{dels_100},{dels_500},{dels_1000},{dels_1500},{dels_GT_1500}\n") 
    output_file.write(new_line)

# Close output file
output_file.close()





