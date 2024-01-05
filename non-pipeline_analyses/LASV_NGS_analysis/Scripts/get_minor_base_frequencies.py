
# Description:
# Parse the mpileup summary to get minor base frequencies,
# strand ratio bias, and indel frequencies.

# Author:
# Caleb Carr


# Imports
import pandas as pd 
import statistics
import math

# This function is used to calculate strand bias
def get_ratio(a, b):
    """ Calculate ratio of two numbers """
    if a == 0 and b == 0:
        return 0
    else:
        return a/(a+b)


# Read input csv file that has alignment data
alignment_data = pd.read_csv(str(snakemake.input))

# Open/create the output file accessed through the snakemake 
# object and allow it to be written to
output_file = open(str(snakemake.output), "w")

# Create header for output file
output_file.write(("Reference,Position,Depth,Major Base,Minor Bases Freq,"
                    "A Freq,C Freq,G Freq,T Freq,"
                    "A Bias,C Bias,G Bias,T Bias,"
                    "Insert Freq,Del Freq,"
                    "InsertsLTE10,InsertsLTE50,InsertsLTE100,InsertsLTE500,InsertsLTE1000,InsertsLTE1500,InsertsGT1500,"
                    "DelsLTE10,DelsLTE50,DelsLTE100,DelsLTE500,DelsLTE1000,DelsLTE1500,DelsGT1500\n"))

# Extract all aligned gene reference names
gene_list = alignment_data['Reference'].unique()
# Extract sample consensus sequence bases
consensus_sample_bases = alignment_data['Consensus Base']
# Extract base counts for forward and reverse strands
A_count_For = alignment_data['A For'] 
C_count_For = alignment_data['C For'] 
G_count_For = alignment_data['G For'] 
T_count_For = alignment_data['T For'] 
A_count_Rev = alignment_data['A Rev'] 
C_count_Rev = alignment_data['C Rev'] 
G_count_Rev = alignment_data['G Rev'] 
T_count_Rev = alignment_data['T Rev'] 
# Extract read coverage for all positions
depth = alignment_data['Depth']
# Extract indel counts for all positions
inserts = alignment_data['Insertions']
delets = alignment_data['Deletions']
# Extract insert bin size counts
insert_10 = alignment_data['InsertsLTE10']
insert_50 = alignment_data['InsertsLTE50']
insert_100 = alignment_data['InsertsLTE100']
insert_500 = alignment_data['InsertsLTE500']
insert_1000 = alignment_data['InsertsLTE1000']
insert_1500 = alignment_data['InsertsLTE1500']
insert_GT_1500 = alignment_data['InsertsGT1500']
# Extract dels bin size counts
dels_10 = alignment_data['DelsLTE10']
dels_50 = alignment_data['DelsLTE50']
dels_100 = alignment_data['DelsLTE100']
dels_500 = alignment_data['DelsLTE500']
dels_1000 = alignment_data['DelsLTE1000']
dels_1500 = alignment_data['DelsLTE1500']
dels_GT_1500 = alignment_data['DelsGT1500']

# Iterate through all the consensus bases
for i in range(len(consensus_sample_bases)):

    # If the consensus base is A, then increment all other  
    # base percentages as the minor base divided by total 
    # base counts.
    if consensus_sample_bases.iloc[i] == 'A':
        total = depth.iloc[i]
        C_percentage = ((C_count_For.iloc[i] + C_count_Rev.iloc[i])/total) * 100
        G_percentage = ((G_count_For.iloc[i] + G_count_Rev.iloc[i])/total) * 100
        T_percentage = ((T_count_For.iloc[i] + T_count_Rev.iloc[i])/total) * 100
        A_percentage = 0 # major base has a 0 percentage
        # total minor base frequencies
        total_percentage = C_percentage + G_percentage + T_percentage 
        # strand bias ratios for each base
        A_bias = get_ratio(A_count_For.iloc[i], A_count_Rev.iloc[i])
        C_bias = get_ratio(C_count_For.iloc[i], C_count_Rev.iloc[i])
        G_bias = get_ratio(G_count_For.iloc[i], G_count_Rev.iloc[i])
        T_bias = get_ratio(T_count_For.iloc[i], T_count_Rev.iloc[i])
        # Total indel frequencies. Bin size frequencies 
        # are calculated within the f string.
        insert_freq = (inserts.iloc[i]/total) * 100
        indel_freq = (delets.iloc[i]/total) * 100
        # Write new line to output file
        new_line = (f"{alignment_data['Reference'].iloc[i]},"
                    f"{alignment_data['Position'].iloc[i]},"
                    f"{total},"
                    f"{consensus_sample_bases.iloc[i]},"
                    f"{total_percentage},"
                    f"{A_percentage},{C_percentage},{G_percentage},{T_percentage},"
                    f"{A_bias},{C_bias},{G_bias},{T_bias},"
                    f"{insert_freq},{indel_freq},"
                    f"{(insert_10.iloc[i]/total)*100},{(insert_50.iloc[i]/total)*100},{(insert_100.iloc[i]/total)*100},"
                    f"{(insert_500.iloc[i]/total)*100},{(insert_1000.iloc[i]/total)*100},{(insert_1500.iloc[i]/total)*100},{(insert_GT_1500.iloc[i]/total)*100},"
                    f"{(dels_10.iloc[i]/total)*100},{(dels_50.iloc[i]/total)*100},{(dels_100.iloc[i]/total)*100},"
                    f"{(dels_500.iloc[i]/total)*100},{(dels_1000.iloc[i]/total)*100},{(dels_1500.iloc[i]/total)*100},{(dels_GT_1500.iloc[i]/total)*100}\n")
        output_file.write(new_line)
        continue

    # If the consensus base is C, then increment all other 
    # base percentages as the minor base divided by total 
    # base counts.
    if consensus_sample_bases.iloc[i] == 'C':
        total = depth.iloc[i]
        A_percentage = ((A_count_For.iloc[i] + A_count_Rev.iloc[i])/total) * 100
        G_percentage = ((G_count_For.iloc[i] + G_count_Rev.iloc[i])/total) * 100
        T_percentage = ((T_count_For.iloc[i] + T_count_Rev.iloc[i])/total) * 100
        C_percentage = 0 # major base has a 0 percentage
        # total minor base frequencies
        total_percentage = A_percentage + G_percentage + T_percentage
        # strand bias ratios for each base
        A_bias = get_ratio(A_count_For.iloc[i], A_count_Rev.iloc[i])
        C_bias = get_ratio(C_count_For.iloc[i], C_count_Rev.iloc[i])
        G_bias = get_ratio(G_count_For.iloc[i], G_count_Rev.iloc[i])
        T_bias = get_ratio(T_count_For.iloc[i], T_count_Rev.iloc[i])
        # Total indel frequencies. Bin size frequencies 
        # are calculated within the f string.
        insert_freq = (inserts.iloc[i]/total) * 100
        indel_freq = (delets.iloc[i]/total) * 100
        # Write new line to output file
        new_line = (f"{alignment_data['Reference'].iloc[i]},"
                    f"{alignment_data['Position'].iloc[i]},"
                    f"{total},"
                    f"{consensus_sample_bases.iloc[i]},"
                    f"{total_percentage},"
                    f"{A_percentage},{C_percentage},{G_percentage},{T_percentage},"
                    f"{A_bias},{C_bias},{G_bias},{T_bias},"
                    f"{insert_freq},{indel_freq},"
                    f"{(insert_10.iloc[i]/total)*100},{(insert_50.iloc[i]/total)*100},{(insert_100.iloc[i]/total)*100},"
                    f"{(insert_500.iloc[i]/total)*100},{(insert_1000.iloc[i]/total)*100},{(insert_1500.iloc[i]/total)*100},{(insert_GT_1500.iloc[i]/total)*100},"
                    f"{(dels_10.iloc[i]/total)*100},{(dels_50.iloc[i]/total)*100},{(dels_100.iloc[i]/total)*100},"
                    f"{(dels_500.iloc[i]/total)*100},{(dels_1000.iloc[i]/total)*100},{(dels_1500.iloc[i]/total)*100},{(dels_GT_1500.iloc[i]/total)*100}\n")
        output_file.write(new_line)
        continue

    # If the consensus base is G, then increment all other 
    # base percentages as the minor base divided by total 
    # base counts.
    if consensus_sample_bases.iloc[i] == 'G':
        total = depth.iloc[i]
        A_percentage = ((A_count_For.iloc[i] + A_count_Rev.iloc[i])/total) * 100
        C_percentage = ((C_count_For.iloc[i] + C_count_Rev.iloc[i])/total) * 100
        T_percentage = ((T_count_For.iloc[i] + T_count_Rev.iloc[i])/total) * 100
        G_percentage = 0 # major base has a 0 percentage
        # total minor base frequencies
        total_percentage = A_percentage + C_percentage + T_percentage
        # strand bias ratios for each base
        A_bias = get_ratio(A_count_For.iloc[i], A_count_Rev.iloc[i])
        C_bias = get_ratio(C_count_For.iloc[i], C_count_Rev.iloc[i])
        G_bias = get_ratio(G_count_For.iloc[i], G_count_Rev.iloc[i])
        T_bias = get_ratio(T_count_For.iloc[i], T_count_Rev.iloc[i])
        # Total indel frequencies. Bin size frequencies 
        # are calculated within the f string.
        insert_freq = (inserts.iloc[i]/total) * 100
        indel_freq = (delets.iloc[i]/total) * 100
        # Write new line to output file
        new_line = (f"{alignment_data['Reference'].iloc[i]},"
                    f"{alignment_data['Position'].iloc[i]},"
                    f"{total},"
                    f"{consensus_sample_bases.iloc[i]},"
                    f"{total_percentage},"
                    f"{A_percentage},{C_percentage},{G_percentage},{T_percentage},"
                    f"{A_bias},{C_bias},{G_bias},{T_bias},"
                    f"{insert_freq},{indel_freq},"
                    f"{(insert_10.iloc[i]/total)*100},{(insert_50.iloc[i]/total)*100},{(insert_100.iloc[i]/total)*100},"
                    f"{(insert_500.iloc[i]/total)*100},{(insert_1000.iloc[i]/total)*100},{(insert_1500.iloc[i]/total)*100},{(insert_GT_1500.iloc[i]/total)*100},"
                    f"{(dels_10.iloc[i]/total)*100},{(dels_50.iloc[i]/total)*100},{(dels_100.iloc[i]/total)*100},"
                    f"{(dels_500.iloc[i]/total)*100},{(dels_1000.iloc[i]/total)*100},{(dels_1500.iloc[i]/total)*100},{(dels_GT_1500.iloc[i]/total)*100}\n")
        output_file.write(new_line)
        continue

    # If the consensus base is T, then increment all other 
    # base percentages as the minor base divided by total 
    # base counts.
    if consensus_sample_bases.iloc[i] == 'T':
        total = depth.iloc[i]
        A_percentage = ((A_count_For.iloc[i] + A_count_Rev.iloc[i])/total) * 100
        C_percentage = ((C_count_For.iloc[i] + C_count_Rev.iloc[i])/total) * 100
        G_percentage = ((G_count_For.iloc[i] + G_count_Rev.iloc[i])/total) * 100
        T_percentage = 0 # major base has a 0 percentage
        # total minor base frequencies
        total_percentage = A_percentage + C_percentage + G_percentage
        # strand bias ratios for each base
        A_bias = get_ratio(A_count_For.iloc[i], A_count_Rev.iloc[i])
        C_bias = get_ratio(C_count_For.iloc[i], C_count_Rev.iloc[i])
        G_bias = get_ratio(G_count_For.iloc[i], G_count_Rev.iloc[i])
        T_bias = get_ratio(T_count_For.iloc[i], T_count_Rev.iloc[i])
        # Total indel frequencies. Bin size frequencies 
        # are calculated within the f string.
        insert_freq = (inserts.iloc[i]/total) * 100
        indel_freq = (delets.iloc[i]/total) * 100
        # Write new line to output file
        new_line = (f"{alignment_data['Reference'].iloc[i]},"
                    f"{alignment_data['Position'].iloc[i]},"
                    f"{total},"
                    f"{consensus_sample_bases.iloc[i]},"
                    f"{total_percentage},"
                    f"{A_percentage},{C_percentage},{G_percentage},{T_percentage},"
                    f"{A_bias},{C_bias},{G_bias},{T_bias},"
                    f"{insert_freq},{indel_freq},"
                    f"{(insert_10.iloc[i]/total)*100},{(insert_50.iloc[i]/total)*100},{(insert_100.iloc[i]/total)*100},"
                    f"{(insert_500.iloc[i]/total)*100},{(insert_1000.iloc[i]/total)*100},{(insert_1500.iloc[i]/total)*100},{(insert_GT_1500.iloc[i]/total)*100},"
                    f"{(dels_10.iloc[i]/total)*100},{(dels_50.iloc[i]/total)*100},{(dels_100.iloc[i]/total)*100},"
                    f"{(dels_500.iloc[i]/total)*100},{(dels_1000.iloc[i]/total)*100},{(dels_1500.iloc[i]/total)*100},{(dels_GT_1500.iloc[i]/total)*100}\n")
        output_file.write(new_line)
        continue

# Close output files
output_file.close()
