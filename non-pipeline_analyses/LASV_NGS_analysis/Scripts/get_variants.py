
# Description:
# Based on the cutoffs listed in the config file,
# iterate through the minor base frequencies and
# extract the minor variants that meet the specified
# requirements.

# Author:
# Caleb Carr


# Imports
import pandas as pd 

# Read input csv file that has alignment data
alignment_data = pd.read_csv(str(snakemake.input))

# Open/create the output file accessed through the snakemake 
# object and allow it to be written to
output_file = open(str(snakemake.output), "w")

# Variant parameters that are specified in the config.yml file. 
# Specified as floats because many are decimal values
base_cutoff = float(snakemake.params.minor_base_cutoff)
insert_cutoff = float(snakemake.params.insert_cutoff)
del_cutoff = float(snakemake.params.del_cutoff)
min_bias = float(snakemake.params.min_ratio_bias)
max_bias = float(snakemake.params.max_ratio_bias)
min_depth = float(snakemake.params.min_depth)


# Create header for output file
output_file.write(("Reference,Position,A Freq,C Freq,G Freq,T Freq,"
                    "Insert Freq,Del Freq,"
                    "InsertsLTE10,InsertsLTE50,InsertsLTE100,InsertsLTE500,InsertsLTE1000,InsertsLTE1500,InsertsGT1500,"
                    "DelsLTE10,DelsLTE50,DelsLTE100,DelsLTE500,DelsLTE1000,DelsLTE1500,DelsGT1500\n"))

# Extract base positions
positions = alignment_data['Position']

# Iterate through all positions
for i in range(len(positions)):

    # Check if position meets the minimum depth requirement
    if alignment_data['Depth'].iloc[i] >= min_depth:

        # Boolean for writing a new line for a variant.
        # Only true when a variant is found and written 
        # to output file.
        write_new_line = False
        # Initializing values to be written to the output file
        # Reference and position
        reference = alignment_data['Reference'].iloc[i]
        position = alignment_data['Position'].iloc[i]
        # Base frequencies
        A_freq = 0
        C_freq = 0
        G_freq = 0
        T_freq = 0
        # Total indel frequencies
        insert_freq = 0
        del_freq = 0
        # Insert bin size frequencies
        insert_10 = 0
        insert_50 = 0
        insert_100 = 0
        insert_500 = 0
        insert_1000 = 0
        insert_1500 = 0
        insert_GT_1500 = 0
        # Dels bin size frequencies
        dels_10 = 0
        dels_50 = 0
        dels_100 = 0
        dels_500 = 0
        dels_1000 = 0
        dels_1500 = 0
        dels_GT_1500 = 0


        # If given base meets the cutoff and strand bias requirements, then write 
        # the base freq to the output file
        if (alignment_data['A Freq'].iloc[i] >= base_cutoff and 
            (alignment_data['A Bias'].iloc[i] > min_bias and alignment_data['A Bias'].iloc[i] < max_bias)):

            # Update boolean for writing a new line and base freq for variant
            write_new_line = True
            A_freq += alignment_data['A Freq'].iloc[i]

        # If given base meets the cutoff and strand bias requirements, then write 
        # the base freq to the output file      
        if (alignment_data['C Freq'].iloc[i] >= base_cutoff and 
            (alignment_data['C Bias'].iloc[i] > min_bias and alignment_data['C Bias'].iloc[i] < max_bias)):

            # Update boolean for writing a new line and base freq for variant
            write_new_line = True
            C_freq += alignment_data['C Freq'].iloc[i]

        # If given base meets the cutoff and strand bias requirements, then write 
        # the base freq to the output file
        if (alignment_data['G Freq'].iloc[i] >= base_cutoff and 
            (alignment_data['G Bias'].iloc[i] > min_bias and alignment_data['G Bias'].iloc[i] < max_bias)):

            # Update boolean for writing a new line and base freq for variant
            write_new_line = True
            G_freq += alignment_data['G Freq'].iloc[i]

        # If given base meets the cutoff and strand bias requirements, then write 
        # the base freq to the output file
        if (alignment_data['T Freq'].iloc[i] >= base_cutoff and 
            (alignment_data['T Bias'].iloc[i] > min_bias and alignment_data['T Bias'].iloc[i] < max_bias)):

            # Update boolean for writing a new line and base freq for variant
            write_new_line = True
            T_freq += alignment_data['T Freq'].iloc[i]

        # If insertion meets the cutoff requirement, then write the
        # total insert frequency and bin size insert frequencies to
        # output file
        if alignment_data['Insert Freq'].iloc[i] >= insert_cutoff:

            # Update boolean for writing a new line and insert freq for variant
            write_new_line = True
            insert_freq += alignment_data['Insert Freq'].iloc[i]
            # Update del bin size frequency if above cutoff
            if alignment_data['InsertsLTE10'].iloc[i] >= insert_cutoff:
                insert_10 += alignment_data['InsertsLTE10'].iloc[i]
            if alignment_data['InsertsLTE50'].iloc[i] >= insert_cutoff:
                insert_50 += alignment_data['InsertsLTE50'].iloc[i]
            if alignment_data['InsertsLTE100'].iloc[i] >= insert_cutoff:
                insert_100 += alignment_data['InsertsLTE100'].iloc[i]
            if alignment_data['InsertsLTE500'].iloc[i] >= insert_cutoff:
                insert_500 += alignment_data['InsertsLTE500'].iloc[i]
            if alignment_data['InsertsLTE1000'].iloc[i] >= insert_cutoff:
                insert_1000 += alignment_data['InsertsLTE1000'].iloc[i]
            if alignment_data['InsertsLTE1500'].iloc[i] >= insert_cutoff:
                insert_1500 += alignment_data['InsertsLTE1500'].iloc[i]
            if alignment_data['InsertsGT1500'].iloc[i] >= insert_cutoff:
                insert_GT_1500 += alignment_data['InsertsGT1500'].iloc[i]
            

        # If deletion meets the cutoff requirement, then write the 
        # total del frequencies and bin size del frequencies to 
        # output file
        if alignment_data['Del Freq'].iloc[i] >= del_cutoff:

            # Update boolean for writing a new line and del freq for variant
            write_new_line = True
            del_freq += alignment_data['Del Freq'].iloc[i]
            # Update del bin size frequency if above cutoff
            if alignment_data['DelsLTE10'].iloc[i] >= insert_cutoff:
                dels_10 += alignment_data['DelsLTE10'].iloc[i]
            if alignment_data['DelsLTE50'].iloc[i] >= insert_cutoff:
                dels_50 += alignment_data['DelsLTE50'].iloc[i]
            if alignment_data['DelsLTE100'].iloc[i] >= insert_cutoff:
                dels_100 += alignment_data['DelsLTE100'].iloc[i]
            if alignment_data['DelsLTE500'].iloc[i] >= insert_cutoff:
                dels_500 += alignment_data['DelsLTE500'].iloc[i]
            if alignment_data['DelsLTE1000'].iloc[i] >= insert_cutoff:
                dels_1000 += alignment_data['DelsLTE1000'].iloc[i]
            if alignment_data['DelsLTE1500'].iloc[i] >= insert_cutoff:
                dels_1500 += alignment_data['DelsLTE1500'].iloc[i]
            if alignment_data['DelsGT1500'].iloc[i] >= insert_cutoff:
                dels_GT_1500 += alignment_data['DelsGT1500'].iloc[i]
            

        # If boolean for writing a new line is true, then a variant 
        # was found and the information is written to the output file
        if write_new_line == True:

            # Write variant information to output file
            new_line = (f"{reference},{position},{A_freq},{C_freq},{G_freq},{T_freq},"
                        f"{insert_freq},{del_freq},"
                        f"{insert_10},{insert_50},{insert_100},{insert_500},{insert_1000},{insert_1500},{insert_GT_1500},"
                        f"{dels_10},{dels_50},{dels_100},{dels_500},{dels_1000},{dels_1500},{dels_GT_1500}\n") 
            output_file.write(new_line)
    

# Close output files
output_file.close()