#Written by Frances shepherd, 2/18/2021
#Description of code:
#Imports: 
#    1. A csv file containing metadata for each sample that will be added to the dataframe
#    2. a string that defines the list of files path where the filtered variant tsv files are located (e.g. "./variant_calling/variants/"). The code will look for the iva output files that are the FILTERED list of variants (i.e. "M1_AK04-LIV-astro-rdrp-D2.filtered.tsv")
#    3. a string that defines the name of the output file of combined variants, with .csv extension (e.g. "combined_filtered_variants.csv")
#
#Outputs a csv file that combines the filtered variant call data along with metadata that further describes the sample from where the reads came, and an average of the variant frequency between the two sequencing replicates. 
#
#
#Usage: ./combine_variant_files.py <metadata_file_name.csv> <path_for_variants> <output_file_name.csv>
#----------------------------------------------------------------------------------------

########################################################################################
#There is no need to mess with the rest of the code below this line
#######################################################################################
import glob
import pandas as pd
import os
import string
import sys

#Read in text file defined in code initiation and assign to variables
metadata = sys.argv[1]
path = sys.argv[2]
output_name = sys.argv[3]


all_files = glob.glob(path+"/*filtered.tsv") #create list of the filtered variants
meta = pd.read_csv(metadata) #read in metadata to add rows to the variant data

data = [] #Create empty list that will be populated by the filtered variant dataframes. The concat command takes a list of dataframes

#Define column names
new_columns =["REGION", "POS", "REF", "ALT", "GFF_FEATURE", "REF_CODON", "REF_AA", "ALT_CODON", "ALT_AA", "REF_DP_repA",
          "REF_RV_repA", "REF_QUAL_repA", "ALT_DP_repA", "ALT_RV_repA", "ALT_QUAL_repA", "ALT_FREQ_repA",
          "TOTAL_DP_repA", "PVAL_repA", "PASS_repA", "REF_DP_repB", "REF_RV_repB", "REF_QUAL_repB", "ALT_DP_repB",
          "ALT_RV_repB", "ALT_QUAL_repB", "ALT_FREQ_repB", "TOTAL_DP_repB", "PVAL_repB", "PASS_repB", "filename"]

for f in all_files: #For each file in the list of tsv files
    frame = pd.read_csv(f,sep="\t",header=None,skiprows = 1) #read the tsv file, split by the tab character, ignore headers and skip the first row which corresponds to row number
    frame['filename'] = os.path.basename(f) #Add a column to each .tsv file that lists the file name so the metadata can refer to it
    data.append(frame) #Append the altered tsv file to the empty list
    

combined_csv = pd.concat(data) #Concat all the tsv dataframes into one
combined_csv.columns = new_columns #Rename column names
new_df = pd.merge(meta, combined_csv, how='inner', on='filename') #Merge the metadata and the filtered variant files from ivar together
new_df['average_variant_freq'] = new_df[['ALT_FREQ_repA', 'ALT_FREQ_repB']].mean(axis=1) #Calculate the average variant frequency from the 2 replicates
new_df.to_csv(output_name, index = False) #Save the file as a new csv file, taking the name specified in the input argument.