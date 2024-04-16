import os
# os.chdir('/Users/luli/Dropbox (UFL)/Tianqi/Lu_Tianqi/Exp 360. 2023-12-16 AQ-seq mouse huntington and human/Mouse/ISOFORM')
os.chdir("/Users/luli/Dropbox (UFL)/Tianqi/Lu_Tianqi/Exp 348.2023-09-19 BPS-AQ-seq miRNA analysis 293T S2/COLLAPSED/20231031/Isoform")
import pandas as pd

f1 = pd.read_csv("20230928_Collapsed_Isoform_miRNA_raw_count.csv") # Loading the DataFrame
f1['Isoform_length'] = f1['miRNA_sequence'].str.len() # Calculate the length of each string in the 'miRNA_sequence' column
f1.set_index(['miRNA_name', 'Isoform_length'], inplace=True) # Set 'miRNA_name' and 'Isoform_length' as the new index
f1.drop(columns=['miRNA_sequence'], inplace=True) # Delete the 'miRNA_sequence' column
compressed_df = f1.groupby(level=['miRNA_name', 'Isoform_length']).sum() #  Group by both indices and sum all other numeric columns
isoform_lengths = range(18, 31)  # 18 to 30 inclusive, Tianqi wants the isofrom length from whole 18nt to 30nt, each the value of each sample is 0.
miRNA_names = compressed_df.index.get_level_values('miRNA_name').unique() # Generate the new index: a product of unique miRNA_names and the desired range of Isoform_lengths
new_index = pd.MultiIndex.from_product([miRNA_names, isoform_lengths], names=['miRNA_name', 'Isoform_length'])
df_reindexed = compressed_df.reindex(new_index, fill_value=0) # Reindex the DataFrame with the new index, filling missing values with 0
df_reindexed.to_csv("20240103_miRNA_each_length_18ntTo30nt_exp348.csv")
