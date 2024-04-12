from random import sample
import pandas as pd

# Load your first DataFrame with short tandem repeats and classification
#classification_df = pd.read_csv("/home/alextu/scratch/nonb_db_preprocessing/nonb_windows/str_motifs_types.csv")

# Load your second DataFrame with start and end positions
translocation_df = pd.read_csv("/home/alextu/scratch/compute_windows/translocations_dfs/NA12878/tts_chr21_all_nonbstructures_positive.csv")


# Filter translocation_df for Short_Tandem_Repeat features
filtered_translocation_df = translocation_df[translocation_df['feature'] == 'Control']
print(filtered_translocation_df.shape)

# Sample 10,000 random windows
sampled_df = filtered_translocation_df.sample(n=10000, random_state=42)  # You can adjust the random_state if needed
print(sampled_df.shape)
# Merge filtered_translocation_df with classification_df
#merged_df = pd.merge(filtered_translocation_df, classification_df, 
#                     left_on=['mapped_start', 'mapped_end'],
#                     right_on=['win_start', 'win_end'],
#                     how='left')

#print(merged_df)
sampled_df.to_csv("/home/alextu/scratch/compute_windows/translocations_dfs/NA12878/tts_chr21_pos_controls10k.csv")