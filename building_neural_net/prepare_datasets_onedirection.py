import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os

# Create an empty list to store DataFrames
dfs = []

# Loop through all CSV files
directory = '/home/alextu/scratch/compute_windows/translocations_dfs/HG00268/done'
for filename in os.listdir(directory):
    if filename.startswith("tts_chr"):
        # Load the DataFrame from the CSV file
        path = os.path.join(directory, filename)
        df = pd.read_csv(path)
        # Append the DataFrame to the list
        dfs.append(df)

# Concatenate all DataFrames into one
df = pd.concat(dfs, ignore_index=True)
# Filter out rows where feature is "G_Quadruplex_Motif" and strand is "-"
#df = df[~((df['feature'] == "G_Quadruplex_Motif") & (df['strand'] == "-"))]

#### FLIP ORIENTATION OF NEGATIVE G QUADS
# Select rows where feature is "G_Quadruplex_Motif" and strand is "-"
rows_to_flip = df[(df['feature'] == "G_Quadruplex_Motif") & (df['strand'] == "-") & (df['direction'] == "opposite")]

# Flip the order of columns for these rows
for index, row in rows_to_flip.iterrows():
    cols_to_flip = ['translocation_time_' + str(i) for i in range(1, 100)]
    df.loc[index, cols_to_flip] = df.loc[index, cols_to_flip].values[::-1]

df_same = df[df['direction'] == "same"]
df_opposite = df[df['direction'] == 'opposite']

# Count features for positive strand
positive_features = df[df['direction'] == 'same']['feature'].value_counts()

# Count features for negative strand
negative_features = df[df['direction'] == 'opposite']['feature'].value_counts()

print("Positive Strand Features:")
print(positive_features)

print("\nNegative Strand Features:")
print(negative_features)

### FILTER BASED ON COVERAGE TO ONLY PULL TTS WITH COVERAGE > 5
# Extract counts of unique ranges
range_counts_same = df_same.groupby(['mapped_start', 'mapped_end']).size().reset_index(name='count')
# Filter out combinations with counts less than 5 - already filtered from extract_tts
#filtered_ranges_same = range_counts_same[range_counts_same['count'] >= 3]
print(range_counts_same)

range_counts_opposite = df_opposite.groupby(['mapped_start', 'mapped_end']).size().reset_index(name='count')
#filtered_ranges_opposite = range_counts_opposite[range_counts_opposite['count'] >= 3]
print(range_counts_opposite)


# Merge the two DataFrames on mapped_start and mapped_end
common_ranges = pd.merge(range_counts_same, range_counts_opposite, on=['mapped_start', 'mapped_end'], suffixes=('_same', '_opposite'))

# Print common ranges
print("Common Ranges:")
print(common_ranges)

# Initialize a new DataFrame to store the results
result_df = pd.DataFrame()
#result_df_same = pd.DataFrame()
#result_df_opposite = pd.DataFrame()

# Define a function to calculate median translocation times for a range and strand
def calculate_median_translocation_times(range_df, flip=False):
    translocation_times = range_df.filter(like='translocation_time_').values
    #if flip:
    #    translocation_times = np.flip(translocation_times, axis=1)
    return np.median(translocation_times, axis=0)

# Iterate over each common range
for index, row in range_counts_same.iterrows():
    # Get the mapped start and end values for the common range
    mapped_start = row['mapped_start']
    mapped_end = row['mapped_end']
    
    # Filter the original DataFrames to get the translocation times for the common range
    same_range_df = df_same[(df_same['mapped_start'] == mapped_start) & (df_same['mapped_end'] == mapped_end)]
    #opposite_range_df = df_opposite[(df_opposite['mapped_start'] == mapped_start) & (df_opposite['mapped_end'] == mapped_end)]
    
    # Get the feature information for the common range
    feature = same_range_df['feature'].iloc[0]  # Assuming feature is the same for all rows in the range
    coverage = row['count']

    # Calculate median translocation times for both strands
    same_medians = calculate_median_translocation_times(same_range_df)
    #opposite_medians = calculate_median_translocation_times(opposite_range_df)
    
    # Concatenate the median translocation times for both strands
    #translocation_times = np.concatenate(same_medians)
    
    # Create dictionaries for forward and reverse median translocation times
    forward_dict = {f'forward_tt_{i}': same_medians[i-1] for i in range(1, 101)}
    #reverse_dict = {f'reverse_tt_{i}': opposite_medians[i-1] for i in range(1, 101)}
    
    # Combine both dictionaries into a single dictionary
    data = {**forward_dict}
    #data_forward = forward_dict
    #data_reverse = reverse_dict
    
    # Add the feature information to the dictionary
    data['feature'] = feature
    data['coverage'] = coverage
    #data_forward['feature'] = feature
    #data_reverse['feature'] = feature
    
    # Create a DataFrame from the dictionary
    result_df = result_df._append(data, ignore_index = True)
    #result_df_same = result_df_same._append(data_forward, ignore_index = True)
    #result_df_opposite = result_df_opposite._append(data_reverse, ignore_index = True)


### BDNA TRAIN VAL TEST SPLIT
# Define the percentages
train_percent_bdna = 0.3
test_percent_bdna = 0.5
val_percent_bdna = 0.2

bdna_df = result_df[result_df['feature'] == "Control"]

# Split the data into training, testing, and validation sets
total_rows = len(bdna_df)
train_size = int(total_rows * train_percent_bdna)
test_size = int(total_rows * test_percent_bdna)
val_size = total_rows - train_size - test_size

train_df_bdna = bdna_df[:train_size]
test_df_bdna = bdna_df[train_size:train_size + test_size]
val_df_bdna = bdna_df[train_size + test_size:]

# Define the file paths for saving
train_file_path_bdna = f"/home/alextu/scratch/compute_windows/prepared_datasets_bdna_oneway/bdna_train.csv"
test_file_path_bdna = f"/home/alextu/scratch/compute_windows/prepared_datasets_bdna_oneway/bdna_test.csv"
val_file_path_bdna = f"/home/alextu/scratch/compute_windows/prepared_datasets_bdna_oneway/bdna_validation.csv"

# Save each set to a CSV file
train_df_bdna.to_csv(train_file_path_bdna, index=False)
test_df_bdna.to_csv(test_file_path_bdna, index=False)
val_df_bdna.to_csv(val_file_path_bdna, index=False)

#bdna_df.to_csv("/home/alextu/scratch/compute_windows/prepared_datasets/full_bdna_dataset.csv")

non_b_types = ['A_Phased_Repeat', 'G_Quadruplex_Motif', 'Inverted_Repeat', 'Mirror_Repeat', 'Direct_Repeat',
               'Short_Tandem_Repeat', 'Z_DNA_Motif']

# Define the percentages
train_percent_nonb = 0.1
test_percent_nonb = 0.1
val_percent_nonb = 0.8

# Iterate over each non-B type feature
for feature in non_b_types:
    # Filter forward and reverse dataframes based on the feature
    filtered_df = result_df[result_df['feature'] == feature]

    # Split the data into training, testing, and validation sets
    total_rows = len(filtered_df)
    train_size = int(total_rows * train_percent_nonb)
    test_size = int(total_rows * test_percent_nonb)
    val_size = total_rows - train_size - test_size

    train_df_nonb = filtered_df[:train_size]
    test_df_nonb = filtered_df[train_size:train_size + test_size]
    val_df_nonb = filtered_df[train_size + test_size:]

    # Define the file paths for saving
    train_file_path_nonb = f"/home/alextu/scratch/compute_windows/prepared_datasets_nonb_dna_oneway/{feature}_train.csv"
    test_file_path_nonb = f"/home/alextu/scratch/compute_windows/prepared_datasets_nonb_dna_oneway/{feature}_test.csv"
    val_file_path_nonb = f"/home/alextu/scratch/compute_windows/prepared_datasets_nonb_dna_oneway/{feature}_validation.csv"

    # Save each set to a CSV file
    train_df_nonb.to_csv(train_file_path_nonb, index=False)
    test_df_nonb.to_csv(test_file_path_nonb, index=False)
    val_df_nonb.to_csv(val_file_path_nonb, index=False)

#Perform Robust scaling
#bdna_forward = result_df_same[result_df_same["feature"] == "Control"]
#bdna_forward = bdna_forward.drop(columns=['feature'])
#bdna_reverse = result_df_opposite[result_df_opposite["feature"] == "Control"]
#bdna_reverse = bdna_reverse.drop(columns=['feature'])

#iqr_fd = np.subtract(*np.percentile(bdna_forward, [75, 25], interpolation='midpoint'))
#iqr_rd = np.subtract(*np.percentile(bdna_reverse, [75, 25], interpolation='midpoint'))
#tfd_mu, tfd_std = np.median(bdna_forward), iqr_fd
#trd_mu, trd_std = np.median(bdna_reverse), iqr_rd

#print(iqr_fd, tfd_mu, iqr_rd, trd_mu)