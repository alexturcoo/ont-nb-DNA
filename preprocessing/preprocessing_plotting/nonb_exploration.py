import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Example data (replace this with your actual data)
path = "/home/alextu/scratch/nonb_db_preprocessing/nonb_windows/all_non_b_motifs_0_based.csv"

#Read in Data
all_chromosomes_data = pd.read_csv(path)

# Remove rows where chr = chrMT
all_chromosomes_data = all_chromosomes_data[all_chromosomes_data['chr'] != 'chrM']
print(all_chromosomes_data.shape) # get the dimensions of the df
#print(all_chromosomes_data[1:10])

# Find the number of each type of feature
feature_counts = all_chromosomes_data['feature'].value_counts()
print("Number of Each type of Motif")
print(feature_counts)

#Find the longest and shortest read
# Find the shortest and longest motif length for each feature
min_max_lengths = {}
for feature in feature_counts.index:
    feature_data = all_chromosomes_data[all_chromosomes_data['feature'] == feature]
    min_length = feature_data['motif_len'].min()
    max_length = feature_data['motif_len'].max()
    min_max_lengths[feature] = (min_length, max_length)

# Print the feature counts along with the shortest and longest motif lengths
for feature, count in feature_counts.items():
    print(f"Feature: {feature}, Count: {count}, Shortest & Longest Lengths: {min_max_lengths[feature]}")

# Finding number of reads over a certain length
# Define the threshold
threshold = 120

# Count the number of motif lengths over the threshold
count_over_threshold = (all_chromosomes_data['motif_len'] > threshold).sum()

print("Threshold=", threshold)
print("Number of motif lengths over the threshold:", count_over_threshold)

# Count the number of structures per chromosome
structures_per_chromosome = all_chromosomes_data.groupby('chr')['feature'].value_counts()
g_quad_per_chromosome = all_chromosomes_data[all_chromosomes_data['feature'] == 'G_Quadruplex_Motif'].groupby('chr')['feature'].value_counts()

print(g_quad_per_chromosome)

# Calculate the percentage of each feature per chromosome
percentage_per_chromosome = all_chromosomes_data.groupby(['chr', 'feature']).size() /all_chromosomes_data.groupby('chr').size()
percentage_per_chromosome = percentage_per_chromosome.unstack()

# Plot the stacked bar plot for structures per chromosome
plt.figure(figsize=(10, 6))
percentage_per_chromosome.plot(kind='bar', stacked=True)
plt.title('Percentage of Structures per Chromosome')
plt.xlabel('Chromosome')
plt.ylabel('Percentage')
plt.xticks(rotation=90, ha='right')
plt.legend(title='Feature', bbox_to_anchor=(1.05, 1), loc='upper left')
plt.grid(axis='y', linestyle='--', alpha=0.7)
plt.tight_layout()
plt.savefig('/home/alextu/scratch/nonb_db_preprocessing/imgs/percentage_of_structures_per_chr_hg38.png')

# Plot the count of G-Quadruplex Structures across chromosomes
plt.figure(figsize=(10, 6))
g_quad_per_chromosome.plot(kind='bar', color='green')
plt.title('Count of G-Quadruplex Structures across chromosomes')
plt.xlabel('Chromosome')
plt.ylabel('Count')
plt.xticks(rotation=90, ha='right')
plt.legend(title='Feature', bbox_to_anchor=(1.05, 1), loc='upper left')
plt.grid(axis='y', linestyle='--', alpha=0.7)
plt.tight_layout()
plt.savefig('/home/alextu/scratch/nonb_db_preprocessing/imgs/num_of_gquad_structures_per_chr_hg38.png')

# Heatmap showing structures in regions
# Create a matrix to represent the presence/absence of features in chromosome regions
# Here, we're creating a matrix with rows as chromosomes and columns as feature types
heatmap_data = all_chromosomes_data.pivot_table(index='chr', columns='feature', aggfunc='size', fill_value=0)

# Create the heatmap
plt.figure(figsize=(14, 12))
sns.heatmap(heatmap_data, cmap='Blues', annot=True, fmt='d', linewidths=.5)

# Set plot labels and title
plt.xlabel('Feature Types')
plt.ylabel('Chromosomes')
plt.title('Heatmap of Structure Types along Chromosomes')

# Show plot
plt.savefig('/home/alextu/scratch/nonb_db_preprocessing/imgs/heatmap_structures_along_chromosomes_hg38.png')

plt.show()