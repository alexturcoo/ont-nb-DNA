import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Example data (replace this with your actual data)
path = "/home/alextu/scratch/nonb_db_preprocessing/nonb_windows/all_non_b_motifs_0_based.csv"

# Paths to the second CSV file
path2 = "/home/alextu/scratch/chm13_nonb_db_preprocessing/nonb_windows/all_non_b_motifs_0_based.csv"

# Read in both CSV files
data1 = pd.read_csv(path)
data2 = pd.read_csv(path2)

# Remove rows where chr = chrMT
data1 = data1[data1['chr'] != 'chrM']
data2 = data2[data2['chr'] != 'chrM']

# Group data by 'chr' and 'feature', then sum the counts
grouped_data1 = data1.groupby(['chr', 'feature']).size().unstack(fill_value=0)
grouped_data2 = data2.groupby(['chr', 'feature']).size().unstack(fill_value=0)

# Compute the difference between the two datasets
heatmap_diff = grouped_data2 - grouped_data1

# Create the heatmap
plt.figure(figsize=(14, 14))
sns.heatmap(heatmap_diff, cmap='coolwarm', annot=True, fmt='d', linewidths=.5)

# Set plot labels and title
plt.xlabel('Feature Types')
plt.ylabel('Chromosomes')
plt.title('Differences in Structure Counts along Chromosomes')

# Show plot
plt.savefig('/home/alextu/scratch/nonb_db_preprocessing/imgs/heatmap_structures_along_chromosomes_hg38vschm13.png')

plt.show()
