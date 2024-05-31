import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os

# Create an empty list to store DataFrames
dfs = []

# Loop through all CSV files
directory = '/home/alextu/scratch/compute_windows/translocations_dfs/HG00268/'
for filename in os.listdir(directory):
    if filename.endswith(".csv"):
        # Load the DataFrame from the CSV file
        path = os.path.join(directory, filename)
        print(path)
        df = pd.read_csv(path)
        # Append the DataFrame to the list
        dfs.append(df)

# Concatenate all DataFrames into one
df = pd.concat(dfs, ignore_index=True)
# Filter out rows where feature is "G_Quadruplex_Motif" and strand is "-"
#df = df[~((df['feature'] == "G_Quadruplex_Motif") & (df['strand'] == "-"))]

# Select rows where feature is "G_Quadruplex_Motif" and strand is "-"
rows_to_flip = df[(df['feature'] == "G_Quadruplex_Motif") & (df['strand'] == "-") & (df['direction'] == "opposite")]

# Flip the order of columns for these rows
for index, row in rows_to_flip.iterrows():
    cols_to_flip = ['translocation_time_' + str(i) for i in range(1, 100)]
    df.loc[index, cols_to_flip] = df.loc[index, cols_to_flip].values[::-1]

# Count features for positive strand
positive_features = df[df['direction'] == 'same']['feature'].value_counts()

# Count features for negative strand
negative_features = df[df['direction'] == 'opposite']['feature'].value_counts()

print("Positive Strand Features:")
print(positive_features)

print("\nNegative Strand Features:")
print(negative_features)

#Define the elements and their names
elements = ['A_Phased_Repeat', 'G_Quadruplex_Motif', 'Inverted_Repeat', 'Mirror_Repeat', 'Direct_Repeat',
            'Short_Tandem_Repeat', 'Z_DNA_Motif']
elements_name = ['A Phased Repeat', 'G Quadruplex', 'Inverted Repeat', 'Mirror Repeat', 'Direct Repeat',
                 'Short Tandem Repeat', 'Z DNA']

#elements = ["G_Quadruplex_Motif"]
#elements_name = ["G Quadruplex"]
quantiles = [0.05, 0.25, 0.5, 0.75, 0.95]
colors = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple']
x = range(1, 100)

# Define the path to save plots
plot_path = '/home/alextu/scratch/compute_windows/imgs/trans_times_quantiles_plots/HG00268'
if not os.path.exists(plot_path):
    os.makedirs(plot_path)

cols = ['translocation_time_' + str(i) for i in range(1, 100)]

direction = ['same', 'opposite']
#df_control = df[df['feature'] == 'Control'].reset_index(drop=True)
df_control_dict = {dir: df[(df['feature'] == 'Control') & (df['direction'] == dir)].reset_index(drop=True) for dir in direction}

for elem_id in range(len(elements)):
    elem_name = elements_name[elem_id]
    elem = elements[elem_id]

    fig, axs = plt.subplots(1, 2, figsize=(20, 7))  # Create subplots with 2 columns

    for strand_id, strand in enumerate(['same', 'opposite']):  # Loop through strands
        df_control = df_control_dict[strand]
        df_elem_strand = df[(df['feature'] == elem) & (df['direction'] == strand)].reset_index(drop=True)

        signal = df_elem_strand.loc[:, cols].to_numpy()
        control_signal = df_control.loc[:, cols].to_numpy()

        ax = axs[strand_id]  # Select the appropriate axis based on strand_id

        if strand == 'opposite' and elem != "G_Quadruplex_Motif":
            signal = np.flip(signal, axis=1)

        for quant_id in range(len(quantiles)):
            quant = quantiles[quant_id]
            color = colors[quant_id]
            signal_quantile = np.quantile(signal, quant, axis=0)
            control_signal_quantile = np.quantile(control_signal, quant, axis=0)

            ax.plot(x, signal_quantile, color=color, label=elem_name + ' - Quantile ' + str(quant))
            ax.plot(x, control_signal_quantile, color=color, linestyle='--', label='Control - Quantile ' + str(quant))

        ax.set_xlabel('Positions in windows', fontsize=15)
        ax.set_ylabel('Value', fontsize=15)
        ax.legend(loc='upper right', bbox_to_anchor=(1.0, 1.0))
        ax.set_title(elem_name + ' - Strand ' + strand + "_HG00268_chrY_chr22", fontsize=20)

    
    plt.tight_layout()
    plt.savefig(os.path.join(plot_path, elem + '_Control_HG00268_chrY_chr22.png'))
    plt.close()