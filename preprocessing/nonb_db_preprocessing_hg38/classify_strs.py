from distutils.spawn import find_executable
from tkinter.messagebox import RETRY
from urllib.parse import non_hierarchical
import pandas as pd

# This script classifies STR Motifs from the non-B Database as mono, di, tri, etc nucleotide repeats...

def find_repeating_unit(sequence):
    """
    Finds the repeating unit in the sequence.
    """
    repeating_unit = ''
    for i in range(1, len(sequence)):
        if sequence[:i] * (len(sequence) // i) + sequence[:len(sequence) % i] == sequence:
            repeating_unit = sequence[:i]
            break
    return repeating_unit

def classify_str(sequence):
    """
    Classifies the STR sequence as mono, di, or trinucleotide repeat.
    """
    repeating_unit = find_repeating_unit(sequence)
    if len(repeating_unit) == 1:
        return "one"
    elif len(repeating_unit) == 2:
        return "two"
    elif len(repeating_unit) == 3:
        return "three"
    elif len(repeating_unit) == 4:
        return "four"
    elif len(repeating_unit) == 5:
        return "five"
    elif len(repeating_unit) == 6:
        return "six"
    elif len(repeating_unit) == 7:
        return "seven"
    elif len(repeating_unit) == 8:
        return "eight"
    elif len(repeating_unit) == 9:
        return "nine"
    else:
        return "Other"

# Example usage:
#motif_1 = "AAAAAAAAAAAAAA"
#motif_2 = "ATAATATAATATAATATAATATAATATAATATAAT"
#motif_3 = "ATCATCATC"

#print(classify_str(motif_1))  # Output: ('Unknown', 0)
#print(classify_str(motif_2))  # Output: ('Dinucleotide', 3)
#print(classify_str(motif_3))  # Output: ('Trinucleotide', 3)

# Example data (replace this with your actual data)
path = "/home/alextu/scratch/nonb_db_preprocessing/nonb_windows/all_non_b_windows_0_based.csv"

#Read in Data
non_b_annotations = pd.read_csv(path)

# Remove rows where chr = chrMT
short_tandem_repeats_nonb = non_b_annotations[non_b_annotations['feature'] == 'Short_Tandem_Repeat']
print(short_tandem_repeats_nonb.shape) # get the dimensions of the df

# Apply the classification function to each motif in the 'motifs' column
short_tandem_repeats_nonb['Classification'] = short_tandem_repeats_nonb['motif_seq'].apply(classify_str)

# Output the classified dataframe
print(short_tandem_repeats_nonb.head(10))

# Count the unique occurrences in the 'classification' column
classification_counts = short_tandem_repeats_nonb['Classification'].value_counts()

# Display the counts
print(classification_counts)

short_tandem_repeats_nonb.to_csv('/home/alextu/scratch/nonb_db_preprocessing/nonb_windows/str_motifs_types.csv', index=False)