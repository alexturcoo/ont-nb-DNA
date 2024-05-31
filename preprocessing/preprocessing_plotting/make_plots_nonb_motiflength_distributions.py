import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import numpy as np

### PLOTTING MOTIF LENGTH OF NON-B DNA FROM non-B DNA DB

path = '/home/alextu/scratch/nonb_db_preprocessing/nonb_windows/all_non_b_motifs_0_based.csv'
chromosome_df = pd.read_csv(path)
chromosome_df = chromosome_df.sort_values(by='start').reset_index(drop=True)
chromosome_df = chromosome_df[chromosome_df['motif_len']<= 150].reset_index(drop=True)
chromosome_df = chromosome_df[chromosome_df['feature']== "Short_Tandem_Repeat"].reset_index(drop=True)
p = sns.displot(chromosome_df, x="motif_len", hue="feature", kind="kde", fill=True)
# p.set(xlabel='Length')
p.set(xlabel='Length', title='Potential Non-B DNA structures\n(Motifs w/ length > 150 are cut out)')
p._legend.set_title('Non-B DNA structures')
plt.tight_layout()
plt.savefig('/home/alextu/scratch/nonb_db_preprocessing/imgs/motif_length_density_l150_hg38_STR.png')

### PLOT BY CHROMOSOME
#p = sns.displot(chromosome_df, x="motif_len", hue="feature", kind="kde", col="chr", fill=True, col_wrap=5)
# p.set_xlabel('Length', fontsize=20)
# plt.show()
# p.set(xlabel='Length', title='Potential Non-B DNA structures\n(motifs with length > 150 (0.12% of motifs)
# are cut out )')
#p.set(xlabel='Length')
#p._legend.set_title('Non-B DNA structures')
#plt.tight_layout()
#plt.savefig('/home/alextu/scratch/visualize_basecall/imgs/motif_length_density_all_non_b_structures_perchr_l100.png')