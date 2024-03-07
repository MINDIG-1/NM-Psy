import os
import pandas as pd
import numpy as np
import mne
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.stats import mannwhitneyu
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.cluster import AgglomerativeClustering
from scipy.spatial.distance import pdist, squareform
from sklearn.metrics import pairwise_distances
from sklearn.cluster import AgglomerativeClustering
from matplotlib import rcParams
from scipy.stats import mannwhitneyu
import matplotlib as mpl
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
from statsmodels.stats.multitest import multipletests 
import random
import joblib

mpl.rcParams['xtick.labelsize'] = 16
mpl.rcParams['ytick.labelsize'] = 16
mpl.rcParams['axes.labelsize'] = 20
mpl.rcParams['axes.titlesize'] = 20


feature = "Spec"
result_root = "/home/ubuntu/NM-Psy/Results"

zs_path_root = "home/ubuntu/NM-Psy/Results/ModelTraining/Zs-Parallel"

zs_folders = [file for file in os.listdir(zs_path_root) if file.startswith("Allch")]

column_names = [f"zs_X{i}" for i in range(1, 125)]

montage = mne.channels.make_standard_montage('GSN-HydroCel-128')
info = mne.create_info(montage.ch_names[0:124], sfreq=200, ch_types='eeg')
info.set_montage(montage)

vmin = 0
vmax = 100
max_p = float('-inf')
max_n = float('-inf')
group_list = [ "HC_te", "ADHD", "ASD", "Anxiety", "Learning"]
zs_file = "zsAllCh.csv"

f_bands = ['Delta', 'Theta', 'Alpha', 'Beta', 'Gamma']
colors = {'posChED': '#0A2B6B', 'negChED': '#730128'}


nperm = 5000
perm_test = 'group' # 'group' or 'spatial'
level = 'spectral' # 'network' or 'region' or 'spectral'
sig_thresh = 0.05
plot_name = ["uncorrected", "fdr"]




def zs_IO(zs_df, column_names, group):
    zs_subset = zs_df[zs_df['group']==group][column_names].values
    IO_pos = np.where(zs_subset > 2, 1, 0)
    IO_neg = np.where(zs_subset < -2, 1, 0)
    return IO_pos, IO_neg

def zs_counts (IO):
    IO = IO[IO.sum(axis=1) != 0]
    if IO.shape[0] != 0:
        l = 100 / IO.shape[0] #subjects
    else:
        l = 0

    subj_ch = IO.sum(axis = 0)*l

    return subj_ch


def add_significance_stars(group_data, reference_group, ax):
    for variable in ['posChED', 'negChED']:
        if variable in group_data.columns and variable in reference_group.columns:
            for group_name in group_data['group'].unique():
                if group_name != 'HC(test)':
                    _, p_value = mannwhitneyu(reference_group[variable], group_data[group_data['group'] == group_name][variable])
                    if p_value < 0.05:
                        print(group_name, variable, p_value)
                        x_position = group_data['group'].unique().tolist().index(group_name)  # Get the x-coordinate for the group
                        if variable=='posChED':
                            x_position -= 0.2
                        else:
                            x_position += 0.2
              
                        y_position = group_data[variable].min() - 4  # Adjust the y-coordinate for better visibility

                        if p_value > 0.01:
                            ax.text(x_position, y_position, '*', ha='center', va='center', fontsize=20, fontweight='bold')
                        else:
                            ax.text(x_position, y_position, '**', ha='center', va='center', fontsize=20, fontweight='bold')
                        _, p_value_g = mannwhitneyu(reference_group[variable], group_data[group_data['group'] == group_name][variable], alternative = 'greater' )

                        if p_value_g < 0.05:
                            print(f'HC is greater than {group_name}  at {variable}')
                        _, p_value_l = mannwhitneyu(reference_group[variable], group_data[group_data['group'] == group_name][variable], alternative = 'less' )
                        if p_value_l < 0.05:
                            print(f'HC is less than {group_name} at {variable}')



def permutation_test(ip, zs_df, column_names, group_list, folder, col_dict, perm_test):
    fb = folder.split('-')[1]
    print(f'Band {fb} --- Permutation {ip}')
    
    nullval_pos_dict = {group: [] for group in group_list}
    nullval_neg_dict = {group: [] for group in group_list}

    if perm_test == 'group':
        labels_shuffled = zs_df['group'].sample(frac=1).reset_index(drop=True)
        zs_df_shuffled = zs_df.copy()
        zs_df_shuffled['group'] = np.array(labels_shuffled)
    elif perm_test == 'spatial':
        zs_df_shuffled = zs_df.copy()
        zs_df_shuffled[column_names] = np.random.permutation(zs_df[column_names].values.T).T

    IO_pos_HC, IO_neg_HC = zs_IO(zs_df_shuffled, column_names, "HC_te")
    HC_temp_pos = zs_counts(IO_pos_HC)
    HC_temp_neg = zs_counts(IO_neg_HC)

    for i, group in enumerate(group_list):
        IO_pos, IO_neg = zs_IO(zs_df_shuffled, column_names, group)

        if level == 'region' or level == 'spectral':
            grp_temp = zs_counts(IO_pos)

        grp_temp = grp_temp - HC_temp_pos
        nullval_pos_dict[group].append(grp_temp)

        if level == 'region' or level == 'spectral':
            grp_temp = zs_counts(IO_neg)

        grp_temp = grp_temp - HC_temp_neg
        nullval_neg_dict[group].append(grp_temp.reshape(1, -1))

    return nullval_pos_dict, nullval_neg_dict


def apply_fdr_correction(p_values):
    _, p_values_corrected, _, _ = multipletests(p_values, method='fdr_bh')
    return p_values_corrected

def find_rank(sorted_arr_pos, true_val_pos):
    indices = []

    for col_idx in range(sorted_arr_pos.shape[1]):

        row_idx = np.argmax(sorted_arr_pos[:, col_idx] == true_val_pos[col_idx])
        indices.append(row_idx)
    return indices
