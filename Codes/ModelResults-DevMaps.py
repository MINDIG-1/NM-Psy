import sys
import os

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from Codes.Utils import *

stage = "ModelResults"

result_path = os.path.join(result_root, stage, "DeviationMaps")
os.makedirs(result_path, exist_ok=True)

atleastOneED = []
EDsubj_ch_list = []
chPsubj_list = []
chPlobe_list = []
atleastOneEDlobe = []
for folder in (zs_folders):
    
    zs_path = f'{zs_path_root}/{folder}/{zs_file}'
    zs_df = pd.read_csv(zs_path)
    zs_df = zs_df[zs_df['group'] != "HC_tr"]
    fig, axes = plt.subplots(nrows=5, ncols=2, figsize=(8, 24))
    os.makedirs(f'{result_path}/{folder}', exist_ok=True)


    for i, group in enumerate(group_list):

        dataset = zs_df[zs_df['group']==group]['dataset'].values
        IO_pos,IO_neg = zs_IO(zs_df, column_names, group)
        temp_pos = round((IO_pos.sum(axis=1) != 0).sum() / IO_pos.shape[0] * 100, 2)
        temp_neg = round((IO_neg.sum(axis=1) != 0).sum() / IO_neg.shape[0] * 100, 2)

        # print((posnegED_group.sum(axis=1)==0).sum(),posnegED_group.shape[0] )
 
        atleastOneED.append(pd.DataFrame({"group": [group], "f_band": [folder[6:]], "pos_atleastOne": [temp_pos], "pos_range": f'{np.median(IO_pos.sum(axis=1), axis=0)} [{IO_pos.sum(axis=1).min()}-{IO_pos.sum(axis=1).max()}]', "neg_atleastOne": [temp_neg],  "neg_range": f'{np.median(IO_neg.sum(axis=1), axis=0)} [{IO_neg.sum(axis=1).min()}-{IO_neg.sum(axis=1).max()}]'}))

        # pos_subj_ch, pos_ch_subj, neg_subj_ch, neg_ch_subj, id, dataset = zs_counts(zs_df, column_names, group, f'{result_path}/{folder}/{group}/')
        pos_subj_ch, pos_ch_subj = zs_counts(IO_pos) , IO_pos.sum(axis=1)
        neg_subj_ch, neg_ch_subj = zs_counts(IO_neg) , IO_neg.sum(axis=1)

        temp_df = pd.DataFrame(pos_subj_ch).T
        temp_df.columns = [f"zs_X{i}" for i in range(1, len(temp_df.columns) + 1)]
        temp_df = temp_df.assign(group=group, f_band=folder[6:], posneg = 'pos')
        EDsubj_ch_list.append(temp_df)
        
        temp_df = pd.DataFrame(neg_subj_ch).T
        temp_df.columns = [f"zs_X{i}" for i in range(1, len(temp_df.columns) + 1)]
        temp_df = temp_df.assign(group=group, f_band=folder[6:], posneg = 'neg')
        EDsubj_ch_list.append(temp_df)
        
        temp_chPsubj = pd.DataFrame({'global_id':id,'posChED':pos_ch_subj, 'negChED':neg_ch_subj, 'group': group, 'dataset':dataset, "f_band": folder[6:]})

        chPsubj_list.append(temp_chPsubj)

        # print(max(vmax, pos_subj_ch.max(), neg_subj_ch.max()))
        print(f'max ch + : {round(pos_subj_ch.max(),2)}, max ch - {round(neg_subj_ch.max(),2)}, group: {group}, fb:{folder[6:]}')
        max_p = max(round(pos_subj_ch.max(),2), max_p)
        max_n = max(round(neg_subj_ch.max(),2), max_n)

        
        im1, cm1 = mne.viz.plot_topomap(pos_subj_ch, pos=info, axes=axes[i, 0], show=False, vlim=[vmin, vmax], cmap='Blues', contours = 0)
        divider1 = make_axes_locatable(axes[i, 0])
        cax1 = divider1.append_axes('right', size='5%')
        cbar1 = plt.colorbar(im1, cax=cax1)
        axes[i, 0].set_title(f"% > 2 {[group]}-{folder[5:]}")

        im2, cm2 = mne.viz.plot_topomap(neg_subj_ch, pos=info, axes=axes[i, 1], show=False, vlim=[vmin, vmax], cmap='Reds', contours = 0)
        divider2 = make_axes_locatable(axes[i, 1])
        cax2 = divider2.append_axes('right', size='5%')
        cbar2 = plt.colorbar(im2, cax=cax2)
        axes[i, 1].set_title(f"% < -2 {[group]}-{folder[5:]}")
    
    fig.savefig(f'{result_path}/{folder}/topo.png')
    plt.close(fig)


print(f'max_p= {max_p}, max_n={max_n}')
df_atleastchOneED = pd.concat(atleastOneED, ignore_index=True)
df_atleastchOneED.to_csv(os.path.join(result_path, 'atleastOne.csv'), index = False)

df_EDsubj_ch = pd.concat(EDsubj_ch_list, ignore_index=True)
# print(df_EDsubj_ch[df_EDsubj_ch['group']=="HC_tr"])
df_EDsubj_ch.to_csv(os.path.join(result_path, 'EDsubjectsPch.csv'), index = False)

chPsubj = pd.concat(chPsubj_list)
chPsubj.to_csv(f'{result_path}/EDchPsubj.csv')

fig, axes = plt.subplots(nrows=5, ncols=1, figsize=(15, 36))
df = chPsubj
df = df[df['group'] != "HC_tr"]

# grouped_df = df.groupby('f_band')
grouped_df = df.groupby(pd.Categorical(df['f_band'], categories=f_bands, ordered=True))

# Define colors for positive and negative data



for (f_band, group), ax in zip(grouped_df, axes.flatten()):
# for (f_band, group) in grouped_df:

    # fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(15, 10))
    melted_data = pd.melt(group[['posChED', 'negChED', 'group']].reset_index(),
                         id_vars=['index', 'group'], var_name='variable', value_name='value')
    
    sns.violinplot(x='group', y='value', hue='variable', split=False, inner='quart', cut = 0, palette=colors,
                   data=melted_data, ax=ax, alpha = 0.3)

    strip = sns.stripplot(x='group', y='value', hue='variable', data=melted_data, palette=colors, alpha=0.5, jitter=True, dodge=True, ax=ax)
    

    print(f'***{f_band}')
    add_significance_stars(group, df[df['group'] == 'HC_te'], ax)


    # Customize legend labels
    handles, labels = strip.get_legend_handles_labels()
    labels = ['+', '-']
    ax.legend(handles, labels)
    
    ax.set_title(f'Frequency Band: {f_band}')
    ax.set_xlabel('Group')
    ax.set_ylabel('Count')

# Adjust layout
plt.tight_layout(rect=[0, 0, 1, 0.97])
plt.savefig(f'{result_path}/violinbChPsubj.png', dpi = 300)