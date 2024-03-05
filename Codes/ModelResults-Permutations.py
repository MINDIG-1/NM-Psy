import sys
import os

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from Codes.Utils import *

stage = "ModelResults"
devmaps_path = os.path.join(result_root, stage, "DeviationMaps")
result_path = os.path.join(result_root, stage, "Permutations") 


if  level == 'region':
    col_dict = {f"zs_X{i}": f"zs_X{i}" for i in range(0, len(column_names))}
    col_list = [f"zs_X{i}" for i in range(1, len(column_names)+1)]
    true_net = pd.read_csv(f'{devmaps_path}/EDsubjectsPch.csv')
    true_net_pos = true_net[true_net['posneg']=='pos']
    true_net_neg = true_net[true_net['posneg']=='neg']
    true_net_pos = true_net_pos[column_names + ['group', 'f_band']]
    true_net_neg = true_net_neg[column_names + ['group', 'f_band']]
elif level == 'spectral':
    col_dict = {f"zs_X{i}": f"zs_X{i}" for i in range(1, len(column_names)+1)}
    col_list = [f"zs_X{i}" for i in range(1, len(column_names)+1)]
    true_net = pd.read_csv(f'{devmaps_path}/EDsubjectsPch.csv')
    true_net_pos = true_net[true_net['posneg']=='pos']
    true_net_neg = true_net[true_net['posneg']=='neg']
    true_net_pos = true_net_pos[column_names + ['group', 'f_band']]
    true_net_neg = true_net_neg[column_names + ['group', 'f_band']]


        


# zs_folders = ['Allch-Gamma']
for folder in zs_folders:
    
    # start_time = time.time()
    # print(start_time)

    fb = folder.split('-')[1]


    nullval_pos_dict = {gp: [] for gp in group_list}
    nullval_neg_dict = {gp: [] for gp in group_list}
    all_nullval_pos = []
    all_nullval_neg = []

    df_pval_pos_all = pd.DataFrame(columns=list(col_dict.keys())+['group', 'f_band'])
    df_pval_neg_all = pd.DataFrame(columns=list(col_dict.keys())+['group', 'f_band'])
    df_pval_pos_fdr_all = pd.DataFrame(columns=list(col_dict.keys())+['group', 'f_band'])
    df_pval_neg_fdr_all = pd.DataFrame(columns=list(col_dict.keys())+['group', 'f_band'])

    true_net_pos_fb = np.array(true_net_pos[true_net_pos['f_band']==fb])
    true_net_neg_fb = np.array(true_net_neg[true_net_neg['f_band']==fb])


    true_HC_pos = true_net_pos_fb[np.where(true_net_pos_fb == "HC_te")[0],:]
    true_HC_neg = true_net_neg_fb[np.where(true_net_neg_fb == "HC_te")[0],:]


    true_net_pos_fb[:, :-2] = true_net_pos_fb[:, :-2] - true_HC_pos[:, :-2]
    true_net_neg_fb[:, :-2] = true_net_neg_fb[:, :-2] - true_HC_neg[:, :-2]


    zs_path = f'{zs_path_root}/{folder}/{zs_file}'
    zs_df = pd.read_csv(zs_path, index_col=False)
    zs_df = zs_df[zs_df['group']!= "HC_tr"]
 

    parallel_result = joblib.Parallel(n_jobs=-1)(
        joblib.delayed(permutation_test)(ip, zs_df, column_names, group_list, folder, col_dict, perm_test)
        for ip in range(2, nperm + 1)
    )

    all_nullval_pos.extend(result[0] for result in parallel_result)
    all_nullval_neg.extend(result[1] for result in parallel_result)


    for group in group_list:
        nullval_pos_dict[group] = [d[group][0].reshape(1, -1) for d in all_nullval_pos]

        nullval_pos_dict[group].append(true_net_pos_fb[np.where(true_net_pos_fb == group)[0],:-2])
        nullval_neg_dict[group] = [d[group][0] for d in all_nullval_neg]

        nullval_neg_dict[group].append(true_net_neg_fb[np.where(true_net_neg_fb == group)[0],:-2])

        df_nullval_pos = np.concatenate(nullval_pos_dict[group])

        df_nullval_neg = np.concatenate(nullval_neg_dict[group])

        df_pval_pos = pd.DataFrame(columns=list(col_dict.keys())+['group', 'f_band'])
        df_pval_neg = pd.DataFrame(columns=list(col_dict.keys())+['group', 'f_band'])
        df_pval_pos.loc[0, 'group'], df_pval_neg.loc[0, 'group'] = group, group
        df_pval_pos.loc[0, 'f_band'], df_pval_neg.loc[0, 'f_band'] = fb, fb

        
        true_val_pos, true_val_neg = true_net_pos_fb[np.where(true_net_pos_fb == group)[0],:-2][0], true_net_neg_fb[np.where(true_net_neg_fb == group)[0],:-2][0]
        sorted_arr_pos_des, sorted_arr_neg_des = np.sort(df_nullval_pos, axis=0)[::-1],  np.sort(df_nullval_neg, axis=0)[::-1]
        sorted_arr_pos_asc, sorted_arr_neg_asc = np.sort(df_nullval_pos, axis=0)[::1],  np.sort(df_nullval_neg, axis=0)[::1]
        
        rank_pos_des = find_rank(sorted_arr_pos_des,true_val_pos)
        rank_neg_des = find_rank(sorted_arr_neg_des,true_val_neg)
        rank_pos_asc = find_rank(sorted_arr_pos_asc,true_val_pos)
        rank_neg_asc = find_rank(sorted_arr_neg_asc,true_val_neg)
        min_rank_pos = []
        min_rank_neg = []

        # Compare ranks for each column and keep the lower rank
        for col_index in range(sorted_arr_pos_des.shape[1]):
            min_rank_pos.append(min(rank_pos_des[col_index], rank_pos_asc[col_index]))
            min_rank_neg.append(min(rank_neg_des[col_index], rank_neg_asc[col_index]))

        min_rank_pos = np.array(min_rank_pos)
        min_rank_neg = np.array(min_rank_neg)

        pval_pos_list = [round(value / nperm, 4) for value in min_rank_pos]
        pval_neg_list = [round(value / nperm, 4) for value in min_rank_neg]
        df_pval_pos.iloc[0,:-2] = pval_pos_list
        df_pval_neg.iloc[0,:-2] = pval_neg_list

        print(f"{group} - pos: {np.sum((np.array(pval_pos_list) <0.05))} - neg: {np.sum(np.array(pval_neg_list) <0.05)}")

        df_pval_pos_all = pd.concat([df_pval_pos_all, df_pval_pos], ignore_index=True)
        df_pval_neg_all = pd.concat([df_pval_neg_all, df_pval_neg], ignore_index=True)
        
        # Apply FDR correction to p-values
        df_pval_pos_fdr = df_pval_pos.copy()
        df_pval_neg_fdr = df_pval_neg.copy()

        df_pval_pos_fdr[column_names] = apply_fdr_correction(pval_pos_list)

        df_pval_neg_fdr[column_names] = apply_fdr_correction(pval_neg_list)
        df_pval_pos_fdr_all = pd.concat([df_pval_pos_fdr_all, df_pval_pos_fdr], ignore_index=True)
        df_pval_neg_fdr_all = pd.concat([df_pval_neg_fdr_all, df_pval_neg_fdr], ignore_index=True)

    if feature == "Spec":
        fig, axes = plt.subplots(nrows=5, ncols=2, figsize=(8, 24))

        fig.suptitle('p_vale < 0.05 Group Permutation - Spectral Features ', fontsize=16)

        for j, sublist in enumerate([[df_pval_pos_all, df_pval_neg_all], [df_pval_pos_fdr_all, df_pval_neg_fdr_all]]):

            sig_df_pos, sig_df_neg = sublist

            for i, group in enumerate(group_list):

                pos_df = sig_df_pos[sig_df_pos['group']==group]
                df_bin_p = (pos_df.iloc[:, :-2] <sig_thresh).astype(int)
                pos = df_bin_p.to_numpy().reshape(-1)
                im1, cm1 = mne.viz.plot_topomap(pos, pos=info, axes=axes[i, 0],vlim=[0, 1], show=False,  cmap='Blues', contours=0)
                divider1 = make_axes_locatable(axes[i, 0])
                axes[i, 0].set_title(f" +ED {[group]}-{folder[5:]}")

                neg_df = sig_df_neg[sig_df_neg['group']==group]
                df_bin_n = (neg_df.iloc[:, :-2] <sig_thresh).astype(int)
                neg = df_bin_n.to_numpy().reshape(-1)

                im1, cm1 = mne.viz.plot_topomap(neg, pos=info, axes=axes[i, 1],vlim=[0, 1], show=False, cmap='PuRd', contours=0)
                divider1 = make_axes_locatable(axes[i, 1])
                axes[i, 1].set_title(f"-ED {[group]}-{folder[5:]}")
            os.makedirs(f'{result_path}/SpecPermutationPlots/{folder}', exist_ok=True)
            fig.savefig(f'{result_path}/SpecPermutationPlots/{folder}/pvalgroup{plot_name[j]}.png')




    pval_save = os.path.join(result_path, folder)
    if not os.path.exists(pval_save):
        os.makedirs(pval_save)
        
    df_pval_pos_all.to_csv(os.path.join(pval_save, f'pval_{level}_{perm_test}perm_pos.csv'), index = False)
    df_pval_neg_all.to_csv(os.path.join(pval_save, f'pval_{level}_{perm_test}perm_neg.csv'), index = False)

    df_pval_pos_fdr_all.to_csv(os.path.join(pval_save, f'pval_{level}_{perm_test}perm_pos_fdr.csv'), index = False)
    df_pval_neg_fdr_all.to_csv(os.path.join(pval_save, f'pval_{level}_{perm_test}perm_neg_fdr.csv'), index = False)