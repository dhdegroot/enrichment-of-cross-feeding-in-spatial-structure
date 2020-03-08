import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

SAVE_FIGURES = True

working_dir = os.getcwd()
G_df = pd.read_csv(os.path.join(working_dir, "results", "saved_results", "2000_simus_G_df.csv"))

n_simus = int(max(G_df.simulation.values)) + 1

# Initialize empty dataframe in which we will store all simulation factor
results_df = pd.DataFrame(
    columns=['simulation', 'freqa', 'freqb', 'freqc', 'CC0', 'CCC_ind', 'CCA_ind', 'CCB_ind', 'adv_cheat',
             'max_benefit', 'max_benefit_lambda', 'range_lambda_min', 'range_lambda_max', 'fold_diff_ab'])

# Run over all simulations
for simu_ind in range(n_simus):
    G_df_simu = G_df[(G_df.simulation == simu_ind) & (G_df.species == 'coop_advantage')]

    # Determine maximal growth benefit
    max_benefit_ind = np.argmax(G_df_simu.G.values)
    max_benefit = G_df_simu.G.iloc[max_benefit_ind]

    # Determine corresponding avg_nt
    max_benefit_avg_nt = G_df_simu['lambda'].iloc[max_benefit_ind]

    # Determine range of avg_nt for which cooperator benefit is positive
    pos_avg_nt = G_df_simu[G_df_simu.G >= 0]['lambda']
    min_pos_avg_nt = min(pos_avg_nt)
    max_pos_avg_nt = max(pos_avg_nt)

    # Determine absolute fold change difference from a to b
    fold_diff_ab = np.abs(np.log2(G_df_simu['freqa'].iloc[0] / G_df_simu['freqb'].iloc[0]))

    # Store everything in results_df with simulation parameters
    results_dict_simu = {'max_benefit': max_benefit, 'max_benefit_lambda': max_benefit_avg_nt,
                         'range_lambda_min': min_pos_avg_nt, 'range_lambda_max': max_pos_avg_nt,
                         'fold_diff_ab': fold_diff_ab}

    for col in results_df.columns:
        if col in G_df_simu.columns:
            results_dict_simu[col] = G_df_simu[col].iloc[0]

    results_df = results_df.append(results_dict_simu, ignore_index=True)

# Make scatterplot with initial fraction of cheaters on x-axis, the fold change diff between a b on y-axis,
# and the growth benefit as colour
plt.figure(1)
ax = sns.scatterplot(x='freqc', hue='fold_diff_ab', y='max_benefit', data=results_df)
if SAVE_FIGURES:
    plt.savefig(os.path.join(working_dir, "results", "cheater_fraction_explains_max_benefit.png"))

# Make plot in which you plot max and min of positive lambda range as a function of max benefit
results_df_small = results_df[['range_lambda_min', 'range_lambda_max', 'max_benefit']]
results_df_small_tidy = pd.melt(results_df_small, ['max_benefit'], var_name='min_max', value_name='lambda')
plt.figure(2)
ax2 = sns.lineplot(x='max_benefit', y='lambda', hue='min_max', data=results_df_small_tidy)
if SAVE_FIGURES:
    plt.savefig(os.path.join(working_dir, "results", "max_benefit_vs_range_lambda"))

# Make plot in which you plot max_benefit, max_benefit_lambda
plt.figure(3)
ax3 = sns.scatterplot(x='max_benefit', y='max_benefit_lambda', data=results_df)
if SAVE_FIGURES:
    plt.savefig(os.path.join(working_dir, "results", "max_benefit_vs_max_benefit_lambda"))

plt.show()
pass
