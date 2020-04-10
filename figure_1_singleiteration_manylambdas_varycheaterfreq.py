from helpers_coop_droplets import *

"""Constants"""
SAVE_LAST_FIGURE = True
SAVE_FIGURES = False
MAKE_EACH_SIMU_FIGURE = False
PRINT_ALOT = False
MAKE_PLOTS = True
working_dir = os.getcwd()
N_SIMUS = 5 * 2 + 1
np.random.seed(3141)
spec_names = ['a', 'b', 'c']

# Get parameters. These can be changed in helpers_coop_droplets file
freqs, CC0, CCA_ind, CCB_ind, CCC_ind, adv_cheat = get_starting_parameters()

first_avg_nt = .05
last_avg_nt = 4
n_avg_nt = 40

# sims_df = pd.DataFrame(columns=['freqs', 'CC0', 'CCC_ind'])
simu_dict_list = []

"""Generate set of parameters for all simulations"""
cheater_freqs = np.linspace(0, 1, (N_SIMUS + 1) / 2)
for ind in range(int(np.ceil((N_SIMUS+1) / 2))):
    freqs_simu_first = np.array([(1 - cheater_freqs[ind]) / 2, (1 - cheater_freqs[ind]) / 2, cheater_freqs[ind]])
    freqs_simu_second = np.array([4 * (1 - cheater_freqs[ind]) / 5, (1 - cheater_freqs[ind]) / 5, cheater_freqs[ind]])

    simu_dict_list.append(
        {'freqs': freqs_simu_first, 'CC0': CC0, 'CCC_ind': CCC_ind, 'CCA_ind': CCA_ind,
         'CCB_ind': CCB_ind, 'adv_cheat': adv_cheat})
    simu_dict_list.append(
        {'freqs': freqs_simu_second, 'CC0': CC0, 'CCC_ind': CCC_ind, 'CCA_ind': CCA_ind,
         'CCB_ind': CCB_ind, 'adv_cheat': adv_cheat})

simu_dict_list.pop(-1)  # Remove first

# avg_nt_range = np.exp(np.linspace(np.log(first_avg_nt), np.log(last_avg_nt), n_avg_nt))
avg_nt_range = np.linspace(first_avg_nt, last_avg_nt, n_avg_nt)
if avg_nt_range.size > 2:
    MAKE_PLOTS = False
    SAVE_FIGURES = False

# Initialize dataframe in which we are going to store growth factor data
G_df = pd.DataFrame(columns=['lambda', 'simulation', 'species', 'G'])
end_freq_df = pd.DataFrame(columns=['lambda', 'simulation', 'species', 'end_freq'])

# Run all different simulations
for simu_ind in range(N_SIMUS):
    print("Starting simulation '{0}'".format(simu_ind))
    simu_dict = simu_dict_list[simu_ind]

    # Calculate growth dataframe for one simulation
    G_df_simu_tidy, end_freq_df_simu_tidy = calc_for_one_simu(simu_ind, simu_dict, avg_nt_range, spec_names=spec_names,
                                                              MAKE_PLOTS=MAKE_PLOTS,
                                                              PRINT_ALOT=PRINT_ALOT, SAVE_FIGURES=SAVE_FIGURES,
                                                              MAKE_EACH_SIMU_FIGURE=MAKE_EACH_SIMU_FIGURE)
    G_df = G_df.append(G_df_simu_tidy)
    end_freq_df = end_freq_df.append(end_freq_df_simu_tidy)

"""Store end_freq_df"""
start_freq_df = pd.DataFrame(columns=['lambda', 'simulation', 'species', 'end_freq'])

for simu_ind in range(N_SIMUS):
    simu_dict = simu_dict_list[simu_ind]['freqs']
    start_freq_df = start_freq_df.append(
        {'lambda': 0, 'simulation': simu_ind, 'species': 'A', 'end_freq': simu_dict[0]},
        ignore_index=True)
    start_freq_df = start_freq_df.append(
        {'lambda': 0, 'simulation': simu_ind, 'species': 'B', 'end_freq': simu_dict[1]},
        ignore_index=True)
    start_freq_df = start_freq_df.append(
        {'lambda': 0, 'simulation': simu_ind, 'species': 'C', 'end_freq': simu_dict[2]},
        ignore_index=True)

end_freq_df = pd.concat([end_freq_df, start_freq_df])

end_freq_df.to_csv(os.path.join(working_dir, "results", "end_freq_df_singlepropagation_manylambdas.csv"), index=False, header=True)
end_freq_df['lambda'] = end_freq_df['lambda'].astype(object)

"""Make stacked bar plots, one for each simulation"""
# Set general plot properties
sns.set_style("white")
colors = sns.color_palette("muted", n_colors=3)
g = sns.FacetGrid(end_freq_df, row="simulation", legend_out=True, height=1, aspect=N_SIMUS)
g = g.map_dataframe(mullerplot, 'lambda', 'end_freq', new_colors=colors)

# Plot 1 - background - "total" (top) series
# g = g.map_dataframe(plot_stacked_bars, 'lambda', 'end_freq', species='sum_all', new_color=colors[0])
# Plot 2 - overlay - "bottom" series
# g = g.map_dataframe(plot_stacked_bars, 'lambda', 'end_freq', species='sum_coops', new_color=colors[1])
# Plot 3
# g = g.map_dataframe(plot_stacked_bars, 'lambda', 'end_freq', species='A', new_color=colors[2])
g.set_titles('')
g.set_axis_labels('Avg number of cells per compartment', '')
sns.despine(left=True)

topbar = plt.Rectangle((0, 0), 1, 1, fc=colors[0], edgecolor='none')
middle_bar = plt.Rectangle((0, 0), 1, 1, fc=colors[1], edgecolor='none')
bottombar = plt.Rectangle((0, 0), 1, 1, fc=colors[2], edgecolor='none')
g.add_legend({'Coop A': topbar, 'Coop B': middle_bar, 'Cheater': bottombar})

plt.savefig(os.path.join(working_dir, "results", "Figure1_singleiteration_manylambda.png"))
