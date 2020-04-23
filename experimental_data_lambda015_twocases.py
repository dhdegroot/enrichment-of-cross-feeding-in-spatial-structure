from helpers_coop_droplets import *

"""Constants"""
working_dir = os.getcwd()
spec_names = ['a', 'b', 'c']
spec_names_caps = ['A', 'B', 'C']
model_types = ['dependent', 'indep. costly']
SAVE_LAST_FIGURE = True
SAVE_FIGURES = False
MAKE_EACH_SIMU_FIGURE = False
PRINT_ALOT = False
MAKE_PLOTS = False
N_ITERATIONS = 4

"""Get parameters"""
# These parameters can be changed in the helpers_coop_droplets.py-file
freqs, CC0, CCA_ind, CCB_ind, CCC_ind, adv_cheat = get_starting_parameters()

simu_dict_list = []
avg_nt = .15  # This nt was picked because it was the lambda that gave the largest benefit to coops on average

CCAind_types = [7.5, 150]
CCBind_types = [0, 0]
CCCind_types = [7.5, 150]

start_freqs = [np.array([.41, .09, .5]), np.array([0.35, 0.43, 0.22])]

"""Generate set of parameters for all simulations"""
for i in range(2):
    simu_dict_list.append(
        {'freqs': start_freqs[i], 'CC0': CC0, 'CCC_ind': CCCind_types[i], 'CCA_ind': CCAind_types[i],
         'CCB_ind': CCBind_types[i], 'adv_cheat': adv_cheat, 'model_type': model_types[i]})

"""Calculate for different cheater starter fractions, if the cooperators or the cheaters will win"""
end_freq_df = pd.DataFrame(
    columns=['lambda', 'simulation', 'species', 'end_freq', 'iteration', 'model_type'])

for simu_ind, simu_dict in enumerate(simu_dict_list):
    end_freq_df_iter = pd.DataFrame(columns=['lambda', 'simulation', 'species', 'end_freq', 'iteration'])

    # Add as first iter the startpoint
    end_freq_df_iter = end_freq_df_iter.append(
        pd.DataFrame.from_dict(
            {'lambda': avg_nt, 'simulation': simu_ind, 'species': spec_names_caps, 'end_freq': simu_dict['freqs'],
             'iteration': 0}))

    print("Starting simulation '{0}'".format(simu_ind))
    freqs_simu = simu_dict['freqs'].copy()
    CCC_ind = simu_dict['CCC_ind']
    CCA_ind = simu_dict['CCA_ind']
    CCB_ind = simu_dict['CCB_ind']
    iter_dict_list = []
    for iter_ind in range(N_ITERATIONS):
        if iter_ind % 5 == 4:
            print("Computing '{0}'th iteration".format(iter_ind + 1))
        iter_dict = {'freqs': freqs_simu.copy(), 'CC0': CC0, 'CCC_ind': CCC_ind, 'CCA_ind': CCA_ind, 'CCB_ind': CCB_ind,
                     'adv_cheat': adv_cheat}
        iter_dict_list.append(iter_dict)

        # Calculate growth dataframe for one simulation
        G_df_simu_tidy, end_freq_df_simu_tidy = calc_for_one_simu(simu_ind, iter_dict, avg_nt, spec_names=spec_names,
                                                                  MAKE_PLOTS=MAKE_PLOTS,
                                                                  PRINT_ALOT=PRINT_ALOT, SAVE_FIGURES=SAVE_FIGURES,
                                                                  MAKE_EACH_SIMU_FIGURE=MAKE_EACH_SIMU_FIGURE)
        end_freq_df_simu_tidy['iteration'] = iter_ind + 1
        end_freq_df_iter = end_freq_df_iter.append(end_freq_df_simu_tidy)

        for spec_ind, spec in enumerate(spec_names_caps):
            freqs_simu[spec_ind] = end_freq_df_simu_tidy[end_freq_df_simu_tidy.species == spec].end_freq.values[0]

    end_freq_df_iter['model_type'] = simu_dict['model_type']

    end_freq_df = end_freq_df.append(end_freq_df_iter)

"""Store end_freq_df"""
end_freq_df.to_csv(os.path.join(working_dir, "results", "experimental_simulations_end_freq_df_lambda015_twocases.csv"),
                   index=False, header=True)

"""Make stacked bar plots, one for each simulation"""
# Set general plot properties
sns.set_style("white")
colors = sns.color_palette("muted", n_colors=3)
g = sns.FacetGrid(end_freq_df, col='model_type', legend_out=True)
g = g.map_dataframe(mullerplot_fig2, 'lambda', 'end_freq', new_colors=colors)

# Plot 1 - background - "total" (top) series
# g = g.map_dataframe(plot_stacked_bars, 'lambda', 'end_freq', species='sum_all', new_color=colors[0])
# Plot 2 - overlay - "bottom" series
# g = g.map_dataframe(plot_stacked_bars, 'lambda', 'end_freq', species='sum_coops', new_color=colors[1])
# Plot 3
# g = g.map_dataframe(plot_stacked_bars, 'lambda', 'end_freq', species='A', new_color=colors[2])
g.set_titles('')
g.set_axis_labels('Number of propagations', '')
sns.despine(left=True)

topbar = plt.Rectangle((0, 0), 1, 1, fc=colors[0], edgecolor='none')
middle_bar = plt.Rectangle((0, 0), 1, 1, fc=colors[1], edgecolor='none')
bottombar = plt.Rectangle((0, 0), 1, 1, fc=colors[2], edgecolor='none')
g.add_legend({'Coop A': topbar, 'Coop B': middle_bar, 'Cheater': bottombar})

plt.savefig(os.path.join(working_dir, "results", "experimental_simulation_lambda015_twocases.png"))

