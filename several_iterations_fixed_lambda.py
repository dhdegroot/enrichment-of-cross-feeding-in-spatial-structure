from helpers_coop_droplets import *

"""Constants"""
working_dir = os.getcwd()
spec_names = ['a', 'b', 'c']
spec_names_caps = ['A', 'B', 'C']
SAVE_LAST_FIGURE = True
SAVE_FIGURES = False
MAKE_EACH_SIMU_FIGURE = False
PRINT_ALOT = False
MAKE_PLOTS = False
N_ITERATIONS = 100

"""Get parameters"""
# These parameters can be changed in the helpers_coop_droplets.py-file
freqs, CC0, CCA_ind, CCB_ind, CCC_ind, adv_cheat = get_starting_parameters()

avg_nt = 1.811  # This nt was picked because it was the lambda that gave the largest benefit to coops on average
start_freqs = np.array([0.5, 0.01, 0.49])

"""Calculate for different cheater starter fractions, if the cooperators or the cheaters will win"""
simu_dict_list = []
G_df = pd.DataFrame(columns=['lambda', 'simulation', 'species', 'G'])

freqs_simu = start_freqs
for iter_ind in range(N_ITERATIONS):
    simu_dict = {'freqs': freqs_simu.copy(), 'CC0': CC0, 'CCC_ind': CCC_ind, 'CCA_ind': CCA_ind, 'CCB_ind': CCB_ind,
                 'adv_cheat': adv_cheat}
    simu_dict_list.append(simu_dict)

    print("Starting simulation '{0}'".format(iter_ind))

    # Calculate growth dataframe for one simulation
    G_df_simu_tidy = calc_for_one_simu(iter_ind, simu_dict, avg_nt, spec_names=spec_names, MAKE_PLOTS=MAKE_PLOTS,
                                       PRINT_ALOT=PRINT_ALOT, SAVE_FIGURES=SAVE_FIGURES,
                                       MAKE_EACH_SIMU_FIGURE=MAKE_EACH_SIMU_FIGURE)
    G_df = G_df.append(G_df_simu_tidy)
    for spec_ind, spec in enumerate(spec_names_caps):
        growth_factor = G_df_simu_tidy[G_df_simu_tidy.species == spec].G.values[0]
        freqs_simu[spec_ind] *= growth_factor

    freqs_simu = freqs_simu/np.sum(freqs_simu)


"""At this point, we are left with G_df in which the average growth rate curves for the different simulations are shown,
and with a simu_dict_list where we have stored the parameters for each simulation"""
G_df['start_freq'] = 0
G_df['CC0'] = G_df['CCC_ind'] = G_df['CCA_ind'] = G_df['CCB_ind'] = G_df['adv_cheat'] = 0
for simu_ind in range(N_ITERATIONS):
    simu_dict = simu_dict_list[simu_ind]
    for key in simu_dict:
        if key == 'freqs':
            G_df.loc[(G_df.simulation == simu_ind) & (G_df.species == 'A'), 'start_freq'] = simu_dict[key][0]
            G_df.loc[(G_df.simulation == simu_ind) & (G_df.species == 'B'), 'start_freq'] = simu_dict[key][1]
            G_df.loc[(G_df.simulation == simu_ind) & (G_df.species == 'C'), 'start_freq'] = simu_dict[key][2]
        else:
            G_df.loc[G_df.simulation == simu_ind, key] = simu_dict[key]

G_df.to_csv(os.path.join(working_dir, "results", "G_df_fixed_lambda_several_iterations.csv"), index=False, header=True)


plt.figure(1)
ax = sns.lineplot(x='simulation', y='start_freq', hue='species', data=G_df[G_df.species != 'coop_advantage'])
ax.grid(b=True)
if SAVE_LAST_FIGURE:
    plt.savefig(os.path.join(working_dir, "results", "species_fraction_during_iterations.png"))
