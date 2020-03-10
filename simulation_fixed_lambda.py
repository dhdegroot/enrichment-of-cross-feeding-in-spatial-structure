from helpers_coop_droplets import *

"""Constants"""
working_dir = os.getcwd()
spec_names = ['a', 'b', 'c']
SAVE_LAST_FIGURE = True
SAVE_FIGURES = False
MAKE_EACH_SIMU_FIGURE = False
PRINT_ALOT = False
MAKE_PLOTS = False

"""Get parameters"""
# These parameters can be changed in the helpers_coop_droplets.py-file
freqs, CC0, CCA_ind, CCB_ind, CCC_ind, adv_cheat = get_starting_parameters()

avg_nt = 1.811  # This nt was picked because it was the lambda that gave the largest benefit to coops on average
N_START_FREQCS = 100

"""Calculate for different cheater starter fractions, if the cooperators or the cheaters will win"""
start_freqcs = np.linspace(0.0001, 0.9999, N_START_FREQCS)
simu_dict_list = []

for freqc in start_freqcs:
    freqs_simu = np.array([0.5 - 0.5 * freqc, 0.5 - 0.5 * freqc, freqc])
    simu_dict_list.append(
        {'freqs': freqs_simu, 'CC0': CC0, 'CCC_ind': CCC_ind, 'CCA_ind': CCA_ind, 'CCB_ind': CCB_ind,
         'adv_cheat': adv_cheat})

G_df = pd.DataFrame(columns=['lambda', 'simulation', 'species', 'G'])

# Run all different simulations
for freqc_ind in range(N_START_FREQCS):
    print("Starting simulation '{0}'".format(freqc_ind))
    simu_dict = simu_dict_list[freqc_ind]

    # Calculate growth dataframe for one simulation
    G_df_simu_tidy = calc_for_one_simu(freqc_ind, simu_dict, avg_nt, spec_names=spec_names, MAKE_PLOTS=MAKE_PLOTS,
                                       PRINT_ALOT=PRINT_ALOT, SAVE_FIGURES=SAVE_FIGURES,
                                       MAKE_EACH_SIMU_FIGURE=MAKE_EACH_SIMU_FIGURE)
    G_df = G_df.append(G_df_simu_tidy)

"""At this point, we are left with G_df in which the average growth rate curves for the different simulations are shown,
and with a simu_dict_list where we have stored the parameters for each simulation"""
G_df['freqa'] = G_df['freqb'] = G_df['freqc'] = 0
G_df['CC0'] = G_df['CCC_ind'] = G_df['CCA_ind'] = G_df['CCB_ind'] = G_df['adv_cheat'] = 0
for simu_ind in range(N_START_FREQCS):
    simu_dict = simu_dict_list[simu_ind]
    for key in simu_dict:
        if key == 'freqs':
            G_df.loc[G_df.simulation == simu_ind, 'freqa'] = simu_dict[key][0]
            G_df.loc[G_df.simulation == simu_ind, 'freqb'] = simu_dict[key][1]
            G_df.loc[G_df.simulation == simu_ind, 'freqc'] = simu_dict[key][2]
        else:
            G_df.loc[G_df.simulation == simu_ind, key] = simu_dict[key]

G_df.to_csv(os.path.join(working_dir, "results", "G_df_fixed_lambda_changing_freqc.csv"), index=False, header=True)

plt.figure(1)
ax = sns.lineplot(x='freqc', y='G', data=G_df[G_df.species == 'coop_advantage'])
ax.grid(b=True)
ax.legend(labels=['Growth benefit cooperation'])
if SAVE_LAST_FIGURE:
    plt.savefig(os.path.join(working_dir, "results", "cheater_fraction_vs_coop_benefit.png"))
