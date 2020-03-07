import os

import seaborn as sns

from helpers_coop_droplets import *

"""Constants"""
SAVE_LAST_FIGURE = True
SAVE_FIGURES = False
SAVE_EACH_SIMU_FIGURE = False
MAKE_EACH_SIMU_FIGURE = False
PRINT_ALOT = False
MAKE_PLOTS = True
working_dir = os.getcwd()
N_SIMUS = 20
np.random.seed(3141)
spec_names = ['a', 'b', 'c']

freqs = np.asarray([5 / 12, 5 / 12, 2 / 12])
CC0 = 750
CCA_ind = 10
CCB_ind = 10
CCC_ind = 50
adv_cheat = 2 / 5
first_avg_nt = .15
last_avg_nt = 8
n_avg_nt = 20

# sims_df = pd.DataFrame(columns=['freqs', 'CC0', 'CCC_ind'])
simu_dict_list = []

"""Generate set of parameters for all simulations"""
for ind in range(N_SIMUS):
    freqs_simu = np.random.rand(3)
    freqs_simu /= sum(freqs_simu)
    diff_between_coops = np.random.normal(loc=0.5, scale=0.1, size=1)[0]
    sum_coops = freqs_simu[0] + freqs_simu[1]
    freqs_simu[0] = diff_between_coops * sum_coops
    freqs_simu[1] = (1 - diff_between_coops) * sum_coops

    CC0_simu = CC0
    CCC_ind_simu = CC0_simu * max(np.random.normal(loc=CCC_ind / CC0, scale=0.2 * (CCC_ind / CC0), size=1), 0)[0]
    simu_dict_list.append(
        {'freqs': freqs_simu, 'CC0': CC0_simu, 'CCC_ind': CCC_ind_simu, 'CCA_ind': CCA_ind, 'CCB_ind': CCB_ind,
         'adv_cheat': adv_cheat})

avg_nt_range = np.exp(np.linspace(np.log(first_avg_nt), np.log(last_avg_nt), n_avg_nt))
if avg_nt_range.size > 2:
    MAKE_PLOTS = False
    SAVE_FIGURES = False

# Initialize dataframe in which we are going to store growth factor data
G_df = pd.DataFrame(columns=['lambda', 'simulation', 'species', 'G'])
for simu_ind in range(N_SIMUS):
    print("Starting simulation '{0}'".format(simu_ind))
    simu_dict = simu_dict_list[simu_ind]

    G_array = np.zeros((avg_nt_range.size, 4))
    G_array[:, 0] = avg_nt_range

    for nt_ind, avg_nt in enumerate(avg_nt_range):
        if nt_ind % 10 == 9:
            print("Starting calculation for the '{0}'th lambda".format(nt_ind+1))

        """Calculate for one choice of avg_nt a droplet_dataframe, and the average growth rate"""
        droplet_df_simu = calc_for_one_avg_nt(simu_dict=simu_dict, avg_nt=avg_nt)

        # Plot probability distribution for each species of being in a droplet with pairs of cooperators
        # And the probability distribution for the growth rates, and calculate average growth in the meantime
        plot_type = 'Fast' if MAKE_PLOTS else 'None'
        for spec_ind, spec in enumerate(spec_names):
            val_freq_coop_frac = np.hstack((droplet_df_simu[['coop_frac']], droplet_df_simu[['p' + spec]]))
            plot_pdf(val_freq_coop_frac, xlabel='fraction cooperators', species_name=spec, plot_type=plot_type)
            if SAVE_FIGURES:
                plt.savefig(os.path.join(working_dir, "results", "pdf_frac_coop_" + spec + ".png"))

            val_freq_G = np.hstack((droplet_df_simu[['mu' + spec]], droplet_df_simu[['p' + spec]]))
            G_array[nt_ind, 1 + spec_ind] = plot_pdf(val_freq_G, xlabel='growth factor', species_name=spec,
                                                     plot_type=plot_type)
            if SAVE_FIGURES:
                plt.savefig(os.path.join(working_dir, "results", "pdf_G_" + spec + ".png"))

    # G_diff_avg_coop_min_cheater = np.zeros((G_array.shape[0],1))
    G_diff_avg_coop_min_cheater = 0.5 * (G_array[:, 1] + G_array[:, 2]) - G_array[:, 3]
    # G_array = np.concatenate((G_array, 0.5 * (G_array[:, 1] + G_array[:, 2]) - G_array[:, 3]))
    G_array = np.concatenate((G_array, simu_ind * np.ones((avg_nt_range.size, 1)).astype(int)), axis=1)
    G_df_simu = pd.DataFrame(data=G_array, columns=['lambda', 'A', 'B', 'C', 'simulation'])
    G_df_simu['coop_advantage'] = G_diff_avg_coop_min_cheater
    G_df_simu_tidy = pd.melt(G_df_simu, ['lambda', 'simulation'], var_name='species', value_name='G')
    G_df = G_df.append(G_df_simu_tidy)

    if MAKE_EACH_SIMU_FIGURE:
        ax = sns.lineplot(x="lambda", y="G", hue="species", data=G_df_simu_tidy)
        if SAVE_EACH_SIMU_FIGURE:
            plt.savefig(os.path.join(working_dir, "results", "growth_dep_on_lambda.png"))

        plt.show()

"""At this point, we are left with G_df in which the average growth rate curves for the different simulations are shown,
and with a simu_dict_list where we have stored the parameters for each simulation"""
ax = sns.lineplot(x="lambda", y="G", hue="simulation", data=G_df[G_df.species == 'coop_advantage'])
if SAVE_LAST_FIGURE:
    plt.savefig(os.path.join(working_dir, "results", "growth_dep_on_lambda_different_simulations.png"))

plt.show()
