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
N_ITERATIONS = 40
N_CCC_ind = 20
N_CCAB_ind = 20

"""Get parameters"""
# These parameters can be changed in the helpers_coop_droplets.py-file
# For the figure, I have overwritten all these parameters later in the file, except for CC0 and adv_cheat. For those
# I have used:
# CC0 = 750
# adv_cheat = 2/5
freqs, CC0, CCA_ind, CCB_ind, CCC_ind, adv_cheat = get_starting_parameters()

avg_nt = 2  # This nt was picked because it was the lambda that gave a large benefit to coops on average
start_freqs = np.array([0.5, 0.01, 0.49])

start_CCC_ind_factor = 0.05
end_CCC_ind_factor = 0.5
CCC_ind_factors = np.exp(np.linspace(np.log(start_CCC_ind_factor), np.log(end_CCC_ind_factor), N_CCC_ind))

start_CCAB_ind_factor = 0
end_CCAB_ind_factor = 1
CCAB_ind_factors = np.linspace(start_CCAB_ind_factor, end_CCAB_ind_factor, N_CCAB_ind)

simu_CC_inds = product(CCC_ind_factors, CCAB_ind_factors)

"""Calculate for different cheater starter fractions, if the cooperators or the cheaters will win"""
G_df = pd.DataFrame(
    columns=['lambda', 'simulation', 'species', 'G', 'iteration', 'start_freq', 'CC0', 'CCC_ind', 'CCA_ind', 'CCB_ind',
             'adv_cheat'])

for simu_ind, CC_ind_factors in enumerate(simu_CC_inds):
    CCC_ind_factor = CC_ind_factors[0]
    CCAB_ind_factor = CC_ind_factors[1]
    G_df_iter = pd.DataFrame(columns=['lambda', 'simulation', 'species', 'G', 'iteration'])
    print("Starting simulation '{0}'".format(simu_ind))
    freqs_simu = start_freqs.copy()
    CCC_ind = CCC_ind_factor * CC0
    CCA_ind = CCB_ind = CCAB_ind_factor * CCC_ind
    iter_dict_list = []
    for iter_ind in range(N_ITERATIONS):
        if iter_ind % 10 == 9:
            print("Computing '{0}'th iteration".format(iter_ind + 1))
        iter_dict = {'freqs': freqs_simu.copy(), 'CC0': CC0, 'CCC_ind': CCC_ind, 'CCA_ind': CCA_ind, 'CCB_ind': CCB_ind,
                     'adv_cheat': adv_cheat}
        iter_dict_list.append(iter_dict)

        # Calculate growth dataframe for one simulation
        G_df_simu_tidy, end_freqs = calc_for_one_simu(simu_ind, iter_dict, avg_nt, spec_names=spec_names,
                                                      MAKE_PLOTS=MAKE_PLOTS,
                                                      PRINT_ALOT=PRINT_ALOT, SAVE_FIGURES=SAVE_FIGURES,
                                                      MAKE_EACH_SIMU_FIGURE=MAKE_EACH_SIMU_FIGURE)
        G_df_simu_tidy['iteration'] = iter_ind
        G_df_iter = G_df_iter.append(G_df_simu_tidy)

        for spec_ind, spec in enumerate(spec_names_caps):
            growth_factor = G_df_simu_tidy[G_df_simu_tidy.species == spec].G.values[0]
            freqs_simu[spec_ind] *= growth_factor

        freqs_simu = freqs_simu / np.sum(freqs_simu)

    G_df_iter['start_freq'] = 0
    G_df_iter['CC0'] = G_df_iter['CCC_ind'] = G_df_iter['CCA_ind'] = G_df_iter['CCB_ind'] = G_df_iter[
        'adv_cheat'] = 0

    for iter_ind in range(N_ITERATIONS):
        iter_dict = iter_dict_list[iter_ind]
        for key in iter_dict:
            if key == 'freqs':
                G_df_iter.loc[(G_df_iter.iteration == iter_ind) & (G_df_iter.species == 'A'), 'start_freq'] = \
                    iter_dict[key][0]
                G_df_iter.loc[(G_df_iter.iteration == iter_ind) & (G_df_iter.species == 'B'), 'start_freq'] = \
                    iter_dict[key][1]
                G_df_iter.loc[(G_df_iter.iteration == iter_ind) & (G_df_iter.species == 'C'), 'start_freq'] = \
                    iter_dict[key][2]
            else:
                G_df_iter.loc[G_df_iter.iteration == iter_ind, key] = iter_dict[key]

    G_df = G_df.append(G_df_iter)

# Add column that gives ratio between CCC_ind and CCAB_ind
G_df['ratio_CCABind_CCCind'] = G_df['CCA_ind'] / G_df['CCC_ind']
G_df['ratio_CCCind_CC0'] = G_df['CCC_ind'] / G_df['CC0']
G_df['ratio_CCABind_CCCind'] = G_df['ratio_CCABind_CCCind'].astype(float)
G_df['ratio_CCCind_CC0'] = G_df['ratio_CCCind_CC0'].astype(float)

G_df.to_csv(os.path.join(working_dir, "results", "figure4_dataG_df_40iterations_lambda2_varyingCCinds.csv"),
            index=False, header=True)

if N_CCAB_ind < 3:
    plt.figure(1)
    ax = sns.lineplot(x='iteration', y='start_freq', hue='ratio_CCCind_CC0', style='ratio_CCABind_CCCind',
                      data=G_df[G_df.species == 'C'])
    ax.grid(b=True)
    if SAVE_LAST_FIGURE:
        plt.savefig(os.path.join(working_dir, "results", "species_fraction_during_iterations.png"))

plt.figure(2)
ax2 = sns.scatterplot(x='ratio_CCCind_CC0', y='ratio_CCABind_CCCind', hue='start_freq',
                      data=G_df[(G_df.species == 'C') & (G_df.iteration == N_ITERATIONS - 1)], s=400)

L = ax2.legend()
L.get_texts()[0].set_text("Cheater freq")
ax2.set_title("Cheater frequency after '{0}' iterations".format(N_ITERATIONS))
ax2.set_xlabel("Background growth cheater/Max. cooperator growth")
ax2.set_ylabel('Background growth cooperator/Background growth cheater')

plt.savefig(os.path.join(working_dir, "results", "Final cheater fraction after iterations.png"))
