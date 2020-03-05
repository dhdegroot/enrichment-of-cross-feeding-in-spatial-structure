import os

import pandas as pd
from scipy.stats import poisson

from helpers_coop_droplets import *

"""Constants"""
SAVE_FIGURES = True
MAKE_PLOTS = False
working_dir = os.getcwd()
avg_nt = 1
freqs = np.asarray([1 / 4, 1 / 4, 1 / 2])
CC0 = 750
CCC_ind = 100

"""First, find the range of the number of cells we expect in droplets"""
nt_range = np.arange(poisson.ppf(0.0001, avg_nt), poisson.ppf(0.9999, avg_nt)).astype(int)
droplet_df = pd.DataFrame(
    columns=['na', 'nb', 'nc', 'prob_droplet', 'pa', 'pb', 'pc', 'coop_frac', 'mua', 'mub', 'muc'])

"""Then create all possible droplet compositions for one number of cells"""
for nt in nt_range:
    if nt % 10 == 0:
        print("Calculated everything for droplets with less than '{0}' cells".format(nt))
    prob_nt = poisson.pmf(nt, avg_nt)
    droplets = findTriplets(np.arange(0, nt + 1), nt)
    for droplet in droplets:
        na, nb, nc = droplet
        prob_droplet = prob_nt * multinomial([na, nb, nc], freqs)

        # Calculate prob of cell of species to be in this droplet (according to maintext-file)
        p_species = np.asarray([na, nb, nc]) / (freqs * avg_nt) * prob_droplet

        # Find fraction of cooperator pairs
        denom = 2 * min(na, nb) + nc
        if denom == 0:
            coop_frac = 0.0
        else:
            coop_frac = (2 * min(na, nb)) / denom

        # Find growth rates for all species in this droplet
        end_ab = CC0 * coop_frac * (1 / 2) * (1 - (nc / nt)) if nt != 0 else 0.
        mua = (1 / na) * end_ab if na != 0 else 0.
        mub = (1 / nb) * end_ab if nb != 0 else 0.
        muc = (1 / nc) * (CC0 * coop_frac + CCC_ind) * (nc / nt) if nc != 0 else 0.

        droplet_df = droplet_df.append(
            {'na': na, 'nb': nb, 'nc': nc, 'prob_droplet': prob_droplet, 'pa': p_species[0], 'pb': p_species[1],
             'pc': p_species[2], 'coop_frac': coop_frac, 'mua': mua, 'mub': mub, 'muc': muc}, ignore_index=True)

# Plot probability distribution for each species of being in a droplet with pairs of cooperators
# Species A
plot_type = 'Fast' if MAKE_PLOTS else 'None'
val_freq_coop_frac_A = np.transpose(np.vstack((droplet_df.coop_frac.values, droplet_df.pa)))
plot_pdf(val_freq_coop_frac_A, xlabel='fraction cooperators', species_name='A', plot_type=plot_type)
if SAVE_FIGURES:
    plt.savefig(os.path.join(working_dir, "results", "pdf_frac_coop_A.png"))

# Species A growth rate
val_freq_G_A = np.transpose(np.vstack((droplet_df.mua.values, droplet_df.pa)))
avg_G_A = plot_pdf(val_freq_G_A, xlabel='growth factor', species_name='A', plot_type=plot_type)
if SAVE_FIGURES:
    plt.savefig(os.path.join(working_dir, "results", "pdf_G_A.png"))

# Species B
val_freq_coop_frac_B = np.transpose(np.vstack((droplet_df.coop_frac.values, droplet_df.pb)))
plot_pdf(val_freq_coop_frac_B, xlabel='fraction cooperators', species_name='B', plot_type=plot_type)
if SAVE_FIGURES:
    plt.savefig(os.path.join(working_dir, "results", "pdf_frac_coop_B.png"))

# Species B growth rate
val_freq_G_B = np.transpose(np.vstack((droplet_df.mub.values, droplet_df.pb)))
avg_G_B = plot_pdf(val_freq_G_B, xlabel='growth factor', species_name='B', plot_type=plot_type)
if SAVE_FIGURES:
    plt.savefig(os.path.join(working_dir, "results", "pdf_G_B.png"))

# Species C
val_freq_coop_frac_C = np.transpose(np.vstack((droplet_df.coop_frac.values, droplet_df.pc)))
plot_pdf(val_freq_coop_frac_C, xlabel='fraction cooperators', species_name='C', plot_type=plot_type)
if SAVE_FIGURES:
    plt.savefig(os.path.join(working_dir, "results", "pdf_frac_coop_C.png"))

# Species C
val_freq_G_C = np.transpose(np.vstack((droplet_df.muc.values, droplet_df.pc)))
avg_G_C = plot_pdf(val_freq_G_C, xlabel='growth_factors', species_name='C', plot_type=plot_type)
if SAVE_FIGURES:
    plt.savefig(os.path.join(working_dir, "results", "pdf_G_C.png"))

print([avg_G_A, avg_G_B, avg_G_C])
pass
