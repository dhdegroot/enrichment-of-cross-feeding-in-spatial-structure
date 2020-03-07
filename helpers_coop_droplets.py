from itertools import product

import numpy as np
import pandas as pd
import scipy
from matplotlib import pyplot as plt
from scipy.special import gammaln
from scipy.stats import poisson


def findTriplets(lst, key):
    def valid(val):
        return sum(val) == key

    return list(filter(valid, list(product(lst, repeat=3))))


def log_factorial(x):
    """Returns the logarithm of x!
    Also accepts lists and NumPy arrays in place of x."""
    return gammaln(np.array(x) + 1)


def multinomial(xs, ps):
    n = sum(xs)
    xs, ps = np.array(xs), np.array(ps)
    result = log_factorial(n) - sum(log_factorial(xs)) + sum(xs * np.log(ps))
    return np.exp(result)


def plot_pdf(val_freq_array, plot_type='Beautiful', xlabel='fraction cooperators', species_name='A'):
    vals = val_freq_array[:, 0]
    freqs = val_freq_array[:, 1]
    mean = np.sum(vals * freqs) / np.sum(freqs)

    bandwidth = 0.1 * vals.std() * vals.size ** (-1 / 5.)
    support = np.linspace(min(vals) - bandwidth, max(vals) + bandwidth, 200)
    supp_width = support[1] - support[0]

    if plot_type == 'Beautiful':
        kernels = []
        for ind, val_i in enumerate(vals):
            if ind % 1000 == 0:
                print("Calculating mu dist: '{0}' out of '{1}'".format(ind, len(vals)))
            kernel = scipy.stats.norm(val_i, bandwidth).pdf(support)
            kernel = kernel * freqs[ind]
            kernels.append(kernel)

        density = np.sum(kernels, axis=0)
        density /= scipy.integrate.trapz(density, support)
    elif plot_type == "Fast":
        density = np.zeros(len(support))
        for ind, supp_i in enumerate(support):
            supp_i_lb = supp_i - supp_width / 2
            supp_i_ub = supp_i + supp_width / 2
            density[ind] = np.sum(freqs[np.where((vals > supp_i_lb) & (vals <= supp_i_ub))])

    if plot_type != 'None':
        plt.figure()
        plt.plot(support, density)
        plt.xlabel(xlabel)
        plt.ylabel('Probability distribution')
        plt.title('PDF for ' + species_name)

    return mean


def calc_for_one_avg_nt(simu_dict, avg_nt=8, PRINT_ALOT=False):
    freqs = simu_dict['freqs']
    CC0 = simu_dict['CC0']
    CCA_ind = simu_dict['CCA_ind']
    CCB_ind = simu_dict['CCB_ind']
    CCC_ind = simu_dict['CCC_ind']
    adv_cheat = simu_dict['adv_cheat']

    """First, find the range of the number of cells we expect in droplets"""
    nt_range = np.arange(poisson.ppf(0.0001, avg_nt), poisson.ppf(0.9999, avg_nt)).astype(int)
    droplet_df = pd.DataFrame(
        columns=['na', 'nb', 'nc', 'prob_droplet', 'pa', 'pb', 'pc', 'coop_frac', 'mua', 'mub', 'muc'])

    """Then create all possible droplet compositions for one number of cells"""
    for nt in nt_range:
        if nt % 10 == 0:
            if PRINT_ALOT:
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
            coop_frac = 0.0 if denom == 0 else (2 * min(na, nb)) / denom

            # Find growth rates for all species in this droplet
            # Growth factor option 1
            # end_ab = CC0 * coop_frac * (1 / 2) * (1 - (nc / nt)) if nt != 0 else 0.
            # mua = (1 / na) * end_ab if na != 0 else 0.
            # mub = (1 / nb) * end_ab if nb != 0 else 0.
            # muc = (1 / nc) * (CC0 * coop_frac + CCC_ind) * (nc / nt) if nc != 0 else 0.

            # Find growth rates for all species in this droplet
            # Growth factor option 2
            # end_ab = CC0 * coop_frac * (1 / 2) * coop_frac if nt != 0 else 0.
            # mua = (1 / na) * end_ab if na != 0 else 0.
            # mub = (1 / nb) * end_ab if nb != 0 else 0.
            # muc = (1 / nc) * (CC0 * coop_frac + CCC_ind) * (1 - coop_frac) if nc != 0 else 0.

            # Find growth rates for all species in this droplet
            # Growth factor option 3 (best option)
            end_fc = (1 + adv_cheat) * (nc / nt) - adv_cheat * (nc / nt) ** 2 if nt != 0 else 0.
            end_fab = (1 / 2) * (1 - end_fc)

            mua = (1 - coop_frac) * (1 / na) * CCA_ind + coop_frac * (1 / na) * CC0 * end_fab if na != 0 else 0.
            mub = (1 - coop_frac) * (1 / nb) * CCB_ind + coop_frac * (1 / nb) * CC0 * end_fab if nb != 0 else 0.
            muc = (1 - coop_frac) * (1 / nc) * CCC_ind + coop_frac * (1 / nc) * CC0 * end_fc if nc != 0 else 0.

            droplet_df = droplet_df.append(
                {'na': na, 'nb': nb, 'nc': nc, 'prob_droplet': prob_droplet, 'pa': p_species[0], 'pb': p_species[1],
                 'pc': p_species[2], 'coop_frac': coop_frac, 'mua': mua, 'mub': mub, 'muc': muc}, ignore_index=True)

    # Replace all nans by zero. These nans occur when no cells are around, so growth rate should be zero then.
    return droplet_df
