from itertools import product

import numpy as np
import scipy
from matplotlib import pyplot as plt
from scipy.special import gammaln


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
