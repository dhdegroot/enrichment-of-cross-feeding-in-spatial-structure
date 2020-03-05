import numpy as np
from scipy.stats import poisson

from helpers_coop_droplets import *

"""Constants"""
avg_nt = 8
freqs = np.asarray([1 / 4, 1 / 4, 1 / 2])
CC0 = 750
CCC_ind = 100

"""First, find the range of the number of cells we expect in droplets"""
nt_range = np.arange(poisson.ppf(0.0001, avg_nt), poisson.ppf(0.9999, avg_nt)).astype(int)
p_species_sum = np.array([0., 0., 0.])

"""Then create all possible droplet compositions for one number of cells"""
for nt in nt_range:
    prob_nt = poisson.pmf(nt, avg_nt)
    droplets = findTriplets(np.arange(0, nt + 1), nt)
    for droplet in droplets:
        na, nb, nc = droplet
        prob_droplet = prob_nt * multinomial([na, nb, nc], freqs)

        # Calculate prob of cell of species to be in this droplet (according to maintext-file)
        p_species = np.asarray([na, nb, nc]) / (freqs * avg_nt) * prob_droplet
        p_species_sum += p_species

        # Find fraction of cooperator pairs
        denom = 2 * min(na, nb) + nc
        if denom == 0:
            frac_coop = 0.0
        else:
            frac_coop = (2 * min(na, nb)) / denom

print(p_species_sum)