import pandas as pd
import os

working_dir = os.getcwd()
G_df = pd.read_csv(os.path.join(working_dir, "results", "saved_results", "200_simus_G_df.csv"))

n_simus = max(G_df.simulation.values)

# Initialize empty dataframe in which we will store all simulation factor

# Run over all simulations
for i in range(n_simus):
    G_df
    # Determine maximal growth benefit

    # Determine corresponding avg_nt

    # Determine range of avg_nt for which cooperator benefit is positive

    # Store everything in results_df with simulation parameters
    

pass
