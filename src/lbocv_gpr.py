"""Leave-ball-out cross validation for Chris GPR data on Bonne Pierre glacier."""
import os
import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist

# Define random state (for reproducibility)
random_state = 42
# Number of leave-ball-out simulations
nsim = 20
# Ball size (radius)
ball_size = 500
# Working directory (just need to change this to run localy)
work_dir = "/home/atom/ongoing/collabs/bonnepierre_wpof_Ogier/"
# Open input file, name columns, subset to first 3 columns
df = pd.read_csv(os.path.join(work_dir, 'BonnePierre_ice_thickness_GPR_resampled_50cmtxt_reprojectLambert93.txt'), sep="\t", header=None)
df.columns = ["x", "y", "gpr", "gpr_err", "gpr_err2"]
df = df[["x", "y", "gpr"]]

# Define function for LBOCV
def lbocv(df: pd.DataFrame, ball_size: float, nsim: int, random_state: int | None) -> list[pd.DataFrame]:
    """
    Leave block-out 2D spatial for point data.

    :param df: Input dataframe.
    :param ball_size: Size of spatial ball (radius of disk).
    :param nsim: Number of simulated experiments.
    :param random_state: Random seed or state.

    :return: List of dataframes with blocks left out.
    """

    # Get extent of coordinates
    min_x = np.min(df.x)
    max_x = np.max(df.x)
    min_y = np.min(df.y)
    max_y = np.max(df.y)

    # Define random_state
    rng = np.random.default_rng(random_state)

    # Get random center points for the leave-block-out
    x_c = rng.uniform(min_x, max_x, size=nsim)
    y_c = rng.uniform(min_y, max_y, size=nsim)

    # N x 2 array of inputs points
    input_coords = np.array([df.x.values, df.y.values]).T

    list_df = []
    # Loop over number of simulations (computing all pairwise distances at once would be too RAM-intensive)
    for i in range(nsim):

        # Get center for this simulation
        input_center = np.array([[x_c[i], y_c[i]]])

        # Get pairwise distance matrix, squeeze as it's only 1D with a single center coordinate
        pw_dists = cdist(input_coords, input_center, metric="euclidean").squeeze()

        # Keep only pairwise distances smaller than ball radius
        ind_block = pw_dists < ball_size
        print(f"Simulation {i+1}: removing {np.sum(ind_block)} points out of {len(ind_block)}.")

        # Copy dataframe and replace by NaNs
        df_copy = df.copy()
        df_copy.loc[ind_block, "gpr"] = np.nan
        list_df.append(df_copy)

    return list_df


# Run LBOCV
list_df = lbocv(df, ball_size, nsim, random_state)

# Save output, one file per simulation, same format as input
for i in range(nsim):
    list_df[i].to_csv(os.path.join(work_dir, f"bonnepierre_lbocv_size_{ball_size}_sim_{i+1}.txt"), sep="\t", index=False, header=False)

# Quick visualization
import matplotlib.pyplot as plt
for i in range(nsim):
    plt.figure()
    df.plot(x="x", y="y", kind="scatter", ax=plt.gca(), color="black", marker="x", label="All points")
    df_sim1 = list_df[i].copy()
    df_sim1 = df_sim1[~np.isfinite(df_sim1["gpr"])]
    df_sim1.plot(x="x", y="y", kind="scatter", ax=plt.gca(), color="red", marker="o", label=f"Ball left-out for simulation {i+1}")
    plt.legend()
    plt.savefig(os.path.join(work_dir, "lbocv_"+str(i+1)+".png"))
