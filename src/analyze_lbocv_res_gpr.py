"""Leave-ball-out cross validation for Chris GPR data on Bonne Pierre glacier."""
import os
import numpy as np
import pandas as pd
from xdem.spatialstats import nd_binning, plot_2d_binning, plot_1d_binning

# Define random state (for reproducibility)
random_state = 42
# Number of leave-ball-out simulations
nsim = 20
# Ball size (radius)
ball_size = 500
# Outline points
fn_with_outline = "/home/atom/ongoing/collabs/bonnepierre_wpof_Ogier/gpr_plusoutline_df.csv"
# Working directory (just need to change this to run localy)
work_dir = "/home/atom/ongoing/collabs/bonnepierre_wpof_Ogier/"
# Open input file, name columns, subset to first 3 columns
from glob import glob
list_fn = glob(os.path.join(work_dir, "bonnepierre_lbocv_size_500_sim_residuals", '*.txt'))
# Keep only outline points
df_outl = pd.read_csv(fn_with_outline)
df_outl = df_outl[df_outl["h"] == 0]
list_df = []
for fn in list_fn:
    # Read simulation dataframe and add number
    df_sim = pd.read_csv(fn, sep="\t")
    df_sim["nsim"] = int(os.path.basename(fn).split("_")[-2])
    # Rename columns for ease of use
    df_sim.columns = ["x", "y", "gpr", "gpr_lbo", "gpr_model", "err", "nsim"]
    # Compute distance to closest GPR point for this LBOCV simulation
    ind_lbocv = np.isnan(df_sim["gpr_lbo"])
    df_sim_out = df_sim[ind_lbocv]
    df_sim_in = df_sim[~ind_lbocv]
    coords_in = np.stack([np.concatenate((df_sim_in["x"].values, df_outl["x"].values)),
                          np.concatenate((df_sim_in["y"].values, df_outl["y"].values))], axis=1)
    coords_out = np.stack([df_sim_out["x"], df_sim_out["y"]], axis=1)
    from scipy.spatial import KDTree
    kdt = KDTree(coords_in)
    dist, _ = kdt.query(coords_out)
    df_sim.loc[ind_lbocv, "dist_closest"] = dist

    list_df.append(df_sim)

# Group all simulations
df = pd.concat(list_df)
# Keep only residuals for LBOCV
ind_lbocv = np.isnan(df["gpr_lbo"])
df = df[ind_lbocv]

# PLOT: Quick look at error distribution across all simulation
import matplotlib.pyplot as plt
plt.figure()
plt.hist(df["err"], bins=100)
plt.xlabel("Thickness error (cm)")
plt.ylabel("Frequency")
plt.show()

# PLOT: Quick look at error magnitude dependency with distance
df["dist_bins"] = pd.cut(df["dist_closest"], bins=10)
df_grouped_d = df.groupby("dist_bins")["err"].std()
dist_mid_values = np.array([d.mid for d in df_grouped_d.index])
# Linear fits in 1D?
dd = np.linspace(0, 350)
p_d = np.polyfit(dist_mid_values, df_grouped_d.values, deg=1)
yy = np.polyval(p_d, dd)
print({f"1-sigma error = {p_d[0]} * distance + {p_d[1]}"})
plt.figure()
plt.scatter(dist_mid_values, df_grouped_d.values, label="Binned")
plt.plot(dd, yy, linestyle="--", color="black", label="Linear fit")
plt.xlabel("Distance to closest obs (m)")
plt.ylabel("STD of thickness errors (m)")
plt.ylim((0, 35))
plt.legend()
plt.savefig(os.path.join(work_dir, "fig_error_dist_d.png"), dpi=400)

# Thickness
df["h_bins"] = pd.cut(df["gpr"], bins=10)
df_grouped_h = df.groupby("h_bins")["err"].std()
h_mid_values = np.array([h.mid for h in df_grouped_h.index])
# Linear fits in 1D?
hh = np.linspace(0, 100)
p_h = np.polyfit(h_mid_values, df_grouped_h.values, deg=1)
yy = np.polyval(p_h, hh)
plt.figure()
plt.scatter(h_mid_values, df_grouped_h.values, label="Binned")
plt.plot(hh, yy, linestyle="--", color="black", label="Linear fit")
plt.xlabel("Thickness (m)")
plt.ylabel("STD of thickness errors (m)")
plt.ylim((0, 30))
plt.legend()
plt.savefig(os.path.join(work_dir, "fig_error_dist_h.png"), dpi=400)

# DO WE NEED 2D BINNING? LET'S LOOK AT VARIABILITY WITH THICKNESS AFTER STANDARDIZATION BY DISTANCE

# Standardization by variability (once fixed with outline)
x = df["dist_closest"]
y = np.polyval(p_d, x)
df["err"] /= y

# Thickness
df["h_bins"] = pd.cut(df["gpr"], bins=10)
df_grouped_h = df.groupby("h_bins")["err"].std()
h_mid_values = np.array([h.mid for h in df_grouped_h.index])
# Linear fits in 1D?
hh = np.linspace(0, 100)
p_h = np.polyfit(h_mid_values, df_grouped_h.values, deg=1)
yy = np.polyval(p_h, hh)
plt.figure()
plt.scatter(h_mid_values, df_grouped_h.values, label="Binned")
plt.plot(hh, yy, linestyle="--", color="black", label="Linear fit")
plt.xlabel("Thickness (m)")
plt.ylabel("STD of thickness errors (m)")
plt.ylim((0, 2))
plt.legend()
plt.savefig(os.path.join(work_dir, "fig_error_dist_h_std.png"), dpi=400)

# NO: WE DON'T NEED 2D BINNING

# Spatial correlation analysis

# We need to look at dependencies between error, hence at individual simulations again
from skgstat import Variogram
list_df_vgm = []
for nsim in np.arange(1, 21):
    df_sim = df[df["nsim"] == nsim]
    # Get coordinates of LBOCV points
    coords = np.stack([df_sim["x"], df_sim["y"]]).T
    values = df_sim["err"].values
    # Estimate variogram on error values
    bin_func = np.linspace(0, 500, 20)
    vgm = Variogram(coords, values=values, normalize=False, bin_func=bin_func, fit_method=None, estimator="dowd")
    bins, exp = vgm.get_empirical()
    count = vgm.bin_count

    # Write to dataframe
    df_vgm_sim = pd.DataFrame()
    df_vgm_sim = df_vgm_sim.assign(exp=exp, bins=bins, count=count)
    list_df_vgm.append(df_vgm_sim)

df_vgm = pd.concat(list_df_vgm)

# Aggregate variograms by lags
df_vgm["exp_times_count"] = df_vgm["exp"] * df_vgm["count"]
df_grouped = df_vgm.groupby("bins", dropna=False)

# We take the sum of variance and counts for a weighted mean by pairwise observations
df_vgm_all = df_grouped[["exp_times_count", "count"]].sum()

# Remove the last spatial lag bin which is always undersampled
df_vgm_all.drop(df_vgm_all.tail(1).index, inplace=True)

# And the sum for the count of pairwise observations
df_vgm_all["exp"] = df_vgm_all["exp_times_count"] / df_vgm_all["count"]
df_vgm_all["lags"] = df_vgm_all.index.values
df_vgm_all["err_exp"] = np.nan

# PLOT: Aggregate variogram
from xdem.spatialstats import plot_variogram, fit_sum_model_variogram
func, df_params = fit_sum_model_variogram(list_models=["Gau"], empirical_variogram=df_vgm_all)
print(f"Variogram range: {df_params["range"].values[0]} m")
plot_variogram(df_vgm_all, list_fit_fun=[func], xlim=(0, 350))
plt.savefig(os.path.join(work_dir, "fig_variogram_std.png"), dpi=400)
