using Pkg
Pkg.activate("/scratch-3/cogier/hydraulic_barriers/")

using CairoMakie
using Statistics

output_dir = "/scratch-3/cogier/data/BonnePierre_output/WWFS_analysis"

# -------- params --------
L   = 100.0                       # transect length (m)
nx  = 201                         # samples
x   = range(0, L; length=nx)
H0  = 50.0                        # baseline thickness (m)
Ai  = 910.0                       # ice density (kg/m^3)
Aw  = 1000.0                      # water density (kg/m^3)
ρi_over_ρw = Ai / Aw              # ~0.910
f_list = [1.0]               # flotation fractions to compare
A_srf = 10.0                      # surface curvature amplitude (m)

# -------- shapes --------
ξ(x) = 2x/L - 1                   # map [0,L] -> [-1,1]
shape_cc(x) = -(ξ(x)^2 - 1)       # concave bowl: center high, edges zero

function shape_zero_mean(xgrid)
    s = [shape_cc(xi) for xi in xgrid]
    μ = mean(s)
    return [(si - μ) for si in s]
end

# shapes
S = shape_zero_mean(x)

# reference convex surface (same in both cases)
zs_conv_ref = H0 .- A_srf .* S      

# beds
zb_flat      = zeros(Float64, nx)
zb_conv_same = -A_srf .* S          # convex bed with same curvature amplitude

# panels
geoms = [
    ("Flat bed + convex surface",                 zb_flat,      zs_conv_ref),
    ("Convex bed + convex surf", zb_conv_same, zs_conv_ref),
]

# head (Shreve potential in water-equivalent height)
h_head(f, zs, zb) = ρi_over_ρw .* f .* (zs .- zb) .+ zb

# -------- plot --------
fig = Figure(resolution=(900, 400))

for (j, (title_txt, zb, zs)) in enumerate(geoms)
    ax = Axis(fig[1, j],
        title = title_txt,
        xlabel = "Distance x (m)",
        ylabel = "Elevation (m)",
        aspect = 1.8
    )

    # surface LAST so it is visible
    lines!(ax, x, zs; color=:gray, linewidth=2,
           label=(j == 1 ? "surface" : ""))

    # compute φ and filled φ (for pockets)
    phi = h_head(1.0, zs, zb)
    phi_filled = copy(phi)
    for i in 2:length(phi)
        phi_filled[i] = max(phi_filled[i-1], phi[i])
    end
    hwp = clamp.(phi_filled .- phi, 0, Inf)

    # overlay heads
    for f in f_list
        h = h_head(f, zs, zb)
        lines!(ax, x, h; label=(j == 1 ? "hydraulic head (f=$(Int(f)))" : ""))
    end


    # water pocket shading
    poly!(ax,
        vcat(x, reverse(x)),
        vcat(zb, reverse(zb .+ hwp));
        color = (:lightblue, 0.5),
        strokewidth = 0
    )
    
    # bed
    lines!(ax, x, zb; color=:black, linewidth=1, label=(j == 1 ? "bed" : ""))

   

    # padding
    yall = vcat(zb, zs, phi, [h_head(f, zs, zb) for f in f_list]...)
    ylims!(ax, minimum(yall) - 2, maximum(yall) + 2)

    if j == 1
        # dummy line for legend entry for pockets
        lines!(ax, [NaN], [NaN]; color=:lightblue, linewidth=6, label="water pocket depth (f=1)")
    end
end

# single legend placed outside on the right of the second subplot
Legend(fig[1, 3], fig.content[1,1];
    framevisible=false,
    tellwidth=false,
    tellheight=false
)

save(joinpath(output_dir, "synthetic_phi_profiles.png"), fig, px_per_unit=3)
println("✅ saved: synthetic_phi_profiles.png")
