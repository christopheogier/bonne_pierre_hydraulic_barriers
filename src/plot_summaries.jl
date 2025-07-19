### plot_summaries.jl
using Pkg
Pkg.activate("/scratch-1/cogier/hydraulic_barriers/")
using CSV, DataFrames
include("plot.jl")  # or "plot_summary_overview.jl" if separated

output_dir = "/scratch-3/cogier/data/BonnePierre_output/WWFS_analysis"

df = CSV.read("/scratch-3/cogier/data/BonnePierre_output/WWFS_analysis/WWFS_lake_summary.csv", DataFrame)

#plot_selected_scenarios(df, output_dir)


# plots

using CairoMakie

# JoG setup
two_column_cm   = 17.8
one_column_cm   = 8.6
font_size_pt    = 10
label_font_size = 12

pt_in_cm = 28.3465;

two_column_pt = two_column_cm * pt_in_cm;
one_column_pt = one_column_cm * pt_in_cm;

dpi        = 600
cm_in_inch = 2.54

px_per_unit = two_column_cm / cm_in_inch * dpi / two_column_pt
pt_per_unit = 1

# custom Makie theme for publications
makie_theme = merge(theme_latexfonts(),
                    Theme(; fontsize=font_size_pt,
                          Axis=(spinewidth=0.5,
                                xtickwidth=0.5,
                                ytickwidth=0.5,
                                xticksize=3,
                                yticksize=3),
                          Colorbar=(spinewidth=0.5, tickwidth=0.5, ticksize=3, size=7),
                          Label=(fontsize=label_font_size, font=:bold),
                          Legend=(rowgap=-8, labelsize=8, framewidth=0.25, padding=(2, 2, 2, 2), margin=(4, 4, 4, 4))))

with_theme(makie_theme) do
    # create a fake glacier thickness and velocity field
    x = LinRange(-7, 7, 501)
    y = LinRange(-10, 10, 501)

    # fake glacier thickness
    Hc = @. 200 * exp(-(x / 04)^2 - ((y' + 1) / 10)^2) + # big blob
            100 * exp(-(x / 20)^2 - ((y' - 4) / 04)^2) + # small blob
            10 * sin(x * 3π / 7) * cos(y' * 3π / 10) -   # wiggle
            100                                          # cutoff

    ice_mask = Hc .< 0 # mask for ice-free areas
    Hc[ice_mask] .= 0  # set ice-free areas to zero for countour

    H = copy(Hc) # ice thickness for visualization
    H[ice_mask] .= NaN

    V = @. 1e-12 * H^6 # fake ice velocity

    # make a figure
    fig = Figure(; size=(two_column_pt, 270))

    # two subplots
    axs = (Axis(fig[1, 1][1, 1]; aspect=DataAspect()),
           Axis(fig[1, 2][1, 1]; aspect=DataAspect()))

    # hide the y-axis of the second subplot
    hideydecorations!.((axs[2],))

    # make grid lines visible and set limits
    for ax in axs
        ax.xgridvisible = true
        ax.ygridvisible = true
        limits!(ax, -7, 7, -10, 10)
    end

    # set labels and titles
    axs[1].ylabel = L"y~\mathrm{(km)}"

    axs[1].xlabel = L"x~\mathrm{(km)}"
    axs[2].xlabel = L"x~\mathrm{(km)}"

    axs[1].title = L"H~\mathrm{(m)}"
    axs[2].title = L"V~\mathrm{(m\,a^{-1})}"

    # plot the glacier thickness and velocity
    hms = (heatmap!(axs[1], x, y, H),
           heatmap!(axs[2], x, y, V))

    # add contour lines for the glacier outline
    contour!(axs[1], x, y, Hc; levels=0.1:0.1, linewidth=0.5, color=:black)
    contour!(axs[2], x, y, Hc; levels=0.1:0.1, linewidth=0.5, color=:black)

    # enable interpolation for smoother picture
    foreach(hms) do h
        h.interpolate = true
        h.rasterize   = px_per_unit
    end

    # plot velocity magnitude in logarithmic scale
    hms[2].colorscale = log10

    hms[1].colormap = Reverse(:ice)
    hms[2].colormap = Reverse(:roma)

    hms[1].colorrange = (0, 200)
    hms[2].colorrange = (1e-2, 100)

    # add colorbars
    cbs = (Colorbar(fig[1, 1][1, 2], hms[1]),
           Colorbar(fig[1, 2][1, 2], hms[2]))

    # add labels for panels
    for (label, idx) in zip('a':'b', ((1, 1), (1, 2)))
        Label(fig[idx..., TopLeft()], string(label); padding=(0, 10, 0, 0))
    end

    # adjust layout
    for l in fig.layout.content[1:2]
        colgap!(l.content, 1, Fixed(10))
    end

    # save both raster and vector versions
    save("myplot.pdf", fig; pt_per_unit, px_per_unit)
    save("myplot.png", fig; pt_per_unit, px_per_unit)

    fig
end

