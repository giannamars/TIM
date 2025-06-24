include("thermostoichwizard.jl")

using .ThermoStoich
using LinearAlgebra
using Plots, LaTeXStrings
gr()
ENV["GR_DPI"] = 600  # Set export DPI
using DataStructures

els = ["CH2O2", "C2H4O2", "C6H6","C7H8", "C3H6O", "C5H8", "CH4O"]
chemFormBiom = [1, 1.8, 0.2, 0.5, 0, 0, 0]

γ_B = 4.2                     # Biomass degree of reduction
ΔH_T = 469.0                 # Heat released per mol e⁻ transferred [kJ/mol e⁻]

results = Dict{String, NamedTuple{(:y_AE, :gamma, :CR), Tuple{Float64, Float64, Float64}}}()

for el in els
    chemical_indices = extract_composition(el) 
    a = chemical_indices[1]  # C
    b = chemical_indices[2]
    c = chemical_indices[3]
    d = chemical_indices[4]
    e = chemical_indices[5]
    f = chemical_indices[6]
    z = 0

    γ_C = (4a + b - 3c - 2d + 5e - 2f - z) / a  # DoR per C-mol

    out = get_lambda(el, chemFormBiom)
    bio = out[4][10]
    sub = out[4][1]

    y = abs(bio / (sub * a))

    # Catabolic heat release (CR)
    numerator = (1 - y * γ_B / γ_C) * γ_C / 4 * ΔH_T
    denominator = (1 - y)
    CR = numerator / denominator

    results[el] = (y_AE = y, gamma = γ_C, CR = CR)
end

y_vals = collect(0.01:0.001:0.99)

# Use OrderedDict to fix the legend (and plot) order
custom_colors = OrderedDict(
    "CH2O2" => "#b35806",
    "C2H4O2" => "#f1a340",
    "C6H6" => "#d8daeb",
    "C7H8" => "#b2abd2",
    "C3H6O" => "#8073ac",
    "C5H8" => "#542788",
    "CH4O" => "#2d004b"
)

custom_names = OrderedDict(
    "CH2O2" => "Formic acid",
    "C2H4O2" => "Acetate",
    "C6H6" => "Benzene",
    "C7H8" => "Toluene",
    "C3H6O" => "Acetone",
    "C5H8" => "Isoprene",
    "CH4O" => "Methanol"
)

# Calculate the golden ratio
golden_ratio = (1 + sqrt(5)) / 2

# Set the plot size (e.g., width = 600 pixels)
width = 600
height = round(Int, width / golden_ratio) # Calculate height based on golden ratio


plt = plot(
    title = "",
    xlabel = "Calorespirometric ratio (CR) " * L"[\mathrm{kJ}\ \mathrm{mol}^{-1}\ \mathrm{CO}_2]",
    ylabel = "Growth yield " * L"[\mathrm{C\!-\!mol}_{\mathrm{B}}\ \mathrm{C\!-\!mol}^{-1}_{\mathrm{VOC}}]",
    size = (500, 540),          # pixels
    xlims = (0, 1500),
    xformatter = :plain,
    legend = (0.7, 0.5),
    legend_title_font = font(8),
    legendtitle = "",
    foreground_color_legend = nothing,
    grid = true,
    guidefontsize = 10,
    tickfontsize = 8,
    legendfontsize = 8,
    titlefontsize = 10,
    margin = 8Plots.mm
)

annotate!(1200, 0.53, text("Open symbols: TEEM", 8))


# Loop in the exact order defined in OrderedDict
for substrate in keys(custom_colors)
    data = results[substrate]
    γ_C = data.gamma
    CR_vals = [((1 - y * γ_B / γ_C) * γ_C / 4 * ΔH_T) / (1 - y) for y in y_vals]
    c = custom_colors[substrate]
    name = custom_names[substrate]

    # Plot solid line with marker, but no label for legend
    plot!(plt, CR_vals, y_vals,
        label = "",
        color = c,
        linestyle = :solid,
        linewidth = 2,
    )

    # Plot dummy point for legend marker
    scatter!(plt, [NaN], [NaN],
        label = name,
        marker = :circle,
        markercolor = :white,
        markerstrokecolor = c,
        markersize = 7,
    )
end

for (substrate, data) in results
    c = get(custom_colors, substrate, :auto)
    scatter!(plt, [data.CR], [data.y_AE],
        label = "",
        markersize = 5,
        marker = :circle,
        markercolor = :white,
        markerstrokecolor = c,
        markerstrokewidth = 1.5,
        linewidth = 1.5,
    )
end

custom_measurements = Dict(
    "CH2O2" => [(CR=195.93, y=0.13), (CR=200.03, y=0.118), (CR=198.67, y=0.122), (CR=218.93, y=0.057), (CR=176.75+50.0, y=0.183), (CR=205.32, y=0.102), (CR=220.11, y=0.053), (CR=183.86+99.0, y=0.164), (CR=190.93, y=0.145), (CR=207.26, y=0.096), (CR=200.56, y=0.116)],
    "C2H4O2" => [(CR=447.35, y=0.48), (CR=452.70, y=0.41), (CR=447.35, y=0.48), (CR=456.92, y=0.34), (CR=476, y=0.27)],
    "C6H6" => [(CR=724.63, y=0.596)],
    "C7H8" => [(CR=736.29, y=0.548)],
    "C3H6O" => [(CR=692.0, y=0.34), (CR=714.68, y=0.40)],
    "CH4O" => [(CR=974.3, y=0.562), (CR=1132.0, y=0.67), (CR=844.78, y=0.401), (CR=932.14, y=0.520), (CR=743.70, y=0.160), (CR=780.38, y=0.267), (CR=802.82, y=0.32)]
)

for (substrate, measurements) in custom_measurements
    c = get(custom_colors, substrate, :auto)

    for point in measurements
        scatter!(plt, [point.CR], [point.y],
            label = "",  # avoid legend duplication
            markersize = 5,
            marker = :utriangle,
            markercolor = c,
            markerstrokecolor = :black,
            markerstrokewidth = 1.0,
            alpha = 0.8  # optional transparency
        )
    end
end

scatter!(plt, [NaN], [NaN],
    label = "Measurements",
    marker = :utriangle,
    markercolor = :white,
    markerstrokecolor = :black,
    markersize = 6
)


# Horizontal line: C-limited
plot!(plt, [0, 520], [1.0, 1.0], color = "black", linewidth = 1.5, linestyle = :solid, label = "")
annotate!(plt, 250, 1.05, text("Energy limited", :black, 8, :center))

# Horizontal line: Energy-limited
plot!(plt, [560, 2000], [1.0, 1.0], color = "black", linewidth = 1.5, linestyle = :solid, label = "")
annotate!(plt, 1250, 1.05, text("Carbon limited", :black, 8, :center))


display(plt)
savefig(plt, "catabolic_yield.svg")