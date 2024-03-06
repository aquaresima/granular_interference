using DrWatson
@quickactivate "granular_media"

using Revise
using speckles
using .HDF5
import YAML
using JLD2
using Plots
using ColorSchemes
using Printf
using LaTeXStrings
using StatsPlots
using StatsBase
using LsqFit
import speckles: shift_mat

##
# Set the location of the folder and import all files.
conf = YAML.load_file(projectdir("conf.yml"))
full_matrices = conf["paths"]["full_matrices"]
analysis_path = conf["paths"]["analysis"]
@assert isdir(analysis_path)
plots_name = conf["paths"]["plots"]

# Read coarse dimensions
conf = YAML.load_file(projectdir("conf.yml"))
coarse = conf["coarse"]
x = coarse["x"]
y = coarse["y"]
times = coarse["times"]

# Set filenames
filename = "coarsen$(x)x$(y)"
coarse_file = joinpath(analysis_path, "$filename.jld2")
temporal_file(root, t) = joinpath(analysis_path, "$(root)_$t.jld2")

# Data
fit_params = conf["fit"]
func = fit_params["func"]
coarse_time = fit_params["coarse_time"]
config = (@strdict func coarse_time) |> dict2ntuple

video_props = load(joinpath(analysis_path, "preprocessing.jld2")) |> dict2ntuple
coarse_matrices = load(coarse_file) |> dict2ntuple
fits_data = load(joinpath(analysis_path, "$(filename)_fits.jld2")) |> dict2ntuple
band_data =
    load(joinpath(analysis_path, savename("$filename", config, "jld2"))) |> dict2ntuple
corr_data = load(joinpath(analysis_path, "$(filename)_correlations.jld2")) |> dict2ntuple


functions = [:stretched, :gaussian, :cosh, :half]
@unpack speeds, matrices = coarse_matrices
@unpack speeds, times = corr_data
@unpack inverse_fits, fits = fits_data
@unpack shearband = band_data
@unpack norm = video_props
keys(fits_data)

# fits_file = joinpath(analysis_path,"coarsen10_fits.jld")
# times_fits =
# load(fits_file,"fits")
# maxima = load(fits_file,"maxima")
# inverse_fits = load(fits_file,"inverse")
# shearband = load(fits_file,"shearband")


##
#=================================================
                    Plots
=================================================#
myplotsdir(args...) = projectdir(plots_name, args...)

gr()
default(
    guidefontsize=18,
    tickfontsize=13,
    titlefontsize=15,
    grid=false,
    colorbar_tickfontsize=13,
    colorbar_titlefontsize=18,
    legendfontsize=12,
)

intensity = ones(256, 1)
for x = 1:256
    intensity[x] = 256 - x
end

metapixel_size = 16.8 / 42
_xticks = (0:80:240, 0:3:9)
xticks = (0:80:240, 0:3:9)
_yticks = (0:105:420, reverse(0:4:16))

## Main text Figure 2: sample of interference pattern
N = 7
if isfile(joinpath(full_matrices, file))
    files = readdir(full_matrices)
    file = files[N]
    _mat = read(h5open(joinpath(full_matrices, file), "r"))["matrix"]
    p = heatmap(
        speckles.shift_mat(_mat, false)[:, :, 1000] .* norm,
        c=:greys,
        xticks=_xticks,
        yticks=_yticks,
        colorbar=false,
        clims=(60, 120),
    )
    plot!(p, xlabel="width (mm)", ylabel="z (mm)", colorbar=false)
    layout = @layout [b a{0.7w} c{0.1w} _]
    plot!(
        twinx(),
        ytick=(0:100:400, reverse(0:100:420)),
        ylims=(-20, 400),
        ylabel="z (px)",
    )
    intensity_y = (range(0, 256, 5), range(60, 120, 5))
    q = heatmap(
        reverse(intensity[1:end, :]),
        ymirror=true,
        yticks=intensity_y,
        c=:greys,
        cbar=false,
        xticks=:none,
        ylabel=L"\langle" * "Pixel intensity" * L"\rangle",
        box=:on,
        grid=:off,
        rotation=-45,
        guidefontsize=13,
        topmargin=10Plots.mm,
        bottommargin=10Plots.mm,
        leftmargin=10Plots.mm,
    )
    pp = plot(
        [plot(frame=:none), p, q]...,
        layout=layout,
        plot_title="speed: ω=$(speeds[N])",
    )
    savefig(pp, myplotsdir("Fig2.pdf"))
    p = heatmap(
        mean(shift_mat(_mat, false), dims=3)[:, :, 1],
        c=:greys,
        xticks=_xticks,
        yticks=_yticks,
        colorbar=false,
        clims=(70, 130),
    )
    layout = @layout [b a{0.7w} c{0.1w} _]
    plot!(p, xlabel="width (mm)", ylabel="z (mm)", colorbar=false)
    plot!(
        twinx(),
        ytick=(0:100:400, reverse(0:100:420)),
        ylims=(-20, 400),
        ylabel="z (px)",
    )
    intensity_y = (range(0, 256, 5), range(70, 130, 5))
    q = heatmap(
        reverse(intensity[1:end, :]),
        ymirror=true,
        yticks=intensity_y,
        c=:greys,
        cbar=false,
        xticks=:none,
        ylabel=L"\langle" * "Pixel intensity" * L"\rangle",
        box=:on,
        grid=:off,
        rotation=-45,
        guidefontsize=13,
        topmargin=10Plots.mm,
        bottommargin=10Plots.mm,
        leftmargin=10Plots.mm,
    )
    pp = plot([plot(frame=:none), p, q]..., layout=layout)
    savefig(pp, myplotsdir("Fig2_not_norm.pdf"))
else
    @warn "Figure 2 skipped. Enable 'run_long' "
end




## Main text Figure 3
# Plot the average intensity of the interference pattern as a function of the rotation speed.
xticks = (0:80:240, 0:3:9)
plots = []
if _run_long
    for (file, speed) in collect(zip(files, speeds))[[1, 3, 6, 9]]
        _mat = read(h5open(joinpath(full_matrices, file), "r"))["matrix"]
        flat = StatsBase.mean(shift_mat(_mat, false), dims=3)[:, :, 1]
        push!(
            plots,
            heatmap(
                norm .* flat,
                title="ω = $speed (s⁻¹)",
                c=:greys,
                titlefontsize=14,
                xticks=_xticks,
                yticks=_yticks,
                clims=(50, 120),
            ),
        )
        if length(plots) < 6
            plot!(colorbar=false)
        end
    end
    plot(plots...)
    layout = @layout [
        grid(2, 2) a{0.07w}
    ]
    p = heatmap(
        reverse(intensity),
        ymirror=true,
        yticks=(range(0, 256, 5), range(50, 120, 5)),
        c=:greys,
        cbar=false,
        xticks=:none,
        ylabel=L"\langle" * "Pixel intensity" * L"\rangle",
        box=:on,
        grid=:off,
    )
    plot!(
        plots[3],
        ylabel="                   z (mm)",
        xlabel="                            width (mm)",
    )
    mean_pixels = plot([plots...; p]..., layout=layout)
    savefig(mean_pixels, myplotsdir("Fig3.pdf"))

    mean_pixels
else
    @warn "Figure 3 skipped. Enable 'run_long' "
end

## Main text Figure 4
_xticks = (1:8:25, 0:3:9)
_yticks = (1:10.2:42, reverse(0:4:16))
N = 7
file = files[N]
# _mat = read(h5open(joinpath(full_matrices,file),"r"))["matrix"]

my_corr = correlations[3, N][:, :, 15]
p = heatmap(
    my_corr,
    # 100 ./ fits_all[2, :, :, 3, 7],
    c=:viridis,
    xticks=_xticks,
    yticks=_yticks,
    colorbar=false,
    clims=(0, 1.0),
)
plot!(p, xlabel="width (mm)", ylabel="z (mm)", colorbar=false)
layout = @layout [b a{0.7w} c{0.1w} _]
plot!(
    twinx(),
    ytick=(0:100:400, reverse(0:10:42)),
    ylims=(-20, 400),
    ylabel="z (metapixel)",
)
intensity_y = (range(0, 256, 5), round.(range(0, 1.0, 5), digits=2))
q = heatmap(
    reverse(intensity[1:end, :]),
    ymirror=true,
    yticks=intensity_y,
    c=:viridis,
    cbar=false,
    xticks=:none,
    # ylabel = "Decorrelation time " * L"(s^{-1})",
    ylabel="Autocorrelation",
    box=:on,
    grid=:off,
    rotation=-0,
    guidefontsize=13,
    topmargin=10Plots.mm,
    bottommargin=10Plots.mm,
    leftmargin=10Plots.mm,
)
# annotate!((2., 1.1), text(L" \times 10^{-2}"))
pp = plot([plot(frame=:none), p, q]..., layout=layout)
# plot!(plot_title="speed: ω=$(speeds[N])")
savefig(pp, myplotsdir("Fig4.pdf"))
pp

## Main text Figure 5
# Plot the correlation function as a function of the rotation speed (A) and depth (B). Also compare correlations of smoothed and non smoothed data (C1, C2).
speeds_z = ones(1, 9)
for x = 1:9
    speeds_z[x] = x
end
cs = cgrad(:roma, 1:9)[1/9:1/9:1]
p = plot()
@unpack correlations = corr_data
for s in eachindex(speeds[1:end])
    s == 1 && continue
    @info "Plotting correlation for speed: $s, coarse_time: 1"
    v = StatsBase.mean(correlations[1, s][25, :, 1:300], dims=1)[:]
    p = plot!(1:300, v, label=false, c=cs[s], lw=4)
end
# plot!(ylims=(0,1.2), size=(800,600))
p = plot(
    p,
    xlabel="Time (s)",
    ylabel="g(z,t)",
    lw=3,
    xticks=(0:90:90*3, 0:6:18),
    guidefontsize=18,
    tickfontsize=13,
    titlefontsize=15,
    grid=false,
    yticks=(0:0.5:1, 0:0.5:1),
    ylims=(-0.1, 1.5),
)
annotate!((-0.15, 0.98), text("A", :center, 20, :black))

metapixel_size
annotate!(
    6 * 15,
    1.45,
    text("depth: $(round(metapixel_size*25, digits=1))mm", :center, 15, :black),
)

speeds_z
heatmap!(
    speeds_z[:, 2:end],
    c=cs[2:end],
    colorbar=false,
    title="ω (s⁻¹)",
    yticks=:none,
    titlefontsize=13,
    inset_subplots=bbox(0.65, 0.75, 0.3, 0.1, :bottom),
    subplot=2,
    axes=false,
    xticks=(1:2:8, speeds[[2:2:9; 9]]),
    xrotation=-30,
)
#

p0 = plot!(size=(800, 600), xlabel="")
zz = ones(1, 42)
for x = 1:42
    zz[x] = 42 - x
end
PALETTE = 42
cs = cgrad(:inferno, 1:PALETTE)[collect(1:PALETTE)./PALETTE]
# speeds
# times
# correlations
v = StatsBase.mean(correlations[1, 4][:, :, 1:300], dims=2)[:, 1, :]
p1 = plot(
    v',
    c=cs',
    # xlims = (, 90 * 3),
    legend=false,
    title="",
    ylims=(-0.1, 1.6),
    yticks=(0:0.5:1, 0:0.5:1),
    xticks=(0:90:90*3, round.(Int, 0:6:6*3)),
    xlabel="Time (s)",
)

plot!(topmargin=5Plots.mm, ylabel="g(z,t)")
annotate!((-0.15, 1), text("B", :center, 20, :black))
annotate!(
    6 * 15,
    1.45,
    text(
        "speed: $(round(parse(Float64,speeds[4]), digits=1)) " * L"s^{-1}",
        :center,
        15,
        :black,
    ),
)
heatmap!(
    zz,
    c=cs,
    colorbar=false,
    yticks=:none,
    xticks=(1:10:42, 0:4:16),
    inset_subplots=bbox(0.65, 0.75, 0.3, 0.1, :bottom),
    subplot=2,
    title="z (mm)",
    titlefontsize=12,
    ylabel="",
)

#
cs = cgrad(:inferno, 1:PALETTE)[collect(1:PALETTE)./PALETTE]
v = StatsBase.mean(correlations[1, 4][:, :, 1:300], dims=2)[:, 1, :]
v2 = StatsBase.mean(correlations[3, 4][:, :, 1:300], dims=2)[:, 1, :]
p3 = plot(
    v',
    c=cs',
    xlims=(-1, 20),
    legend=false,
    title="",
    ylims=(-0.1, 1.5),
    # ylims = (-0.1, 1.6),
    yticks=(0:0.5:1, 0:0.5:1),
    ylabel="g(z,t)",
    xticks=(0:15:15, round.(Int, 0:15/15:15/15)),
    xlabel="Time (s)",
)
annotate!(p3, (-0.35, 1), text("C1", :center, 20, :black))
annotate!(p3, 10, 1.45, text("raw data", :center, 15, :black))
p4 = plot(
    v2',
    c=cs',
    xlims=(-1, 20),
    legend=false,
    title="",
    ylims=(-0.1, 1.5),
    yticks=(0:0.5:1, 0:0.5:1),
    ylabel="g(z,t)",
    xticks=(0:15:90, round.(Int, 0:15/15:90/15)),
    xlabel="Time (s)",
)
annotate!(p4, -5, 1.45, text("C2", :center, 20, :black))
annotate!(p4, 11, 1.45, text("τ = 0.3 s", :center, 15, :black))


p = plot(p3, p4, layout=(1, 2))
decay = plot(p0, p1, p, layout=(3, 1))

plot!(size=(600, 800), leftmargin=5Plots.mm)
savefig(decay, myplotsdir("Fig5.pdf"))
decay

## Main text Figure 6
PALETTE = 9
cs = cgrad(:roma, 1:PALETTE)[collect(1:PALETTE)./PALETTE]
_xticks = (0:10.5:42, 0:4:16)
mat = inverse_fits[2, :, 3, :]
p = plot()
argmax(mat)
xx = axes(mat, 1)
for p = 2:9
    # plot!(reverse(1:42), 100*mat[:,p], c=cs[p], legend=false)
    scatter!(reverse(xx), 100 * mat[:, p], c=cs[p], legend=false, msc=cs[p], ms=6)
    # plot!(ylims=(0,0.20))
end
fit_all = plot!(ylabel="ν" * L"_{C} (s^{-1})", xlabel="z (mm)", title=" ")
ll = maximum(mat * 100) * 1.09
annotate!((-1.0, ll, text(L" \times 10^{-2}")))
plot!(xticks=_xticks, xlabel="z (mm)") ## z (mPixel)
plot!(leftmargin=5Plots.mm)
heatmap!(
    speeds_z[:, 1:8],
    c=cs[2:9],
    colorbar=false,
    title="ω (s⁻¹)",
    yticks=:none,
    titlefontsize=13,
    inset_subplots=bbox(0.65, 0.8, 0.3, 0.07, :bottom),
    subplot=2,
    axes=false,
    xticks=(1:2:8, speeds[[2:2:9;]]),
    xrotation=-30,
)
savefig(p, myplotsdir("Fig6.pdf"))
p

## Main text Figure 7
fits_data.inverse_fits
ss = parse.(Float64, speeds[2:end])
mat = inverse_fits[2, :, 3, 2:end]

xx = axes(mat, 1)
speeds[2:end]
mat
p1 = heatmap(
    speeds[2:end],
    xx,
    100 * mat,
    c=:redsblues,
    xlabel="ω (s⁻¹)",
    ylabel="z (mm)",
    yticks=(0:10:40, reverse(0:4:16)),
    colorbar_title="\nν" * L"_{C} (s^{-1})",
    cbarfontsize=12,
    title=" ",
    margin=5Plots.mm,
    rightmargin=14Plots.mm,
    xrotation=-45,
)
annotate!((9.5, 48, text(L" \times 10^{-2}")))
annotate!(-1.2, 50, text("A", :center, 20, :black))

ss = parse.(Float64, speeds[2:end])
mat = inverse_fits[4, :, 1, 2:end]
p2 = heatmap(
    speeds[2:end],
    xx,
    mat,
    c=:redsblues,
    xlabel="ω (s⁻¹)",
    ylabel="z (mm)",
    yticks=(0:10:40, reverse(0:4:16)),
    colorbar_title="\nν" * L"_{H} (s^{-1})",
    cbarfontsize=12,
    title=" ",
    margin=5Plots.mm,
    rightmargin=14Plots.mm,
    xrotation=-45,
)
annotate!(-1.2, 50, text("B", :center, 20, :black))
# annotate!((10, 45, text(" x 10²")))
p = plot(p1, p2, layout=(2, 1), size=(600, 800))
savefig(p, myplotsdir("Fig7.pdf"))
p


## Main text Figure 8
my_speeds = [parse(Float64, x) for x in speeds]
@unpack func, coarse_time = config
t = findfirst(x -> x == coarse_time, times)
f = findfirst(x -> String(x) == func, functions)
cs = cgrad(:roma, 1:PALETTE)[collect(1:PALETTE)./PALETTE]
taus = [inverse_fits[f, 42-round(Int, shearband.max[x]), t, x+1] for x = 1:8]

mp = metapixel_size / 2
xs_max = [shearband.max[x] * metapixel_size - mp for x = 1:8]
xs_min = [shearband.min[x] * metapixel_size - mp for x = 1:8]
p0 = scatter(
    my_speeds[2:end],
    100 .* taus[1:end],
    label="",
    c=:black,
    ylabel=L"\nu_C^{max} (s^{-1})",
    ms=12,
    xticks=:none,
    topmargin=10Plots.mm,
)
annotate!((-0.28, 1.1), text("A", :center, 20, :black))
annotate!(((-0.1, 1.1), text(L" \times 10^{-2}")))
scatter!(
    [],
    [],
    shape=:circle,
    label="",
    c=:black,
    msc=:black,
    ms=12,
    fg_legend=:transparent,
    legend=:topleft,
    ylims=(0, 3),
)
p1 = scatter(my_speeds[2:end], xs_max, label="", c=:black, ms=12, shape=:diamond)
scatter!(
    p1,
    my_speeds[2:end],
    xs_min,
    label="",
    yerror=shearband.err[2:end] * metapixel_size,
    c=:black,
    ylabel="z (mm) ",
    ms=12,
    shape=:star,
    xlabel="speed " * "ω (s⁻¹)",
)
scatter!(
    [],
    [],
    shape=:diamond,
    label=" " * L"Z_{max}",
    c=:black,
    msc=:black,
    ms=12,
)
scatter!(
    [],
    [],
    shape=:star,
    label=" " * L"Z_{B}",
    c=:black,
    msc=:black,
    ms=12,
    fg_legend=:transparent,
    legend=:topleft,
)
plot!(ylims=(-1, 18), yticks=(0:4:16))
annotate!((-0.28, 1.1), text("B", :center, 20, :black))
p = plot(
    p1,
    layout=(1, 1),
    xrotation=-45,
    size=(00, 200),
    rightmargin=14Plots.mm,
    xlims=(-0.2, 1.7),
    bottommargin=10Plots.mm,
    leftmargin=15Plots.mm,
    legendfontsize=16,
)
plot!(size=(600, 500))
savefig(p, myplotsdir("Fig8.pdf"))
p


## Figure 9
PALETTE = 9
cs = cgrad(:roma, 1:PALETTE)[collect(1:PALETTE)./PALETTE]
mat = inverse_fits[2, :, 5, :]
inverse_fits
gg_max = []
hh_max = []
p = plot()
shift = 0.0
ylims = (0.1, 1.2)

for speed in eachindex(speeds)[2:end]
    hline!([shift * speed], label="", c=:black, alpha=0.4, ls=:dash)
end
plots = []
for speed in eachindex(speeds)[2:end]
    # println(speed)
    hs = inverse_fits[4, :, 1, speed]
    gs = inverse_fits[2, :, 3, speed]
    @show speeds[speed], argmax(hs), argmax(gs)
    push!(hh_max, (xx[end] - argmax(hs)) * metapixel_size + mp)
    push!(gg_max, (xx[end] - argmax(gs)) * metapixel_size + mp)
end
speeds

hh_max

plot()
@. linear(x, p) = p[1] + p[2] * x
fit_hhgg = curve_fit(linear, hh_max, gg_max, [1.0, 1.0])
scatter!(
    hh_max,
    gg_max,
    label="",
    xlabel="H",
    ylabel="G",
    c=cs[2:end],
    msc=:auto,
    ms=10,
)
hh_max
plot!(
    ylabel=L"Z_{max} " * " from " * L"\nu_H" * "(mm)",
    xlabel=L"Z_{max} " * " from " * L"\nu_C" * "(mm)",
    c=:black,
    ls=:dash,
    margin=5Plots.mm,
    rightmargin=14Plots.mm,
    legend=:topleft,
    bottommargin=10Plots.mm,
    topmargin=10Plots.mm,
)
p = plot!(x -> linear(x, fit_hhgg.param), label="", c=:black, ls=:dash)

heatmap!(
    speeds_z[:, 2:end],
    c=cs[2:9],
    colorbar=false,
    title="ω (s⁻¹)",
    yticks=:none,
    titlefontsize=13,
    inset_subplots=bbox(0.6, 0.38, 0.3, 0.07, :bottom),
    subplot=2,
    axes=false,
    xticks=(1:2:8, speeds[[2:2:9; 9]]),
    xrotation=-30,
)
savefig(p, myplotsdir("Fig9.pdf"))
p
