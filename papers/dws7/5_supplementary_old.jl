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
using NPZ

##
# Set the location of the folder and import all files.
conf = YAML.load_file(projectdir("conf.yml"))
full_matrices = conf["paths"]["full_matrices"]
analysis_path = conf["paths"]["analysis"]
@assert isdir(analysis_path)

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
keys(fits_data)
speeds

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

gr()
default(
    guidefontsize = 18,
    tickfontsize = 13,
    titlefontsize = 15,
    grid = false,
    colorbar_tickfontsize = 13,
    colorbar_titlefontsize = 18,
    legendfontsize = 12,
)
xticks = (0:80:240, 0:3:9)
yticks = (0:105:420, reverse(0:4:16))
intensity = ones(256, 1)
for x = 1:256
    intensity[x] = 256 - x
end
import speckles: shift_mat

metapixel_size = 16 / 42

##
files = readdir(full_matrices)
file = string(joinpath(analysis_path, filename), ".jld2")
speeds
9_000 / 60 / 15
matrices = []
for x = 1:9
    m = read(h5open(full_matrices * "/" * files[x], "r")["matrix"])
    @info "Original size: $(size(m)), in minutes: $(size(m,3)/15/60)"
    _mat = load(file)["matrices"][x]
    @info "Processed size: $(size(_mat)), in minutes: $(size(_mat,3)/15/60)"
    s = size(_mat, 3) / 15 / 60
    push!(matrices, _mat)
end
##

##
plot()
for x = 1:9
    plot!(mean(matrices[x] .+ 20 * x, dims = (1, 2))[1, 1, :])
end
plot!(size = (400, 1500), margin = 10Plots.mm)


## Supplementary Figures
## Normalization matrix, computed by averaging the pixel intensity per rows and columns. 
## This method has been choosen because it leverages the symmetries of the laser beam.
@unpack norm = video_props
norm = norm ./ maximum(norm)
# heatmap(norm)
nn = heatmap(
    norm,
    c = :grays,
    cbar = false,
    title = L"\langle" * "Average pixel intensity" * L"\rangle",
    box = :on,
    grid = :off,
    rotation = -45,
    guidefontsize = 13,
    margin = 10Plots.mm,
)
savefig(nn, plotsdir("SIFig1.pdf"))
nn
##

## Data in Numpy format
# NPZ.npzwrite(datadir("inverse_fits.npz"), Dict("inverse_fits"=>inverse_fits))
# print(size(inverse_fits))
# print(speeds)

@unpack variance = video_props
default(
    guidefontsize = 18,
    legendfontsize = 15,
    fg_color_legend = :transparent,
    tickfontsize = 13,
)
p = plot(
    variance .^ 2,
    xlabel = "Metapixel size (pixels)",
    ylabel = "< σ ",
    lw = 4,
    c = :black,
    label = "Still wheel video",
)
plot!(x -> 1 / x * (maximum(variance)^2), ls = :dash, lw = 4, label = "1/x")
# plot!(maximum(variance)*exp.(-collect(0:25)./10), xlabel="Metapixel size (pixels)", ylabel="Total variance", ls=:dash, lw=4, label="Exp. decay with\ncharacteristic lenght 10 MP")
plot!(lw = 4)# yscale=:log,xscale=:log, legend=false)
savefig(p, plotsdir("SIFig2.pdf"))
p


#
v0 = read(h5open(string(full_matrices, "/V_0.00.h5")))["matrix"]

##
x  =v0
p1 = v0_mean = mean(x, dims=3)[:,:,1] |> heatmap
p2 = histogram(x[1:10:end])
pa = plot(p1,p2, size=(800,400), title="Base")

x = speckles.shift_mat(v0, true)
p1 = v0_mean = mean(x, dims=3)[:,:,1] |> heatmap
p2 = histogram(x[1:10:end])
pb = plot(p1,p2, size=(800,400), title="Abs and add")

x = speckles.shift_mat(v0, false)
p1 = mean(x, dims=3)[:,:,1] |> heatmap
p2 = histogram(x[1:10:end])
pc = plot(p1,p2, size=(800,400), title="Abs only")

a = Int16.(v0)
a[a.<0] .= a[a.<0] .+254
p1 = mean(a, dims=3)[:,:,1] |> heatmap
plot(pa,pb,pc, layout=(3,1), legend=false, size=(800,1200))
##

scatter(hs[2:end],gs[2:end], label="", xlabel="H", ylabel="G")
#

heatmap!(z, c=cs[1:end], colorbar=false, title="ω (s⁻¹)", yticks=:none, titlefontsize=13,
    inset_subplots = bbox(0.65, 0.8, 0.3, 0.07, :bottom), subplot=2, axes=false, xticks=(1:2:10, speeds[[2:2:9;9]]))
savefig(p, joinpath(results_path,"fig7_maxtau.pdf"))

cs = cgrad(:roma, 1:PALETTE)[collect(1:PALETTE) ./PALETTE]
taus_x = [42 - argmax(inverse_fits[:,2,x,5]) for x in 1:8]
taus = 100 .*[maximum(inverse_fits[:,2,x,5]) for x in 1:8]
p = scatter(taus_x, taus, label="", c=cs, msc=cs, ylabel="maximum ν (s⁻¹)", ms=12, xlabel="z (mm)")
plot!(xticks=(0:10:20, 0:4:8), legend=:topleft, xlims=(0:8))
heatmap!(z, c=cs[1:end], colorbar=false, title="ω (s⁻¹)", yticks=:none, titlefontsize=13,
    inset_subplots = bbox(0.65, 0.8, 0.3, 0.07, :bottom), subplot=2, axes=false, xticks=(1:2:10, speeds[[2:2:9;9]]))
savefig(p, joinpath(results_path,"fig7_maxtauinv.pdf"))
plot!(title=" ")
annotate!((0.65, 6.45, text(" x 10²")))

# scatter(1 ./inverse_fits[:,2,:,1])


#
p = groupedbar(xticks=(1:8,string.(speeds)),[-shearband.max -shearband.min], bar_position=:stack,c=[:white cs[6] ], lc=:white, legend=false, ylabel= "z (mm)", xlabel="ω (s⁻¹)", guidefontsize=18, tickfontsize=13, grid=false, ylims=(-42,0),yticks=(reverse(-40:10:0), 0:4:16))
savefig(p, joinpath(results_path,"fig8_bars.pdf"))
##
gr()
layout = @layout [
			a{0.9w} _;_ c{0.0001h}
			]
p = plot(parse.(Float64,speeds[2:end]), shearband.max, seriestypes=[:line, :scatter], c=:red, msc=:red , ms=8, labels=["" "z maximum"],  ylabel="z maximum (mm)", legend=:topleft)
p = plot!(twinx(), parse.(Float64,speeds[2:end]), shearband.min, ylabel="z shear (mm)", label=["" "z shear"], seriestypes=[:line, :scatter], c=:black, msc=:black, ms=8,  legend=:bottomright,xlabel="ω (s⁻¹)")
p =plot(p,plot(frame=:none),layout=layout)
savefig(p, joinpath(results_path,"fig8.pdf"))

##
p = scatter(parse.(Float64,speeds[2:end]), maxima[2,:,5], label="Gaussian Exp. time avg", xlabel="ω (s⁻¹)", ylabel="z (mm)", c=:black, ms=8, shape=:square)
p = scatter!(parse.(Float64,speeds[2:end]), maxima[4,:,1], label="Half width Raw data", xlabel="ω (s⁻¹)", c=:red, msc=:red, ms=8)
plot!(yticks=(30:5:40, reverse(0:2:4)), legendfontsize=12, title=" ", )


savefig(p, joinpath(results_path,"fig9.pdf"))
p
##
p = scatter(maxima[2,:,5], maxima[4,:,1], ylabel="Half width Raw data", xlabel="Gaussian Exp. time avg", c=:black, msc=:black, ms=8, label="")

savefig(p, joinpath(results_path,"fig9_regression.pdf"))
