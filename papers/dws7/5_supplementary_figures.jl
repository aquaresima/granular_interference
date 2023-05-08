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
full_matrices  = conf["paths"]["full_matrices"]
analysis_path  = conf["paths"]["analysis"]
@assert isdir(analysis_path)

# Read coarse dimensions
conf = YAML.load_file(projectdir("conf.yml"))
coarse = conf["coarse"]
x = coarse["x"]
y = coarse["y"]
times = coarse["times"]

# Set filenames
filename = "coarsen$(x)x$(y)"
coarse_file = joinpath(analysis_path,"$filename.jld2") 
temporal_file(root,t) = joinpath(analysis_path,"$(root)_$t.jld2")

# Data
fit_params = conf["fit"]
func = fit_params["func"]
coarse_time = fit_params["coarse_time"]
config = (@strdict func coarse_time) |> dict2ntuple

video_props = load(joinpath(analysis_path,"videos_properties.jld2")) |> dict2ntuple
coarse_matrices = load(coarse_file) |> dict2ntuple
fits_data = load(joinpath(analysis_path,"$(filename)_fits.jld2")) |> dict2ntuple
band_data = load(joinpath(analysis_path,savename("$filename", config, "jld2"))) |> dict2ntuple
corr_data = load(joinpath(analysis_path,"$(filename)_correlations.jld2")) |> dict2ntuple


functions = [:stretched, :gaussian, :cosh, :half] 
@unpack speeds, matrices = coarse_matrices
@unpack speeds, times = corr_data
@unpack inverse_fits, fits= fits_data
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
default(guidefontsize=18, tickfontsize=13, titlefontsize=15, grid=false, colorbar_tickfontsize=13, colorbar_titlefontsize=18, legendfontsize=12)
xticks= (0:80:240, 0:3:9)
yticks = (0:105:420, reverse(0:4:16))
intensity = ones(256,1)
for x in 1:256
    intensity[x] = 256-x
end
import speckles: shift_mat

metapixel_size = 16/42

##
files = readdir(full_matrices)
file = string(joinpath(analysis_path,filename),".jld2")
speeds
9_000/60/15
matrices =[]
for x in 1:9
    m = read(h5open(full_matrices*"/"*files[x], "r")["matrix"])
    @info "Original size: $(size(m)), in minutes: $(size(m,3)/15/60)"
    _mat = load(file)["matrices"][x]
    @info "Processed size: $(size(_mat)), in minutes: $(size(_mat,3)/15/60)"
    s = size(_mat,3)/15/60
    push!(matrices, _mat)
end
##


## Supplementary Figures
## Normalization matrix, computed by averaging the pixel intensity per rows and columns. 
## This method has been choosen because it leverages the symmetries of the laser beam.
@unpack norm = video_props
norm = norm ./maximum(norm)
# heatmap(norm)
nn = heatmap(norm , c=:grays, cbar=false, title=L"\langle"*"Average pixel intensity"*L"\rangle", box=:on,grid=:off, rotation=-45, guidefontsize=13, margin=10Plots.mm)
savefig(nn, plotsdir("SIFig1.pdf"))
nn
##

## Data in Numpy format
# NPZ.npzwrite(datadir("inverse_fits.npz"), Dict("inverse_fits"=>inverse_fits))
# print(size(inverse_fits))
# print(speeds)

@unpack variance = video_props
default(guidefontsize=18, legendfontsize=15, fg_color_legend=:transparent, tickfontsize=13)
p =plot(variance, xlabel="Metapixel size (pixels)", ylabel="< σ ", lw = 4, c=:black, label = "Still wheel video")
plot!(maximum(variance)*exp.(-collect(0:25)./10), xlabel="Metapixel size (pixels)", ylabel="Total variance", ls=:dash, lw=4, label="Exp. decay with τ=10")
plot!(lw = 4)
savefig(p, plotsdir("SIFig2.pdf"))
p

## Figure SI3
PALETTE = 9
cs = cgrad(:roma, 1:PALETTE)[collect(1:PALETTE) ./PALETTE]

plot()
for s in 2:9
    plot!(reverse(inverse_fits[2,:,3,s]), c=cs[s], label="")
    scatter!([shearband.max[s]],[inverse_fits[2,43-round(Int,shearband.max[s]),3,s]], c=cs[s], shape=:diamond, ms=12, msc=:auto, label="")
    scatter!([shearband.min[s]],[inverse_fits[2,43-round(Int,shearband.min[s]),3,s]], c=cs[s], shape=:star, ms=12, msc=:auto, label="")
end 
scatter!([[],[]], label=["max" "min"], c=:black, ms=3, shape=[:diamond :star])
p = plot!(legend=:topright, xlabel="Metapixel", ylabel=L"\nu_c (s^{-1})")
heatmap!(speeds_z, c=cs[1:9], colorbar=false, title="ω (s⁻¹)", yticks=:none, titlefontsize=13,
    inset_subplots = bbox(0.65, 0.5, 0.3, 0.07, :bottom), subplot=2, axes=false, xticks=(1:2:10, speeds[[2:2:9;9]]), xrotation=-30,
)
savefig(p, plotsdir("SIFig3.pdf"))
p

## Figure SI4
PALETTE = 9
cs = cgrad(:roma, 1:PALETTE)[collect(1:PALETTE) ./PALETTE]
mat = inverse_fits[2,:, 5, :]
inverse_fits
gg_max = []
hh_max = []
p = plot()
shift=0.0
ylims=(0.1, 1.2)

for speed in eachindex(speeds)[2:end]
    hline!([shift*speed], label="", c=:black, alpha=0.4, ls=:dash)
end
plots=[]
for speed in eachindex(speeds)[2:end]
    # println(speed)
    hs =  inverse_fits[4,:,1,speed]
    gs =  inverse_fits[2,:,3,speed]
    push!(hh_max,(43- argmax(hs))*metapixel_size)
    push!(gg_max,(43- argmax(gs))*metapixel_size)
    p = scatter(hs,gs.+ shift*speed, label="", xlabel="Half width "*L"\nu_H (s^{-1})", ylabel="Gaussian exp "*L"\nu_H (s^{-1})", c =cs[speed],  msc=:auto, ms=8 )
    scatter!([maximum(hs)],[maximum(gs)].+shift*speed, label="",shape=:star, c =:darkred, ms=18, msc=:auto, leftmargin=10Plots.mm)
    plot!(ylims = (0, 0.038))
    annotate!(minimum(hs),0.03, text("$(speeds[speed]) ω/s", :left, 20, :black))
    # plot!(twinx(p), ylims=ylims, yticks=(0.2:0.1:0.9,speeds[2:end]), ylabel="ω (s⁻¹)", yrotation=-45, )
    # plot!(twinx(p), ylims=ylims, yticks=(0.2:0.1:0.9,speeds[2:end]), ylabel="ω (s⁻¹)", yrotation=-45, )
    push!(plots,p)

end

p = plot(plots...,layout=(4,2), legend=false, size=(1200,1300), left_margin=10Plots.mm, bottom_margin=10Plots.mm,)
savefig(p, plotsdir("SIFig4.pdf"))
p
##
p= plot()
q= plot()
for x in 1:8
    hs =  reverse(fits[4,:,1,x+1])
    plot!(p, hs, c=cs[x+1], label="")
    hs =  reverse(inverse_fits[4,:,1,x+1])
    plot!(q, hs, c=cs[x+1], label="")
end
p = plot!(p,ylabel="Half width "*L"\tau_H (s)", xlabel="Metapixel", legend=:topright, left_margin=10Plots.mm, bottom_margin=10Plots.mm, )
heatmap!(speeds_z, c=cs[1:9], colorbar=false, title="ω (s⁻¹)", yticks=:none, titlefontsize=13,
    inset_subplots = bbox(0.2, 0.8, 0.3, 0.07, :bottom), subplot=2, axes=false, xticks=(1:2:10, speeds[[2:2:9;9]]), xrotation=-30,
)
q = plot!(q, xlabel="Metapixel", legend=:topright, left_margin=10Plots.mm, bottom_margin=10Plots.mm, )
plot!(twinx(),ylabel="Half width "*L"\nu_H (1/s)",rightmargin=10Plots.mm)
p = plot(p,q, size=(1100,400))
savefig(p, plotsdir("SIFig5.pdf"))
p

##

## Figure SI6
## Verify that the image has an exponential distribution in light-intensity to verify the measure regard actual spekles. : Goodman Statistical Optics

speeds
p = [histogram(m[:], alpha=1., lw=4,c=:black, bins=-0:5:254, title="ω = $s (s⁻¹)") 
    for (s,m) in collect(zip(speeds,matrices))[[2,4,6,8]]]
plot!(p[3], xlabel="                                               Pixel intensity (255 bits)")
plot(p..., layout=(2,2), legend=false, size=(800,600))
intensity = plot!( legend=false)
savefig(intensity,plotsdir("SIFig6.pdf"))
intensity

