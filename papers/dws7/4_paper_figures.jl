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

##
# Set the location of the folder and import all files.
conf = YAML.load_file(projectdir("conf.yml"))
full_matrices  = conf["paths"]["full_matrices"]
analysis_path  = conf["paths"]["analysis"]
# @assert isdir(full_matrices) 
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

# video_props = load(joinpath(full_matrices,"videos_properties.jld2")) |> dict2ntuple
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

## Supplementary Figure 2
## Normalization matrix, computed by averaging the pixel intensity per rows and columns. 
## This method has been choosen because it leverages the symmetries of the laser beam.
@unpack norm = video_props
norm = norm ./maximum(norm)
# heatmap(norm)
nn = heatmap(norm , c=:grays, cbar=false, title=L"\langle"*"Average pixel intensity"*L"\rangle", box=:on,grid=:off, rotation=-45, guidefontsize=13, margin=10Plots.mm)
savefig(nn, plotsdir("SIFig1.pdf"))
nn

## Main text Figure 1
## Plot 1 sample of interference pattern
files = readdir(full_matrices)
file = files[2]
_mat = read(h5open(joinpath(full_matrices,file),"r"))["matrix"]
p = heatmap(mean(shift_mat(_mat), dims=3)[:,:,1], c=:roma, xticks=xticks, yticks=yticks, clims=(80,256), colorbar=false)
layout = @layout [b a{0.7w} c{0.1w} _]
plot!(p, xlabel= "width (mm)", ylabel="z (mm)", colorbar=false)
plot!(twinx(), ytick=(0:0.25:1, reverse(0:100:400)), ylims=(0,1), ylabel="z (px)")
intensity_y = (range(0,256,5), range(80,256,5))
q = heatmap(reverse(intensity[1:end,:]) , ymirror=true, yticks=intensity_y, c=:roma, cbar=false, xticks=:none, ylabel=L"\langle"*"Pixel intensity"*L"\rangle", box=:on,grid=:off, rotation=-45, guidefontsize=13, topmargin=10Plots.mm, bottommargin=10Plots.mm, leftmargin=10Plots.mm)
pp = plot([plot(frame=:none), p, q ]..., layout=layout)
savefig(pp, plotsdir("Fig2.pdf"))
q

## Main text Figure 3
## Plot the average intensity of the interference pattern as a function of the rotation speed.
plots = []
for (file,speed) in collect(zip(files, speeds))[[1,3,6,9]]
	_mat = read(h5open(joinpath(full_matrices,file),"r"))["matrix"]
    flat =StatsBase.mean(shift_mat(_mat), dims=3)[:,:,1]
    push!(plots, heatmap(norm .* flat, title="ω = $speed (s⁻¹)", c=:roma,titlefontsize=14,xticks=xticks,yticks=yticks, clims=(80,256)))
	if length(plots) <6
		plot!(colorbar=false)
	end
end
plot(plots...)
layout = @layout [
 		grid(2,2) a{0.07w}
 ]
p = heatmap(reverse(intensity), ymirror=true, yticks=(range(0,256,5), range(80,240,5)), c=:roma, cbar=false, xticks=:none, ylabel=L"\langle"*"Pixel intensity"*L"\rangle", box=:on,grid=:off)
plot!(plots[3], ylabel="                   z (mm)", xlabel="                            width (mm)")
mean_pixels = plot([plots...;p]..., layout=layout)
savefig(mean_pixels,plotsdir("Fig3.pdf"))


## Main text Figure 4
## Plot the correlation function as a function of the rotation speed (A) and depth (B). Also compare correlations of smoothed and non smoothed data (C1, C2).
speeds_z = ones(1,9)
for x in 1:9
    speeds_z[x] = x
end
cs = cgrad(:roma, 1:9)[1/9:1/9:1]
p = plot()
@unpack correlations = corr_data
for s in eachindex(speeds[2:end])
    s == 1 && continue
    @info "Plotting correlation for speed: $s, coarse_time: 1"
	v = StatsBase.mean(correlations[1,s][25,:,1:300], dims=1)[:]
	p=plot!(1:300,v, label=false,  c=cs[s], lw=4)
end
# plot!(ylims=(0,1.2), size=(800,600))
p = plot(p, xlabel="Time (s)", ylabel="g(z,t)", lw=3, xticks=(0:90:90*3, 0:6:18), guidefontsize=18, tickfontsize=13,titlefontsize=15, grid=false, yticks=(0:0.5:1, 0:0.5:1), ylims=(0,1.5))
annotate!(-50, 1.45, text("A", :center, 23, :black))

annotate!(6*15, 1.45, text("depth: $(round(metapixel_size*25, digits=1))mm", :center, 15, :black))
heatmap!(speeds_z, c=cs[1:end], colorbar=false, title="ω (s⁻¹)", yticks=:none, titlefontsize=13,
    inset_subplots = bbox(0.65, 0.75, 0.3, 0.1, :bottom), subplot=2, axes=false, xticks=(1:2:9, speeds[[2:2:9;9]]), xrotation=-30)
#

p0 = plot!(size=(800,600), xlabel="",)
zz = ones(1,42)
for x in 1:42
    zz[x] = 42-x
end
PALETTE = 42
cs = cgrad(:inferno, 1:PALETTE)[collect(1:PALETTE) ./PALETTE]
speeds
times
correlations
v = StatsBase.mean(correlations[1,4][:,:,1:300], dims=2)[:,1,:]

p1 =plot(v', c=cs', xlims=(1,), legend=false, title="", ylims=(0,1.6), yticks=(0:0.5:1, 0:0.5:1), xticks=(0:90:90*3, round.(Int,0:6:18) ), xlabel="Time (s)")
plot!(topmargin=5Plots.mm, ylabel="g(z,t)" )
annotate!(-50, 1.5, text("B", :center, 23, :black))
annotate!(6*15, 1.45, text("speed: $(round(parse(Float64,speeds[4]), digits=1)) "*L"s^{-1}", :center, 15, :black))
heatmap!(zz, c=cs, colorbar=false,  yticks=:none,xticks=(1:10:42,0:4:16),
    inset_subplots = bbox(0.65, 0.75, 0.3, 0.1, :bottom), subplot=2, title= "z (mm)", titlefontsize=12, ylabel="")

#
cs = cgrad(:inferno, 1:PALETTE)[collect(1:PALETTE) ./PALETTE]
v = StatsBase.mean(correlations[1,4][:,:,1:300], dims=2)[:,1,:]
v2 = StatsBase.mean(correlations[3,4][:,:,1:300], dims=2)[:,1,:]
p3 = plot(v', c=cs', xlims=(1,50), legend=false, title="", ylims=(0,1.5), yticks=(0:0.5:1, 0:0.5:1), ylabel="g(z,t)",  xticks=(0:45:135, round.(Int,0:45/15:135/15) ), xlabel="Time (s)" )
annotate!(p3,-15, 1.45, text("C1", :center, 23, :black))
annotate!(p3, 30, 1.45, text("raw data", :center, 15, :black))
p4 = plot(v2', c=cs', xlims=(1,50), legend=false, title="", ylims=(0,1.5), yticks=(0:0.5:1, 0:0.5:1), ylabel="g(z,t)",  xticks=(0:45:135, round.(Int,0:45/15:135/15) ), xlabel="Time (s)" )
annotate!(p4, -15, 1.45, text("C2", :center, 23, :black))
annotate!(p4, 30, 1.45, text("5 frames (0.3 s) ", :center, 15, :black))


p = plot(p3, p4, layout=(1,2))
decay = plot(p0,p1,p,layout=(3,1))

plot!(size=(600,800), leftmargin=5Plots.mm)
savefig(decay, plotsdir("Fig4.pdf"))
decay


## Main text Figure 5
PALETTE = 9
cs = cgrad(:roma, 1:PALETTE)[collect(1:PALETTE) ./PALETTE]
mat = inverse_fits[2,:, 3, :]
p = plot()
for p in 2:9
    # plot!(reverse(1:42), 100*mat[:,p], c=cs[p], legend=false)
    scatter!(reverse(1:42), 100*mat[:,p], c=cs[p], legend=false, msc=cs[p], ms=6)
    # plot!(ylims=(0,0.20))
end
fit_all = plot!(ylabel="ν"*L"_{G} (s^{-1})", xlabel="z (mm)", title=" ")
annotate!((1.2, 17.9, text(" x 10²")))

# plot!(xticks=(0:10:40, 0:4:16)) ## z (mm)
plot!(xticks=(0:10:40, 0:10:40), xlabel="z (mm)") ## z (mPixel)
heatmap!(speeds_z, c=cs[1:9], colorbar=false, title="ω (s⁻¹)", yticks=:none, titlefontsize=13,
    inset_subplots = bbox(0.65, 0.8, 0.3, 0.07, :bottom), subplot=2, axes=false, xticks=(1:2:10, speeds[[2:2:9;9]]), xrotation=-30,
)
savefig(p, plotsdir("Fig5.pdf"))
p

## Main text Figure 6/8
fits_data.inverse_fits
gr()
ss = parse.(Float64,speeds[2:end])
mat = inverse_fits[2,:, 3, 2:end]
p1 = heatmap(speeds[2:end], 1:42, 100*mat, c=:redsblues,
             xlabel="ω (s⁻¹)", ylabel="z (mm)", yticks=(0:10:40, reverse(0:4:16)), 
            colorbar_title="\nν"*L"_{G} (s^{-1})", cbarfontsize=12, title=" ", margin=5Plots.mm, rightmargin=14Plots.mm, xrotation=-45)
annotate!((10, 45, text(" x 10²")))
annotate!(-1.2, 50, text("A", :center, 23, :black))

ss = parse.(Float64,speeds[2:end])
mat = inverse_fits[4,:, 5, 2:end]
p2 = heatmap(speeds[2:end], 1:42, mat, c=:redsblues,
             xlabel="ω (s⁻¹)", ylabel="z (mm)", yticks=(0:10:40, reverse(0:4:16)), 
            colorbar_title="\nν"*L"_{H} (s^{-1})", cbarfontsize=12, title=" ", margin=5Plots.mm, rightmargin=14Plots.mm, xrotation=-45)
annotate!(-1.2, 50, text("B", :center, 23, :black))
# annotate!((10, 45, text(" x 10²")))
p = plot(p1,p2, layout=(2,1), size=(600,800))
savefig(p, plotsdir("Fig6.pdf"))
p

## Main text Figure 7
my_speeds = [parse(Float64,x) for x in speeds]
@unpack func, coarse_time= config
t = findfirst(x->x==coarse_time, times)
f = findfirst(x->String(x)==func, functions)
cs = cgrad(:roma, 1:PALETTE)[collect(1:PALETTE) ./PALETTE]
taus = [inverse_fits[f, 42-round(Int,shearband.max[x]),t,x+1] for x in 1:8]
xs_max = [shearband.max[x]*metapixel_size for x in 2:9]
xs_min = [shearband.min[x]*metapixel_size for x in 2:9]
p0 = scatter(my_speeds[2:end],100 .* taus[1:end], label="",c=:black, ylabel=L"\nu_G (s^{-1})", ms=12, xticks=:none, ylims=(5,25), topmargin=10Plots.mm,)
annotate!(-0.55, 27.5, text("A", :center, 23, :black))
annotate!(-0.55, 2, text("B", :center, 23, :black))
scatter!([],[],shape=:diamond, label=" "*L"Z_{max}", c=:black, msc=:black, ms=12)
scatter!([],[],shape=:star, label=" "*L"Z_{B}", c=:black, msc=:black, ms=12, fg_legend=:transparent, legend=:topleft)
scatter!([],[],shape=:circle, label=" "*L"\nu^{max}_G (s^{-1})", c=:black, msc=:black, ms=12, fg_legend=:transparent, legend=:topleft,)
p1 = scatter(my_speeds[2:end], xs_max, label="", c=:black, ms=12, shape=:diamond, ylims=(0,16))
scatter!(my_speeds[2:end], xs_min, label="", yerror=shearband.err[2:end]*metapixel_size, c=:black, ylabel="z (mm) ", ms=12, shape=:star, xlabel="speed "*"ω (s⁻¹)", ylims=(0,15))
p = plot(p0,p1, layout=(2,1), xrotation=-45, size=(800,400),  rightmargin=14Plots.mm, xlims=(-0.2, 1.7), bottommargin=10Plots.mm, leftmargin=10Plots.mm,
legendfontsize=16)
plot!(size=(600,600))
savefig(p, plotsdir("Fig7.pdf"))
p
##

## Figure 9
PALETTE = 9
cs = cgrad(:roma, 1:PALETTE)[collect(1:PALETTE) ./PALETTE]
mat = inverse_fits[2,:, 5, :]
inverse_fits
gg_max = []
hh_max = []
p = plot()
shift=0.1
ylims=(0.1, 1.2)

for speed in eachindex(speeds)[2:end]
    hline!([shift*speed], label="", c=:black, alpha=0.4, ls=:dash)
end
for speed in eachindex(speeds)[2:end]
    # println(speed)
    hs =  inverse_fits[4,:,1,speed]
    gs =  inverse_fits[2,:,3,speed]
    push!(hh_max,(43- argmax(hs))*metapixel_size)
    push!(gg_max,(43- argmax(gs))*metapixel_size)
    scatter!(hs,gs.+ shift*speed, label="", xlabel="Half width "*L"\nu_H (s^{-1})", ylabel="Gaussian exp "*L"\nu_H (s^{-1})", c =cs[speed],  msc=:auto, ms=8, ylims=ylims)
    scatter!([maximum(hs)],[maximum(gs)].+shift*speed, label="",shape=:square, c =:darkred, ms=4, msc=:auto, yticks=:none, leftmargin=10Plots.mm)

end

plot!(twinx(p), ylims=ylims, yticks=(0.2:0.1:0.9,speeds[2:end]), ylabel="ω (s⁻¹)", yrotation=-45, )
p

savefig(p, plotsdir("SIFig3.pdf"))

p =plot([0,10],[0,10], label="", c=:black, ls=:dash)
scatter!(hh_max, gg_max, label="", xlabel="H", ylabel="G", ylims=(10,17), xlims=(10,17), c=cs, msc=:auto, ms=10)
plot!(ylabel=L"Z_{max} "*" from "*L"\nu_H"*"(mm)", xlabel=L"Z_{max} "*" from "*L"\nu_G"*"(mm)", ylims=(0,8), xlims=(0,8), c=:black, ls=:dash, margin=5Plots.mm, rightmargin=14Plots.mm, legend=:topleft, bottommargin=10Plots.mm, topmargin=10Plots.mm)

heatmap!(speeds_z, c=cs[1:9], colorbar=false, title="ω (s⁻¹)", yticks=:none, titlefontsize=13,
    inset_subplots = bbox(0.2, 0.78, 0.3, 0.07, :bottom), subplot=2, axes=false, xticks=(1:2:10, speeds[[2:2:9;9]]), xrotation=-30,
)
savefig(p, plotsdir("Fig9.pdf"))
p
# scatter(hs[2:end],gs[2:end], label="", xlabel="H", ylabel="G")
##

# heatmap!(z, c=cs[1:end], colorbar=false, title="ω (s⁻¹)", yticks=:none, titlefontsize=13,
#     inset_subplots = bbox(0.65, 0.8, 0.3, 0.07, :bottom), subplot=2, axes=false, xticks=(1:2:10, speeds[[2:2:9;9]]))
# savefig(p, joinpath(results_path,"fig7_maxtau.pdf"))

# cs = cgrad(:roma, 1:PALETTE)[collect(1:PALETTE) ./PALETTE]
# taus_x = [42 - argmax(inverse_fits[:,2,x,5]) for x in 1:8]
# taus = 100 .*[maximum(inverse_fits[:,2,x,5]) for x in 1:8]
# p = scatter(taus_x, taus, label="", c=cs, msc=cs, ylabel="maximum ν (s⁻¹)", ms=12, xlabel="z (mm)")
# plot!(xticks=(0:10:20, 0:4:8), legend=:topleft, xlims=(0:8))
# heatmap!(z, c=cs[1:end], colorbar=false, title="ω (s⁻¹)", yticks=:none, titlefontsize=13,
#     inset_subplots = bbox(0.65, 0.8, 0.3, 0.07, :bottom), subplot=2, axes=false, xticks=(1:2:10, speeds[[2:2:9;9]]))
# savefig(p, joinpath(results_path,"fig7_maxtauinv.pdf"))
# plot!(title=" ")
# annotate!((0.65, 6.45, text(" x 10²")))

# # scatter(1 ./inverse_fits[:,2,:,1])


##
# p = groupedbar(xticks=(1:8,string.(speeds)),[-shearband.max -shearband.min], bar_position=:stack,c=[:white cs[6] ], lc=:white, legend=false, ylabel= "z (mm)", xlabel="ω (s⁻¹)", guidefontsize=18, tickfontsize=13, grid=false, ylims=(-42,0),yticks=(reverse(-40:10:0), 0:4:16))
# savefig(p, joinpath(results_path,"fig8_bars.pdf"))
# ##
# gr()
# layout = @layout [
# 			a{0.9w} _;_ c{0.0001h}
# 			]
# p = plot(parse.(Float64,speeds[2:end]), shearband.max, seriestypes=[:line, :scatter], c=:red, msc=:red , ms=8, labels=["" "z maximum"],  ylabel="z maximum (mm)", legend=:topleft)
# p = plot!(twinx(), parse.(Float64,speeds[2:end]), shearband.min, ylabel="z shear (mm)", label=["" "z shear"], seriestypes=[:line, :scatter], c=:black, msc=:black, ms=8,  legend=:bottomright,xlabel="ω (s⁻¹)")
# p =plot(p,plot(frame=:none),layout=layout)
# savefig(p, joinpath(results_path,"fig8.pdf"))

# ##
# p = scatter(parse.(Float64,speeds[2:end]), maxima[2,:,5], label="Gaussian Exp. time avg", xlabel="ω (s⁻¹)", ylabel="z (mm)", c=:black, ms=8, shape=:square)
# p = scatter!(parse.(Float64,speeds[2:end]), maxima[4,:,1], label="Half width Raw data", xlabel="ω (s⁻¹)", c=:red, msc=:red, ms=8)
# plot!(yticks=(30:5:40, reverse(0:2:4)), legendfontsize=12, title=" ", )


# savefig(p, joinpath(results_path,"fig9.pdf"))
# p
# ##
# p = scatter(maxima[2,:,5], maxima[4,:,1], ylabel="Half width Raw data", xlabel="Gaussian Exp. time avg", c=:black, msc=:black, ms=8, label="")

# savefig(p, joinpath(results_path,"fig9_regression.pdf"))


##
PALETTE = 9
cs = cgrad(:roma, 1:PALETTE)[collect(1:PALETTE) ./PALETTE]

shearband.err
s=4
plot()
for s in 2:9
    plot!(reverse(inverse_fits[2,:,3,s]), c=cs[s])
    scatter!([shearband.max[s]],[inverse_fits[2,43-round(Int,shearband.max[s]),3,s]], c=cs[s], shape=:diamond, ms=12, msc=:auto)
    scatter!([shearband.min[s]],[inverse_fits[2,43-round(Int,shearband.min[s]),3,s]], c=cs[s], shape=:star, ms=12, msc=:auto)
end 
p = plot!(legend=false, xlabel="Metapixel", ylabel=L"\nu_G (s^{-1})")
savefig(p, plotsdir("SIFig2.pdf"))
p

##