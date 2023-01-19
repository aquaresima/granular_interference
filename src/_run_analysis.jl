include("Speckles.jl")
using .Speckles
using HDF5
using LinearAlgebra
using Plots
using StatsBase
using ColorSchemes
using Printf
using LaTeXStrings
using RollingFunctions
using ProgressBars

default(guidefontsize=18, tickfontsize=13, titlefontsize=15)
#%%md Load files

# include("make_matrices.jl")
full_matrix = "/run/media/cocconat/data/granulari/full_matrices"
analysis_path = "/run/media/cocconat/data/granulari/analysis/"
results_path = "final_plots"
!isdir(results_path) && mkdir(results_path)

temporal_file = joinpath(analysis_path,"coarsen10_time.h5")
corr_file = joinpath(analysis_path,"coarsen10_correlation.h5")
norm_file = joinpath(analysis_path,"norm.h5")
times = ["1", "2","5", "10", "20", "40"];

_norm = read(h5open(norm_file,"r")["norm"])
rawdata_file = h5open(joinpath(analysis_path,"coarsen10.h5"),"r")

all_matrix = read(rawdata_file)
speeds = sort(collect(keys(all_matrix)))

#=================================================
                    Analysis
=================================================#

## Plot average pixel amplitude
plots = []



(root,dirs,files) =  first(walkdir(full_matrix))
for (file,speed) in collect(zip(files, speeds))[[1,2,5,9]]
	_mat = read(h5open(joinpath(root,file),"r"))["matrix"]
    flat = StatsBase.mean(_mat, dims=3)[:,:,1]
    push!(plots, heatmap(flat, title="Speed:"*speed*L"\omega / s",colorbar=false, c=:haline,titlefontsize=10, xticks=nothing))
end


plot!(plots[4], ylabel="Depth (px)")
mean_pixels = plot(plots...);
savefig(mean_pixels,joinpath(results_path,"fig2.pdf"))
plot!()

##
# Verify that the image has an exponential distribution in light-intensity to verify the measure regard actual spekles. : Goodman Statistical Optics
speeds
l = @layout [
    b{0.5w} c{0.5w}
	d{0.5w} e{0.5w}
]
plot(histogram( all_matrix["0.0"][:], norm=true , alpha=1., bins=-0:5:200, title="speed: 0.0"),
histogram( all_matrix["1.0"][:], norm=true, alpha=1., bins=0:2:200      , title="speed: 0.5"),
histogram( all_matrix["0.5"][:], norm=true, alpha=1., bins=0:2:200      , title="speed: 1.0",xlabel="Pixel Intensity", ylabel="Pixel density", ),
histogram( all_matrix["1.5"][:], norm=true, alpha=1., bins=0:2:200      , title="speed: 1.5"),
layout=l, lw=3)

intensity = plot!( legend=false)
savefig(intensity, joinpath(results_path,"fig_intensity.pdf"))
# plot!(title="Light ")
plot!(guidefontsize=18, tickfontsize=13)

## Plot vertical correlation matrices

# vcorr_file = joinpath(analysis_path,"coarsen10_vertical_correlation.h5")
#
# vcorr_mats = read(h5open(corr_file))
#
# for T in [1,10,50,100,300]
#     c = cgrad(:roma, 1:30)[0.03:0.03:1]
#     plots = []
#     for (n, speed) in enumerate(speeds)
#         push!(plots,plot(corr_mats["10"][n][:,:,T]',c=c', labels=false,  colorbar=true, title="Speed:"*speed,titlefontsize=11))
#     end
#     plot!(plots[8],xlabel="depth (pixel)")
#     plot!(plots[4],ylabel="Cross-correlation" )
#     # heatmap(zeros(30,2),c=c, clims=(1,30))
#     p = plot(plots...)
#     savefig(p,"Vertical_corr"*string(T)*".pdf")
# end
##
corr_mats = read(h5open(corr_file))
c = cgrad(:roma, 1:9)[1/9:1/9:1]
# plots = []
# for t in times
p = plot()
speeds
for (n,speed) in enumerate(speeds[2:end])
	v = StatsBase.mean(corr_mats["1"][speed][20,:,1:300], dims=1)[:]
	p=plot!(1:300,v, label=false,  xlims=(1,30), c=c[n], lw=3)
end
p = plot(p, xlabel="Time (s)", ylabel="Correlation", lw=3, xticks=(0:7.5:30, 0:0.5:2), guidefontsize=18, tickfontsize=13,titlefontsize=15, title="Correlation decay depth 8mm", grid=false)
savefig(p, joinpath(results_path,"fig3.pdf"))
p
##
    # for (n, speed) in enumerate(speeds[1:end])
    #         # for (z,T) in enumerate(1:300)
    #         # s = parse(Float64, speed)
    #         # p = scatter!(p,ones(30)*s,corr_mats["10"][n][:,1,z], label=false)
    #         if t !== "40"
    #         elseif t == "40"
    #             p=plot!(1:300,v,  xlims=(1,30), c=c[n], label=speed)
    #         end
    #         plot!(title="Temporal mean: $t", titlefontsize=11)
    #         # ,c=c', labels=false,  colorbar=true, title="Speed:"*speed,titlefontsize=11))
    # end
    # plot([])
# plot!(plots[5], xlabel="Auto-Correlation time (frames)")
# plot!(plots[4], ylabel="Correlation")
# annotate!()


p = plot(plots...)
plot!(guidefontsize=18, tickfontsize=13, xrotation=45, grid=false)

savefig(p, joinpath(results_path,"Sample_correlation_row20.pdf"))
# , joinpath(results_path,"Sample_correlation_row20.pdf"))
plot!()

##
#
z = ones(1,42)
for x in 1:42
    z[x] = 42-x
end
PALETTE = 42
cs = cgrad(:roma, 1:PALETTE)[collect(1:PALETTE) ./PALETTE]
v = StatsBase.mean(corr_mats["1"]["0.05"][:,:,1:300], dims=2)[:,1,:]

p1 =plot(v', c=cs', xlims=(1,135), legend=false, title="Raw correlations", ylabel="Correlation", xaxis=false, ylims=(0,1.6), yticks=(0:0.2:1, 0:0.2:1) )

heatmap!(z, c=cs, colorbar=false, ylabel="", yticks=:none,
    inset_subplots = bbox(0.55, 0.6, 0.4, 0.1, :bottom), xticks=(0:10:40, 0:4:16), subplot=2, title= "Pixel Depth (mm)")

v = StatsBase.mean(corr_mats["20"]["0.05"][:,:,1:300], dims=2)[:,1,:]
p2 = plot(v', c=cs', xlims=(1,135), legend=false, title="Time averaged", ylims=(0,1.5), yticks=(0:0.2:1, 0:0.2:1), xticks=(0:45:135, round.(0:45/15:135/15, digits=2) ), xlabel="Time (s)" )
# heatmap!(z, c=cs, colorbar=false, ylabel="", yticks=:none,
    # inset_subplots = bbox(0.55, 0.8, 0.4, 0.1, :bottom), subplot=2,
# )
decay = plot(p1,p2,layout=(2,1))

plot!(guidefontsize=18, tickfontsize=13, xrotation=45, grid=false)
savefig(decay, joinpath(results_path,"fig4.pdf"))
## Fit exponentials
using LsqFit
using HypothesisTests
@. gauss_exp(x,p) =exp(-(x/p[1])^2)
@. stretch_exp(x,p) = exp(-(x/abs(p[1]))^abs(p[2]))
@. fit_cosh(x,p) = 2-cosh(x/p[1])
function halfheight(data)
        x = findfirst(x->x<data[1]/2., data)
        return  isnothing(x) ? [length(data)] : [x]
end


##
"""
Fit decays correlation with three curves:
	- Gaussian_Exponential
	- Hyperbolic Cosine
	- Stretched exponential
Parameters:
----------
	matrix: correlation matrix
	y: pixel depth
	px: number of points per fit

Return
------

"""

function makefit_example(;y, tt, px) # depth, average, points
    mysample = mean(corr_mats[string(tt)]["1.0"], dims=2)[:,1,1:80]
    a  = curve_fit(gauss_exp, 1:px, mysample[y,1:px], [1.]).param
    b  = curve_fit(fit_cosh, 1:px, mysample[y,1:px], [1.]).param
    c  = curve_fit(stretch_exp, 1:px, mysample[y,1:px], [3.,1.]).param
    @show c
    plot(mysample[y,1:px], label="data", lw=3)
    plot!(gauss_exp.(1:px,a), label="gauss", lw=3, ls=:dot)
    plot!(fit_cosh.(1:px,b), label="cosh", lw=3, ls=:dash)
    plot!([stretch_exp(x,c) for x in 1:px], label="stretch", lw=3, ls=:dash)
    # y2 = halfheight(sample[y,:])
    # scatter!([y2],[sample[y,y2] ], label="Half height", ms=10)
    fit_sample = plot!(legendfontsize=13, legend=:bottomleft, title= "avg: $tt, px: $px, y: $y")
    Eg = gauss_exp.(1:px,a)
    Ec = fit_cosh.(1:px,b)
    Es = [stretch_exp(x,c) for x in 1:px]
    O = mysample[y,1:px]
    Xc = sum((O-Ec).^2)/length(Ec)
    Xg = sum((O-Eg).^2)/length(Eg)
    Xs = sum((O-Es).^2)/length(Es)
    return fit_sample, Xc, Xg, Xs
     # O, Eg, Ec
end

function get_best_fit(y)
    plots=[]
    pxs = [5,10,20] #points
    tts = [1, 5,20,40] #coarse
    X = zeros(3, 4,  3) # pxs, tts, depth
    for n in 1:length(pxs)
        for m in 1:length(tts)
            px = pxs[n]
            tt = tts[m]
            # p, Xc, Xg = makefit_example(y=10, tt=tt, px=px)
            # X[n,m,1,1] = Xg
            # X[n,m,1,2] = Xc
            # push!(plots,p)
            p, Xc, Xg, Xs = makefit_example(y=y, tt=tt, px=px)
            X[n,m,1] = Xg
            X[n,m,2] = Xc
            (n*m > 1) && (plot!(p, legend=false))
            push!(plots,p)
            # p, Xc, Xg = makefit_example(y=30, tt=tt, px=px)
            # X[n,m,3,1] = Xg
            # X[n,m,3,2] = Xc
        end
    end
    return plot(plots..., titlefontsize=8,legendfontsize=6)
end

p = get_best_fit(10)
savefig(p, joinpath(@__DIR__,results_path,"bestfit_depth: 10.pdf"))
p = get_best_fit(20)
savefig(p, joinpath(@__DIR__,results_path,"bestfit_depth: 20.pdf"))
p =get_best_fit(30)
savefig(p, joinpath(@__DIR__,results_path,"bestfit_depth: 30.pdf"))
p =get_best_fit(40)
savefig(p, joinpath(@__DIR__,results_path,"bestfit_depth: 40.pdf"))

# plot(heatmap(X[:,2,:,1], title="Gaussian"), heatmap(X[:,2,:,2], title="Cosh"))
##
function fit_correlations(matrices, coarse_t::String)
    corr_mat = corr_mats[coarse_t]
    tt = parse(Int,coarse_t)
    tt = 10 ## HARDCODED
    speed_fits = Dict()
    corr_depths = []
    for s in speeds[2:end]
        speed_fits[s] = []
        mysample = mean(corr_mats[coarse_t][s], dims=2)[:,1,:]
        for y in 1:size(mysample)[1] # number of rows
            ## run all fits
            fits = (
            stretched = curve_fit(stretch_exp, 3tt:4tt, mysample[y,tt:2tt], [tt,1.]).param,
            gaussian  = curve_fit(gauss_exp, 1:tt, mysample[y,1:tt], [1.]).param,
            cosh = curve_fit(fit_cosh, 1:tt, mysample[y,1:tt], [1.]).param,
            half      = halfheight(mysample[y,:])
            )
            push!(speed_fits[s],fits)
        end
    end
    return speed_fits
end


times_fits = Dict()
for t in times
    speed_fits = fit_correlations(corr_mats, t)
    push!(times_fits, (t=>speed_fits))
end
##
times_fits
##
times_fits
plots = []
PALETTE = 9
cs = cgrad(:roma, 1:PALETTE)[collect(1:PALETTE) ./PALETTE]
p1 = plot()
p2 = plot()
p3 = plot()
p4 = plot()
coarse = "1"
for (n,s) in enumerate(speeds[2:end])
    p1 = scatter!(p1,[x.half for x in times_fits["1"][s]],    [x.half for x in times_fits[coarse][s]], label=false, c=cs[n], title="half")
    p2 = scatter!(p2,[x.stretched for x in times_fits["1"][s]],    [x.stretched for x in times_fits[coarse][s]], label=false, c=cs[n], title="stretched")
    p3 = scatter!(p3,[x.gaussian for x in times_fits["1"][s]],    [x.gaussian for x in times_fits[coarse][s]], label=false, c=cs[n], title="gaussian")
    p4 = scatter!(p4,[x.cosh for x in times_fits["1"][s]],    [x.cosh for x in times_fits[coarse][s]], label=false, c=cs[n], title="cosh")
end
p = plot(p1,p2,p3,p4, xlabel="Raw data", ylabel="Time-averaged data")

savefig(p, joinpath(@__DIR__,results_path,"correlation raw-coarse $coarse.pdf"))
p# plot(plots...)

## Correlate stretch fit of raw data with gaussian on averaged
times_fits
plots = []
PALETTE = 9
cs = cgrad(:roma, 1:PALETTE)[collect(1:PALETTE) ./PALETTE]
p1 = plot()

coarse = "1"
for (n,s) in enumerate(speeds[2:end])
    p1 = scatter!(p1,[x.stretched[1:1] for x in times_fits["1"][s]],    [x.gaussian for x in times_fits[coarse][s]], label=false, c=cs[n], title="half")
end
plot!(p1, ylabel="Gaussian coarse grain", xlabel="Stretched Raw")
# p = plot(p1,p2,p3,p4, xlabel="Raw data", ylabel="Time-averaged data")
#
# savefig(p, joinpath(@__DIR__,results_path,"correlation raw-coarse $coarse.pdf"))

##
PALETTE = 9
cs = cgrad(:roma, 1:PALETTE)[collect(1:PALETTE) ./PALETTE]
time_plots = []
for t in times
    s_p = plot(legend=false, title="Stretched")
    g_p = plot(legend=false, title="Gaussian")
    c_p = plot(legend=false, title="Cosh")
    h_p = plot(legend=false, title="Half")
    for (s,c) in zip(speeds[2:end],cs)
        data = times_fits[t][s]
        for (n,d) in enumerate(data)
            scatter!(s_p, [42-n],d.stretched[1:1] , c=c, ylabel="stretch fit")
            scatter!(g_p, [42-n],d.gaussian , c=c, ylabel="gauss fit")
            scatter!(c_p, [42-n],d.cosh , c=c, ylabel="cosh fit")
            scatter!(h_p, [42-n],d.half , c=c, ylabel="half-height time")
        end
    end
    push!(time_plots,(c_p,g_p, h_p, s_p))
end
# plot(time_plots...)


##

z = ones(1,9)
for x in 1:9
    z[x] = x
end
fit_plots = plot(time_plots[1]..., xlabel="(px)")
heatmap!(z, c=cs, colorbar=false, xlabel="speed", yticks=:none,
    inset_subplots = bbox(0.1, 0.85, 0.15, 0.07, :bottom), subplot=5, axes=false
)
savefig(fit_plots, joinpath(@__DIR__,results_path,"Fit_all_coarse20.pdf"))
plot!(fit_plots)
## Find maxima

maxima = Dict()
times
fits_values = zeros(42,4,8,6) # rows, function, speed, coarse
maxima = zeros(4,8,6)
for (k,t) in enumerate(times)
    for (j,s) in enumerate(speeds[2:end])
        data = times_fits[t][s]
        for (m,f) in enumerate([:stretched, :gaussian, :cosh, :half] )
            values = Vector{Real}()
            for (n,d) in enumerate(data)
                push!(values,1. /getfield(d,f)[1])
                fits_values[n,m,j,k] = 1. / getfield(d,f)[1]
            end
            # values = runmean(values,10)
            maxima[m,j,k] = argmax(values)
        end
    end
end

maxima

##
plot()
PALETTE = 9
cs = cgrad(:roma, 1:PALETTE)[collect(1:PALETTE) ./PALETTE]
plots =[plot(),plot(),plot(),plot()]

coarse = 1
for p in 1:8
    for i in 1:4
        plot!(plots[i],reverse(1:42), fits_values[:,i,p,coarse], c=cs[p], legend=false)
        scatter!(plots[i],reverse(1:42), fits_values[:,i,p,coarse], c=cs[p], legend=false)
        plot!(ylims=(0,0.20))
        # plot!(plots[i],reverse(1:42), runmean(fits_values[:,i,p,3]/ maximum(fits_values[:,i,p,3]),1), c=cs[p], legend=false)
    end
end
fit_all = plot(plots..., title=["Stretched" "Gaussian" "Cosh" "Halfheight"], ylabel="tau ^-1", xlabel="Depth (px)")

heatmap!(z, c=cs, colorbar=false, xlabel="speed", yticks=:none,
    inset_subplots = bbox(0.3, 0.85, 0.15, 0.07, :bottom), subplot=5, axes=false
)
savefig(fit_all, joinpath(@__DIR__,results_path,"Fit_all_coarse$(times[coarse])_inversetau.pdf"))
plot(fit_all)
##
function plot_fits(t)
    PALETTE = 4
    cs = cgrad(:roma, 1:PALETTE)[collect(1:PALETTE) ./PALETTE]
    plot(parse.(Float64,speeds[2:end]), maxima[1,:,t] , label="Stretched", c=cs[1])
    plot!(parse.(Float64,speeds[2:end]), maxima[2,:,t], label="Gaussian", c=cs[2])
    plot!(parse.(Float64,speeds[2:end]), maxima[3,:,t], label="Cosh", c=cs[3])
    plot!(parse.(Float64,speeds[2:end]), maxima[4,:,t], label="Half", c=cs[4])
    scatter!(parse.(Float64,speeds[2:end]), maxima[1,:,t], c = cs[1], label="")
    scatter!(parse.(Float64,speeds[2:end]), maxima[2,:,t], c = cs[2], label="")
    scatter!(parse.(Float64,speeds[2:end]), maxima[3,:,t], c = cs[3], label="")
    scatter!(parse.(Float64,speeds[2:end]), maxima[4,:,t], c = cs[4], label="")
end

plots =[]
for (t, tt) in enumerate(times)
    if t > 1
        p = plot!(plot_fits(t), title="Time Coarse $tt", legend=false)
    else
        p = plot!(plot_fits(t), title="Time Coarse $tt", legendfontsize=10, legend=:bottomright)
    end
    plot!(p, ylabel="Depth of maximum", xlabel="Speed (rad/s)", xrotation=45, guidfontsize=14)
    push!(plots,p)
end

final_fit = plot(plots..., legendfontsize=5)

savefig(final_fit, joinpath(@__DIR__,results_path,"Depth_vs_speed.pdf"))

plot(final_fit)

fits_values  # rows, function, speed, coarse

##
PALETTE = 8
cs = cgrad(:roma, 1:PALETTE)[collect(1:PALETTE) ./PALETTE]
plots =[plot(),plot(),plot(),plot()]
plot()

func = 2
mobilmean = 15
coarse = 5

banda = zeros(Int,8,3)
plot()
for (s, c) in zip(1:8,cs)
    new_val = fits_values[:,func,s,coarse]
    plot!(reverse(1:42), new_val, c=c, lw=2)
    # plot!(reverse(1:42), runmean(new_val, mobilmean), c=c, lw=2)
    deriv = diff(reverse(fits_values[:,func,s,coarse])*10)
    mymin = argmin(runmean(deriv,mobilmean))
    mymax = argmax(reverse(fits_values[:,func,s,coarse]))
    scatter!([mymin], [reverse(new_val)[mymin]], c=c, ms=8)
    # scatter!([mymax], [reverse(new_val)[mymax]], c=c, ms=8)
    banda[s,1] = mymin
    banda[s,2] = mymax
    # banda[s,3] = (mymax-mymin)/2 + mymin
end
plot!(legend=false)


##
plot()
for (s, c) in zip(1:8,cs)
    plot!(reverse(1:42),fits_values[:,func, s, coarse], c=c, lw=2)
   scatter!([banda[s,1]],[fits_values[43-banda[s,1], func, s,coarse]],ms=15, shape=:diamond, c=c)
   scatter!([banda[s,2]],[fits_values[43-banda[s,2], func, s,coarse]],ms=15, shape=:star5, c=c)
end
q = plot!(legend=false)
savefig(q,"max_min$(times[coarse]).pdf")
q
##


using StatsPlots

p = groupedbar(xticks=(1:8,string.(speeds)),[-banda[:,2] -banda[:,1]], bar_position=:stack,c=[:white cs[6] ], lc=:white, legend=false, ylabel= "shear band", xlabel="Plate speed (ω/s)", guidefontsize=18, tickfontsize=13)

savefig(p,"shearband$(times[coarse]).pdf")

p
