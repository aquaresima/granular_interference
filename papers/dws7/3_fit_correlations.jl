using DrWatson
@quickactivate "granular_media"

using Revise
using speckles
import YAML
using JLD2
using UnPack

using LsqFit
using HypothesisTests
using Statistics
using RollingFunctions

# Set the location of the folder and import all files.
conf = YAML.load_file(projectdir("conf.yml"))
full_matrices  = conf["paths"]["full_matrices"]
analysis_path  = conf["do_norm"] ? conf["paths"]["analysis"] : conf["paths"]["analysis"] * "_no_norm" 
@assert isdir(full_matrices) 
@assert isdir(analysis_path)

# Set coarse dimensions
coarse = conf["coarse"]
x = coarse["x"]
y = coarse["y"]
times = coarse["times"]

# Set filename]
filename = "coarsen$(x)x$(y)"
corr_file = "$(filename)_correlations.jld2"
temporal_file(root,t) = "$(root)_$t.jld2"
temporal_file(filename,10)

## Fit exponentials
functions = [:stretched, :gaussian, :cosh, :half] 
N_points = 10


@. gauss_exp(x,p) =exp(-(x/p[1])^2)
@. stretch_exp(x,p) = exp(-(x/abs(p[1]))^abs(p[2]))
@. fit_cosh(x,p) = 2-cosh(x/p[1])
function halfheight(data)
        x = findfirst(x->x<data[1]/2., data)
        x = isnothing(x) ? [length(data)] : [x]
        return  x/15 # 15 fps, return the time in seconds
end
function get_fits(data)

    global N_points
    tt = N_points
    return (
    stretched = curve_fit(stretch_exp, 3tt:4tt, data[tt:2tt], [tt,1.]).param,
    gaussian  = curve_fit(gauss_exp,   1:tt, data[1:tt], [1.]).param,
    cosh      = curve_fit(fit_cosh,    1:tt, data[1:tt], [1.]).param,
    half      = halfheight(data)
    )
end

filename
joinpath(analysis_path, corr_file)
data = load(joinpath(analysis_path, corr_file)) |> tosymboldict |> dict2ntuple;
data
config = @strdict data=data functions=functions
fits,_ = produce_or_load(path, config,filename="$(filename)_fits") do config
    data = config["data"]
    functions = config["functions"]
    @unpack correlations, speeds, times = data 
    ys = size(correlations[1,1],1)
    inverse_fits = Array{Float64,4}(undef, size(functions,1), ys,size(times,1), size(speeds,1))
    fits = Array{Float64,4}(undef, size(functions,1), ys,size(times,1), size(speeds,1))
    for s in eachindex(speeds)
        for t in eachindex(times)
            mysample = mean(correlations[t,s], dims=2)[:,1,:]
            @info "Run fit $(speeds[s]) $(times[t]), matrix: $(size(mysample))"
            for y in 1:size(mysample)[1] # number of rows
                _fits = @views get_fits(mysample[y,:])
                ## run all fits
                # fits[t,s,y] = _fits
                inverse_fits[:,y, t,s] = [ 1. / getfield(_fits,f)[1] for f in functions]
                fits[:,y, t,s] = [  getfield(_fits,f)[1] for f in functions]
            end
        end
    end
    return @strdict fits inverse_fits times speeds
end

fits["fits"]

data = load(joinpath(analysis_path,"$(filename)_fits.jld2")) |> tosymboldict |> dict2ntuple;
@unpack fits, inverse_fits, times, speeds = data

## Get bandwidth for specific values
fit_params = conf["fit"]
func = fit_params["func"]
coarse_time = fit_params["coarse_time"]
config = (@strdict func coarse_time) |> dict2ntuple

name = savename(filename, config, "jld2")
isfile(name) && rm(name)
data,_ = produce_or_load(analysis_path,config,prefix=filename) do config
    @unpack func, coarse_time = config
    band = zeros(9,3)
    for s in 1:9
        t = findfirst(x->x==coarse_time, times)
        f = findfirst(x->String(x)==func, functions)
        local data = inverse_fits[f,:,t,s]

        mymax = argmax(rollmean(reverse(data),1))
        band[s,1] = mymax

        means = []
        for mobilmean in 4:15
            zz =rollmean(reverse(data),mobilmean)
            my_diff = diff(zz[argmax(zz)+1:end])
            push!(means,argmax(my_diff)+argmax(zz)+1)
        end
        band[s,2] = mean(means)
        band[s,3] = std(means)
    end
    shearband =(min=band[:,2], max=band[:,1], err=band[:,3])
    @show shearband.min
    return @strdict shearband=shearband func=func coarse_time=coarse_time mobilmean=mobilmean
end
@unpack shearband = data

##
plot()
for s in 2:9
    zz =rollmean(reverse(inverse_fits[2,:,3,s]),8)
    my_diff = diff(zz[argmax(zz)+1:end])
    zz_diff = argmax(my_diff)
    @show  maximum(diff(zz[argmax(zz)+1:end]))
    plot!(zz,c=cs[s])
    @show zz_diff
    scatter!([zz_diff+argmax(zz)+1],[zz[zz_diff+argmax(zz)]],c=cs[s])
end
plot!(legend=false)
##



##
PALETTE = 9
cs = cgrad(:roma, 1:PALETTE)[collect(1:PALETTE) ./PALETTE]
plot()
for s in 1:9
    scatter!([x.half for x in fits[3,s,:]],
        [x.gaussian for x in fits[3,s,:]], label="",c=cs[s],msc=:auto,ms=4)
end
plot!()

# label="half")


tt= 10
gaussian  = curve_fit(gauss_exp,   1:tt, data.correlations[4][1:tt], [1.])

data.correlations[4][1:tt]-gauss_exp.(1:tt, gaussian.param[1])