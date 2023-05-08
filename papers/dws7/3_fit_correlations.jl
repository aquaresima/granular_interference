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
frame_rate = conf["frame_rate"]

# Set filename]
filename = "coarsen$(x)x$(y)"
corr_file = "$(filename)_correlations.jld2"
temporal_file(root,t) = "$(root)_$t.jld2"
temporal_file(filename,10)

## Fit exponentials
functions = [:stretched, :gaussian, :cosh, :half] 
N_points = 10


τ_window = 5
@. gauss_exp(x,p) =exp(-(2x/sqrt(p[1]*τ_window))^2)
@. stretch_exp(x,p) = exp(-(x/abs(p[1]))^abs(p[2]))
@. fit_cosh(x,p) = 2-cosh(x/p[1])
function halfheight(data)
        x = findfirst(x->x<data[1]*1/2, data)
        x = isnothing(x) ? [length(data)] : [x]
        return  x/frame_rate # 15 fps, return the time in seconds
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

# mean(correlations[1,1], dims=2)[:,1,:]
data = load(joinpath(analysis_path, corr_file)) |> tosymboldict |> dict2ntuple;
config = @strdict data=data functions=functions
fits,_ = produce_or_load(analysis_path, config,filename="$(filename)_fits", force=true) do config
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

data = load(joinpath(analysis_path,"$(filename)_fits.jld2")) |> tosymboldict |> dict2ntuple;
@unpack fits, inverse_fits, times, speeds = data

## Get bandwidth for specific values
fit_params = conf["fit"]
func = fit_params["func"]
coarse_time = fit_params["coarse_time"]
config = (@strdict func coarse_time) |> dict2ntuple

name = savename(joinpath(analysis_path,filename), config, "jld2")
isfile(name) && rm(name)
data,_ = produce_or_load(analysis_path,config,prefix=filename, force=true) do config
    @unpack func, coarse_time = config
    band = zeros(Int, 9,5) # 9 speeds, 5 values, max, min, mean, errmin, errmax
    vals=Vector{Vector{Float64}}(undef, 9)
    diffs=Vector{Vector{Float64}}(undef, 9)
    for s in 1:9
        @info "Speed $(speeds[s])"
        t = findfirst(x->x==coarse_time, times)
        f = findfirst(x->String(x)==func, functions)
        local data = inverse_fits[f,:,t,s]

        means = []
        for mobilmean in 4:10

            zz =runmean(reverse(data),1)
            # make derivative
            max_zz = argmax(zz) 
            my_diff = diff(zz[max_zz+1:end])
            average = runmean(my_diff,mobilmean)[3:end-3]

        # find first positive value
            zz_diff =  begin 
                z = findfirst(x->x>-10e-5,average)
                z = isnothing(z) ? argmax(my_diff) : z
            end
            @info "mobilmean $(mobilmean), zz_diff $(zz_diff+max_zz+2)"

            push!(means,zz_diff+max_zz+1)
        end

        band[s,1] = argmax(reverse(data))
        band[s,3] = minimum(means)
        band[s,4] = maximum(means)
        band[s,2] = round(Int,mean(means))
        band[s,5] = round(Int,std(means))
        @info "Bandwidth $(band[s,:]), err: $(band[s,5])"
        vals[s] = reverse(data)
        # diffs[s] = average
    end
    shearband =(max=band[:,1],min=band[:,2], errmin=band[:,3], errmax=band[:,4], err=band[:,5], values=vals, diffs=diffs)
    @show shearband.min
    return @strdict shearband=shearband func=func coarse_time=coarse_time
end
@unpack shearband = data


# plot!()
# ##
# plot()
# for s in 2:9
#     zz =rollmean(reverse(inverse_fits[2,:,3,s]),4)
#     my_diff = diff(zz[argmax(zz)+1:end])
#     zz_diff =  begin 
#         z = findfirst(x->x>0,my_diff)
#         z = isnothing(z) ? argmax(my_diff) : z
#     end
#     # zz_diff_ = findfirst(x->x>0,my_diff)
#     @show  maximum(diff(zz[argmax(zz)+1:end]))
#     # plot!(my_diff,c=cs[s])
#     plot!(zz,c=cs[s])
#     @show zz_diff
#     scatter!([zz_diff+argmax(zz)+1],[zz[zz_diff+argmax(zz)]],c=cs[s])
# end
# plot!(legend=false)
# ##



fits
times[3]
functions
# label="half")

