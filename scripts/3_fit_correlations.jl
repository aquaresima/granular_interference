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
full_matrices = conf["paths"]["full_matrices"]
analysis_path =
    conf["do_norm"] ? conf["paths"]["analysis"] : conf["paths"]["analysis"] * "_no_norm"
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
temporal_file(root, t) = "$(root)_$t.jld2"
temporal_file(filename, 10)

## Fit exponentials
functions = [:stretched, :gaussian, :cosh, :half]
N_points = 10


τ_window = 5
@. gauss_exp(x, p) = exp(-(2x / sqrt(p[1] * τ_window))^2)
@. stretch_exp(x, p) = exp(-(x / abs(p[1]))^abs(p[2]))
@. fit_cosh(x, p) = 2 - cosh(x / p[1])
function halfheight(data)
    x = findfirst(x -> x < data[1] * 1 / 2, data)
    x = isnothing(x) ? [length(data)] : [x]
    return x / frame_rate # 15 fps, return the time in seconds
end
function get_fits(data)
    global N_points
    tt = N_points
    return (
        stretched=curve_fit(stretch_exp, 3tt:4tt, data[tt:2tt], [tt, 1.0]).param,
        gaussian=curve_fit(gauss_exp, 1:tt, data[1:tt], [1.0]).param,
        cosh=curve_fit(fit_cosh, 1:tt, data[1:tt], [1.0]).param,
        half=halfheight(data),
    )
end

corr_file
data = load(joinpath(analysis_path, corr_file)) |> tosymboldict |> dict2ntuple;
config = @strdict data = data functions = functions
fits, _ = produce_or_load(
    analysis_path,
    config,
    filename="$(filename)_fits",
    force=true,
) do config
    data = config["data"]
    functions = config["functions"]
    @unpack correlations, speeds, times = data
    ys = size(correlations[1, 1], 1)
    xs = size(correlations[1, 1], 2)
    inverse_fits =
        Array{Float64,4}(undef, size(functions, 1), ys, size(times, 1), size(speeds, 1))
    fits =
        Array{Float64,4}(undef, size(functions, 1), ys, size(times, 1), size(speeds, 1))
    fits_all = Array{Float64,5}(
        undef,
        size(functions, 1),
        ys,
        xs,
        size(times, 1),
        size(speeds, 1),
    )
    for s in eachindex(speeds)
        for t in eachindex(times)
            mysample = mean(correlations[t, s], dims=2)[:, 1, :]
            mysample_all = correlations[t, s]
            @info "Run fit $(speeds[s]) $(times[t]), matrix: $(size(mysample))"
            for y = 1:size(mysample)[1] # number of rows
                _fits = @views get_fits(mysample[y, :])
                inverse_fits[:, y, t, s] =
                    [1.0 / getfield(_fits, f)[1] for f in functions]
                fits[:, y, t, s] = [getfield(_fits, f)[1] for f in functions]
            end
            for y = 1:size(mysample_all)[1] # number of rows
                for x = 1:size(mysample_all)[2] # number of rows
                    _fits = @views get_fits(mysample_all[y, x, :])
                    inverse_fits[:, y, t, s] =
                        [1.0 / getfield(_fits, f)[1] for f in functions]
                    fits_all[:, y, x, t, s] = [getfield(_fits, f)[1] for f in functions]
                end
            end
        end
    end
    return @strdict fits inverse_fits times speeds fits_all
end

data = load(joinpath(analysis_path, "$(filename)_fits.jld2")) |> tosymboldict |> dict2ntuple;
@unpack fits, fits_all, inverse_fits, times, speeds = data


##
# using Plots
# plot([heatmap(1 ./fits_all[2,:,:,5,x], c=:amp, cbar=false) for x in 1:9]..., plot_title="Inverse_Tau(x,y)")
# plot([heatmap(fits_all[2,:,:,5,x], c=:amp, cbar=false) for x in 1:9]..., plot_title="Tau(x,y)")
# plot(inverse_fits[2,:,:,2])
# plot()    )
# my_data = 1./fits_all[[2,4],:,:,5,:]
# NPZ.save(my_data, "here.npz")
##

## Get bandwidth for specific values
fit_params = conf["fit"]
func = fit_params["func"]
coarse_time = fit_params["coarse_time"]
config = (@strdict func coarse_time) |> dict2ntuple

name = savename(joinpath(analysis_path, filename), config, "jld2")
isfile(name) && rm(name)
data, _ = produce_or_load(analysis_path, config, prefix=filename, force=true) do config
    @unpack func, coarse_time = config
    band = zeros(Int, 8, 5) # 9 speeds, 5 values, max, min, mean, errmin, errmax
    vals = Vector{Vector{Float64}}(undef, 8)
    diffs = Vector{Vector{Float64}}(undef, 8)
    for s = 1:8
        @info "Speed $(speeds[s])"
        t = findfirst(x -> x == coarse_time, times)
        f = findfirst(x -> String(x) == func, functions)
        local data = inverse_fits[f, :, t, s+1]

        means = []
        for mobilmean = 3:13

            zz = runmean(reverse(data), mobilmean)
            # make derivative
            max_zz = argmax(zz)
            my_diff = diff(zz[max_zz+1:end])
            average = runmean(my_diff, mobilmean)[3:end-3]

            # find first positive value
            zz_diff = begin
                z = findfirst(x -> x > -10e-6, average)
                z = isnothing(z) ? argmax(my_diff) : z
            end
            @info "mobilmean $(mobilmean), zz_diff $(zz_diff+max_zz+2)"

            push!(means, zz_diff + max_zz + 1)
        end

        band[s, 1] = argmax(reverse(data))
        band[s, 3] = minimum(means)
        band[s, 4] = maximum(means)
        band[s, 2] = round(Int, mean(means))
        band[s, 5] = round(Int, std(means))
        @info "Bandwidth $(band[s,:]), err: $(band[s,5])"
        # vals[s] = my_diff
        # diffs[s] = average
    end
    shearband = (
        max=band[:, 1],
        min=band[:, 2],
        errmin=band[:, 3],
        errmax=band[:, 4],
        err=band[:, 5],
        values=vals,
        diffs=diffs,
    )
    @show shearband.min
    return @strdict shearband = shearband func = func coarse_time = coarse_time #vals=vals
end

@unpack shearband = data

shearband.err
