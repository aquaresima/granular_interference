@quickactivate "granular_media"

using RollingFunctions
using speckles
using Plots
using HDF5
using NPZ
import YAML
using RollingFunctions

#
paths = YAML.load_file(projectdir("paths.yml"))
full_matrices  = paths["full_matrices"]

# Analysis of the temporal variance in order to select the best coarse-grain values. (slow ∼ 5 mins)
data, _ = produce_or_load(full_matrices, filename=joinpath(full_matrices,"videos_properties.jld2")) do full_matrices
    file = joinpath(full_matrices,"V_0.00.h5")
    still_video = read(h5open(file))
    @info "Temporal variance analysis"
    vars = speckles.measure_temporal_variance(still_video)
    files = [joinpath(full_matrices,f)  for f in filter(endswith(".h5"),readdir(full_matrices))]
    full = read(h5open(joinpath(files[1]),"r")) |> x->x["matrix"][:,:,1]
    @show size(full)
    norm_matrix =zeros(size(full))
    @info "Compute the normalization factor"
    for file in files
         _fid = read(h5open(joinpath(file),"r"))
        full = mean(_fid["matrix"], dims=3)[:,:,1]
        average = copy(full)
        target_mean = mean(full)
        average[full .< (target_mean -50)]
        height = mean(average, dims=2)[:,1]
        width = mean(average, dims=1)[1,:]
        height_norm = height/ maximum(height)|> x -> 1 ./x |> x-> runmean(x,40) |> x-> x/maximum(x)
        width_norm =  width/maximum(width) |> x-> 1 ./x |> x-> runmean(x,40)    |> x-> x/maximum(x)
        norm_matrix += width_norm' .* ones(size(full)) .* height_norm
    end
    norm_matrix = norm_matrix/maximum(norm_matrix) 
    return @strdict norm=norm_matrix variance = vars
end

