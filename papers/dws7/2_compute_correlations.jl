using DrWatson
@quickactivate "granular_media"

using Revise
using speckles
import YAML
using JLD2
using HDF5
using ProgressBars
import speckles: shift_mat, coarsen

conf = YAML.load_file(projectdir("conf.yml"))
# Set the location of the folder and import all files.
paths = conf["paths"]
full_matrices  = paths["full_matrices"]
analysis_path  = conf["do_norm"] ? paths["analysis"] : paths["analysis"] * "_no_norm" |> mkpath
@assert isdir(full_matrices) 

do_norm = conf["do_norm"]
if do_norm 
    norm = load(joinpath(full_matrices,"videos_properties.jld2"))["norm"]
else
    norm = 1.
end

# Set coarse dimensions
coarse = conf["coarse"]
x = coarse["x"]
y = coarse["y"]
frames = conf["frames"]
times = coarse["times"]

# Set filenames
filename ="coarsen$(x)x$(y)"
temporal_file(root,t) = "$(root)_$t"
temporal_file(filename,1)
##


## %% md Spatial coarse grain of each matrix. Metapixels have size 10x10 pixels, from variance analysis (slow ∼ 5 mins)
@info "Coarse grain and compute time-averaged heatmaps"
# data, _ = produce_or_load(full_matrices, filename=coarse_file) do full_matrices

config = @strdict full_matrices 
data,_ = produce_or_load(analysis_path, config,filename=filename) do config
    files = [joinpath(full_matrices,f)  for f in filter(endswith(".h5"),readdir(full_matrices))]
    lk = Threads.ReentrantLock()
    matrices = Vector{Array{Float32,3}}(undef,length(files))
    speeds = Vector{String}(undef,length(files))
    # Threads.@threads  
    for n in eachindex(files)
        file = files[n]
        @info "$file to $( Threads.threadid())"
        _fid = lock(lk) do
            return  read(h5open(joinpath(file),"r"))
        end
        if n>1
            mat = shift_mat(_fid["matrix"])[:,:,1:frames]
        else
            mat = shift_mat(_fid["matrix"])[:,:,1:4500] ## noise only lasts 5min
        end
        @info "$file loaded"
        matrices[n] = Float32.(coarsen(norm .* mat, x=x, y=y))
        if do_norm
            @info "Normalizing"
                matrices[n] = matrices[n] 
        end
        speeds[n] = string(_fid["speed"])
        @info "$file coarsened"
    end
    return @strdict matrices speeds x y
end

data["matrices"]

## Temporal coarsegrain
data = tosymboldict(load(joinpath(analysis_path,filename*".jld2")) |> dict2ntuple;)
@unpack x,y, speeds, matrices= data
for t in ProgressBar(times)
    file = temporal_file(filename,t)
    @show file
    temp_coarse =Vector{Array{Int16,3}}(undef,length(speeds))
    data, _ = produce_or_load(analysis_path, matrices, filename=file) do matrices 
        for k in eachindex(speeds)
            @debug "time: $t, speed: $(speeds[k])"
            m = speckles.smoothen_t(round.(Int16,matrices[k]), t= t)
            temp_coarse[k] = m
        end
        return @strdict matrices=temp_coarse speeds=speeds times=t
    end
end;



##
#====================================
        Correlate signal
====================================#
corr_file = joinpath(analysis_path,"$(filename)_correlations")
data, _ = produce_or_load(times,filename=corr_file) do times
    correlations = Matrix{Array{Float64,3}}(undef,length(times), length(speeds))
    iter = ProgressBar(eachindex(times))
    for t in iter
        file = joinpath(analysis_path,"$(filename)_$(times[t]).jld2")
        set_postfix(iter,Time="$(times[t])")
        temp_matrix =load(file)["matrices"]
        temp_coarse =Vector{Array{Int16,3}}(undef,length(speeds))
        for k in eachindex(speeds)
            @debug  speeds[k], times[t]
            m = speckles.compute_correlation(temp_matrix[k])
            correlations[t,k] =  m
        end
    end
    return @strdict correlations speeds=speeds times=times
end;
##
