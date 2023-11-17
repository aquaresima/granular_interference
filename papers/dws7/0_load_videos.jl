using DrWatson
@quickactivate "granular_media"

using HDF5
using NPZ
import YAML

#%% md
# Load all the matrices. The matrices are exported from the recorded videos and saved to a folder.
#Each matrix has two spatial dimensions (the screen) and a temporal dimension, the first two reflect the pixels of the image, the third corresponds to the video frames.
# Our video has 250x420 pixels and 15 f/s
##
paths = YAML.load_file(projectdir("paths.yml"))
videos = paths["video_npy"]

@info "Load all the npy matrices"
try
    (root, dirs, files) = first(walkdir(videos))
    speed_folders = filter(dir -> startswith(dir, "V"), dirs)
    for dir in speed_folders
        file = joinpath(root, dir, "image_matrix.npy")
        speed = parse(Float64, dir[3:end])
        try
            _data = npzread(file)
            target = joinpath(root, dir * ".h5")
            isfile(target) && rm(target)
            h5open(target, "w") do fid
                fid["matrix"] = _data
                fid["speed"] = speed
            end
            println("loaded:", speed)
        catch
            println("failed:", speed)
        end
    end
catch
    @error "No data found in $videos"
end
