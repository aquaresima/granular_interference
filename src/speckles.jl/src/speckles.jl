module speckles
    using HDF5
    using LinearAlgebra
    using StatsBase
    using Statistics
    using RollingFunctions
    using ProgressBars
    using Logging
	include("functions.jl")

    export Logging, ProgressBars, RollingFunctions, HDF5

end # module speckles