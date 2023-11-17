"""
Coarse matrices by reducing the dimensions on the x, y or t axes.
"""
function coarsen_xy(matrix::Array{Int16,3}; x::Int64, y::Int64)
    _x, _y, _t = floor.(Int,size(matrix) ./ [x,y, 1])
    _matrix = zeros(_x,_y,_t)
    @inbounds @fastmath for n in 1:_x
        for m in 1:_y
            _matrix[n,m, :] .= mean(matrix[1+(n -1)*x:n*x, 1+(m-1)*y:m*y, :], dims=(1,2))[1,1,:]
        end
    end
    return _matrix
end
function coarsen_xy(matrix::Array{Float64,3}; x::Int64, y::Int64)
    _x, _y, _t = floor.(Int,size(matrix) ./ [x,y, 1])
    _matrix = zeros(_x,_y,_t)
    @inbounds @fastmath for n in 1:_x
        for m in 1:_y
            _matrix[n,m, :] .= mean(matrix[1+(n -1)*x:n*x, 1+(m-1)*y:m*y, :], dims=(1,2))[1,1,:]
        end
    end
    return _matrix
end

function coarsen_xy(matrix::Array{Float64,2}; x::Int64, y::Int64)
    _x, _y = floor.(Int,size(matrix) ./ [x,y])
    _matrix = zeros(_x,_y)
    @inbounds @fastmath for n in 1:_x
        for m in 1:_y
            _matrix[n,m] = mean(matrix[1+(n -1)*x:n*x, 1+(m-1)*y:m*y])
        end
    end
    return _matrix
end

## Run a chunk-average function that normalize the signal
function coarsen_t(matrix::Array{Int16,3}; t::Int)
    _x, _y, _t = floor.(Int,size(matrix) ./ [1.,1., t])
    _matrix = zeros(_x,_y,_t)
    @inbounds @fastmath for n in 1:_t
        _matrix[:,:,n] .= mean(matrix[:, :, 1+(n-1)*t:n*t], dims=3)[:,:,1]
    end
    return _matrix
end

## Run a smoothing function that remove the noise
function smoothen_t(matrix::Array{Int16,3}; t::Int)
    _x, _y, _t = size(matrix)
    z = rollmean(matrix[1,1,:], t)
    _matrix = zeros(_x,_y,size(z)...)
    @inbounds @fastmath for n in 1:_x
        for m in 1:_y
            z = rollmean(matrix[n,m,:], t)
            try
                _matrix[n,m,:] .= z
            catch
                @warn size(z), size(_matrix[n,m,:])
            end
        end
    end
    return round.(Int16, _matrix)
end

function coarsen(matrix::Array{Int16,3}; x=nothing, y=nothing, t=nothing)
    isnothing(t) && @time return coarsen_xy(matrix, x=x, y=y)
    isnothing(x) && @time return coarsen_t(matrix, t=t)
end
function coarsen(matrix::Array{Float64,3}; x=nothing, y=nothing, t=nothing)
    isnothing(t) && @time return coarsen_xy(matrix, x=x, y=y)
    isnothing(x) && @time return coarsen_t(matrix, t=t)
end
function coarsen(matrix::Array{Float64,2}; x=nothing, y=nothing, t=nothing)
    isnothing(t) && @time return coarsen_xy(matrix, x=x, y=y)
    isnothing(x) && @time return coarsen_t(matrix, t=t)
end

"""
Measure the stasndard deviation of the matrix over the timein a progressive coarsen of the x and y axis.
The returned matrix suggests the best coarsening to be applied in order to reduce the temporal noise.
"""
function time_variance(matrix::Array{Int8,3})
    _x, _y, _ = floor.(Int,size(matrix) ./3)
    variances = zeros(_x,_y)
    Threads.@threads for n in 3:_x
    end
    return variances
end
function compute_variance(matrix::Array{Float64,3})
end


function compute_correlation(matrix, lags=300)
    lags = minimum([lags, size(matrix)[3]-1])
    corr_mat = zeros(size(matrix)[1:2]...,lags)
    for m in 1:size(corr_mat)[2]
        for n in 1:size(corr_mat)[1]
            corr_mat[n,m,:] .= autocor(matrix[n,m,:],0:lags-1, demean=true)
        end
    end
    return corr_mat
end

function correlate_vertically(matrix, lags=300; max_depth=30)
    lags = minimum([lags, size(matrix)[3]-1])
    n_cols = size(matrix)[2]
    cross_depth = size(matrix)[1] - max_depth
    corr_mat = zeros(max_depth, cross_depth+1, lags)
    Threads.@threads for h in 1:max_depth
        for n in 0:cross_depth
            _average = zeros(lags)
            for col in 1:n_cols
                _average += crosscor(matrix[h,col,:], matrix[h+n,col,:], 0:lags-1, demean=true)
            end
            corr_mat[h,1+n,:] .= _average/n_cols
        end
    end
    return corr_mat
end

"""
Return the temporal variance of each (meta)pixel.
Compute coarse-grained matrix for metapixel sizes 1->25 and compute the variance.
"""
function measure_temporal_variance(data)
    matrix = round.(Int16,(data["matrix"]))
    vars = zeros(25)
    Threads.@threads for x in 1:25
        m = coarsen(matrix, x=x,y=x)
        v = mean(std(m, dims=3)) 
        vars[x] = v
    end
    return vars
end

function measure_snr(data)
    matrix = round.(Int16,(data["matrix"]))
    I_mean = zeros(25)
    I_each = zeros(25)
    Threads.@threads for x in 1:25
        m = coarsen(matrix, x=x,y=x)
        I = [mean(m[i,j,:]) for i in axes(m,1) for j in axes(m,2)]
        σ = [std(m[i,j,:]) for i in axes(m,1) for j in axes(m,2)]
        I_mean[x] = mean(I)/mean(σ)
        I_each[x] = mean(I ./ σ)
    end
    return I_mean, I_each
end


## Clean recordings
function get_spectrum(signal)
    # Number of points
    signal = signal .- mean(signal)
    N = length(signal)
    # Sample period
    Ts = 1 / (N)
    # Start time
    t0 = 0
    tmax = N * Ts
    # time coordinate
    t = t0+Ts:Ts:tmax
    length(t)
    # Fourier Transform of it
    F = fft(signal) |> fftshift
    # F[argmax(abs.(F)) ] = 0.
    freqs = fftfreq(length(t), 1.0/Ts) |> fftshift
    return F, freqs
end


function findlocalmaxima(signal::Vector)
    inds = Int[]
    if length(signal)>1
        if signal[1]>signal[2]
            push!(inds,1)
        end
        for i=2:length(signal)-1
            if signal[i-1]<signal[i]>signal[i+1]
                push!(inds,i)
            end
        end
        if signal[end]>signal[end-1]
            push!(inds,length(signal))
        end
    end
    inds
    end


function find_isolated_maxima(signal::Vector, coef=40)
    maxima = findlocalmaxima(signal)
    isolated =Vector{Int}()
    for n in 2:length(maxima)-1
        low = coef*signal[rand(maxima)]
        upper = coef*signal[rand(maxima)]
        # upper = 2mean(signal[maxima[n]+1: maxima[n]+4])
        if signal[maxima[n]] > low
            if signal[maxima[n]] > upper
                push!(isolated, maxima[n])
            end
        end
    end
    return isolated
end
function get_maxima(signal)
    alls = Vector{Int64}()
    isolated = [1]
    while !isempty(isolated)
        isolated = find_isolated_maxima(signal)
        alls = vcat(alls, isolated)
        signal[isolated] .= 0
    end
    return alls
end

function clean_recordings(matrix)
    signal = matrix[1,1,:]
    full_f = zeros(length(signal))
    clean_matrix = zeros(size(matrix)...)
    freqs = 0
    for x in 1:42
        for y in 1:25
            signal = matrix[x,y,:]
            f, freqs = get_spectrum(signal)
            full_f += f
            global freqs = freqs
        end
    end
    full_f ./(42*25)

    alls = get_maxima(abs.(full_f))
    @show length(alls)

    for x in 1:42
        for y in 1:25
            signal = matrix[x,y,:]
            f, freqs = get_spectrum(signal)
            f[alls] .= 0
            clean_matrix[x,y,:] = abs.(ifft(f))
        end
    end
    return clean_matrix
end

function shift_mat(mat, add=false)
	a = Array{Int16}(mat)
    if add
        a[a.<0] .= abs.(a[a.<0]) .+127
    else
        a[a.<0] .= a[a.<0] .+ 254
    end
	return a
end


export coarsen