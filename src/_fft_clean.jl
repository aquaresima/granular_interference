
## Remove external signals with FFT not used for the present analysis
#
# # an example
# coarse10_file = joinpath(analysis_path,"coarsen10.h5")
# matrix = read(h5open(coarse10_file,"r"))["1.0"]
# corr = Speckles.compute_correlation(matrix)
# v = StatsBase.mean(corr[:,:,1:300], dims=2)[:,1,:]
# p1 = plot(v', c=cs', xlims=(1,250), legend=false, title="Raw")
# p3 = plot(v', c=cs', xlims=(1,50), legend=false, title="Raw zoom")
#
# clean = clean_recordings(matrix)
# corr = Speckles.compute_correlation(clean)
# v = StatsBase.mean(corr[:,:,1:300], dims=2)[:,1,:]
# p2 = plot(v', c=cs', xlims=(1,250), legend=false, title="Cleaned")
# p4 = plot(v', c=cs', xlims=(1,50), legend=false, title="Cleaned zoom")
# plot(p1,p2,p3,p4)

# ##
# # heatmap(mean(matrix, dims=3)[:,:,1])
#
#
# freq_domain = plot(freqs, abs.(full_f), title = "Spectrum", xlim=(-1000, +1000))
# # freq_domain = scatter!(freqs[maxima], abs.(full_f)[maxima], title = "Spectrum", xlim=(-1000, +1000))
# freq_domain = scatter!(freqs[alls], abs.(full_f[alls]), title = "Spectrum", xlim=(-1000, +1000))
# # freq_domain = scatter!(freqs[isolated_2], abs.(full_f[isolated_2]), title = "Spectrum", xlim=(-1000, +1000))
# # freq_domain = scatter(freqs[isolated], abs.(full_f[isolated]), title = "Spectrum")
#
# # isolated




# full_f[isolated] .= 0
# freq_domain = plot(freqs, abs.(full_f), title = "Spectrum", xlim=(-1000, +1000))
# plot(abs.(ifft(full_f)))
# # plots
# time_domain = plot(t, signal, title = "Signal")
# freq_domain = plot(freqs, abs.(full_f), title = "Spectrum", xlim=(-1000, +1000))
# plot(time_domain, freq_domain, layout = 2)
# savefig("Wave.pdf")

##
coarse10_file = joinpath(analysis_path,"coarsen10.h5")
all_matrix = read(h5open(coarse10_file,"r"))
coarse10_clean_file = joinpath(analysis_path,"coarsen10_clean.h5")
isfile(coarse10_clean_file) && rm(coarse10_clean_file)
@printf "Remove steady noise"

fid = h5open(coarse10_clean_file,"w")
for speed in speeds
    @show speed
    # fid[speed] = clean_recordings(all_matrix[speed])
    fid[speed] = all_matrix[speed]
end
close(fid)

h5open()

## Vertical cross-correlation
#
# corr_mats = Dict(t=>[] for t in times)
# corr_file = joinpath(analysis_path,"coarsen10_vertical_correlation.h5")
# isfile(corr_file) && rm(corr_file)
# fid =h5open(corr_file, "w")
#
# for t in times
#     temporal_file = joinpath(analysis_path,"coarsen10_time$t.h5")
#     all_matrix =read(h5open(temporal_file, "r"))
#     fid_t = g_create(fid, t)
#     for k in speeds
#         @show t, k
#         m = Speckles.correlate_vertically(all_matrix[k])
#         push!(corr_mats[t], m)
#         fid_t[k] = m
#     end
# end
# close(fid)

h5open(fid)
