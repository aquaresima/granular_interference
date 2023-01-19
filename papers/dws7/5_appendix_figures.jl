## Verify that the image has an exponential distribution in light-intensity to verify the measure regard actual spekles. : Goodman Statistical Optics
speeds
p = [histogram(m[:], norm=true , alpha=1., bins=-0:5:254, title="ω = $s (s⁻¹)") 
    for (s,m) in collect(zip(speeds,matrices))[[1,3,6,9]]]
plot!(p[3], xlabel="                                               Pixel intensity (255 bits)")
plot(p..., layout=(2,2), legend=false, size=(800,600))
intensity = plot!( legend=false)
savefig(intensity,plotsdir("Fig_A1_intensity.pdf"))


## Correlate stretch fit of raw data with gaussian on averaged
plots = []
PALETTE = 9
cs = cgrad(:roma, 1:PALETTE)[collect(1:PALETTE) ./PALETTE]
p1 = plot()

coarse = "1"
for (n,s) in enumerate(speeds[2:end])
    p1 = scatter!(p1,[x.stretched[1:1] for x in times_fits["1"][s]],    [x.gaussian for x in times_fits[coarse][s]], label=false, c=cs[n], title="", msc=cs[n], ms=6)
end
p = plot!(p1, ylabel="Gaussian Exp. - time avg", xlabel="Stretched Exp. - raw data")
#
savefig(p, joinpath(results_path,"figA1.pdf"))
p

##
PALETTE = 9
cs = cgrad(:roma, 1:PALETTE)[collect(1:PALETTE) ./PALETTE]
time_plots = []

z = ones(1,9)
for x in 1:9
    z[x] = x
end
# for t in times
t = "20"
g_p = plot(legend=false, title="Gaussian")
for (s,c) in zip(speeds[2:end],cs)
    data = times_fits[t][s]
    for (n,d) in enumerate(data)
        scatter!(g_p, [42-n],d.gaussian , c=c, ylabel="Gaussian Exp. τ (ms)", xticks=(0:10:40, 0:4:16),  msc=c, ms=6 )
    end
end
#     push!(time_plots,g_p)
# end
# plot(time_plots...)

# fit_plots = plot(time_plots[5], xlabel="z (mm)")
plot!(g_p, xlabel="z (mm)")
heatmap!(z[:,2:end], c=cs[2:end], colorbar=false, title="ω (s⁻¹)", yticks=:none, titlefontsize=13,
    inset_subplots = bbox(0.15, 0.80, 0.3, 0.07, :bottom), subplot=2, axes=false,
	 xticks=(2:2:9, speeds[2:2:9]), xrotation=-30)
savefig(g_p, joinpath(results_path,"figA2.pdf"))
plot!(g_p)
##