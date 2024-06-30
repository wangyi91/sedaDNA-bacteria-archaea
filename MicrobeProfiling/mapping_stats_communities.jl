#!/usr/bin/env julia

# This script is to plot mapping statistics for taxa of interest 
using DataFrames
using CSVFiles, CSV
using CairoMakie, StatsPlots
using Pipe: @pipe

# load Is, cliq_idx, nsps... from ../NetworkAnalysis/key_communities_plot.jl by running the first 46 lines
#Iall=[[] for _=1:6];
#Iall[1:5] =Is; Iall[6]=I;

push!(titles,"") # update titles
push!(nsps, 0) # update number of species

tb=otudata


# Choose one of the four for plotting
## dmg vs ani
y="read_ani_mean";
y_range = 92:0.05:100 # ani_range
ylab="ANI %"; tit = "ani";

## dmg vs bratio
y="breadth_exp_ratio"
y_range = 0.75:0.001:1 # bratio_range
ylab="Breadth ratio"; tit="bratio";

## dmg vs cov_evenness
y="cov_evenness"
y_range = 0:0.001:0.813 # cov_evenness_range
ylab="Coverage evenness"; tit="covevenness";

## dmg vs coverage_covered_mean
y="coverage_covered_mean";
y_range = 1:0.01:5.3 # coverage_covered_mean_range
ylab="Mean covered depth"; tit="covcovedmean";

using KernelDensity
dmg_range = 0:0.005:0.4
dens = kde((Vector{Float64}(tb.damage), Vector{Float64}(tb[!,y])));
bg = pdf(dens,dmg_range,y_range);

# plot damage vs ANI
using CairoMakie
fig = Figure(resolution=(800,600))
axs = [Axis(fig[m,n], limits=(0,0.3,minimum(y_range),maximum(y_range)), title=permutedims(reshape(titles, (3,2)))[m,n]*"\n$(permutedims(reshape(nsps,(3,2)))[m,n]) species") for m in 1:2, n in 1:3]
ylabel = Label(fig[1:2, 0], ylab, rotation = pi/2)
xlabel = Label(fig[3,1:3], "DNA damage")

for m in 1:2, n in 1:3
    i=m*m+n-1; if i==6 break end
    clqs = collect.(cliper[cliq_idx])[Is[i]]
    clq = reduce(vcat, clqs) |> sort |> unique
    splist = collect(skipmissing(tax_at_rank.(names(netw)[clq], "species")))
    sub = @pipe otudata |> filter(:tax_name => n-> n in splist, _)
    # subset mapping statistics to taxa of interest
    tab = @pipe tb |> filter(:reference=>r->r in unique(sub.reference), _)
    # plot background (all taxa)
    CairoMakie.contourf!(axs[m,n], dmg_range, y_range, bg, levels = 10, colormap=:Greys_4)
    # plot species in community-group
    CairoMakie.scatter!(axs[m,n], Vector{Float64}(tab.damage), tab[!,y], 
                        color=tab.Middle_depth./maximum(tb.Middle_depth), colormap=cgrad(:starrynight, rev=true), colorrange = (0, 1), 
                        markersize=5, alpha=0.7, strokewidth=0.1, strokecolor=(:black,0.3))
end
save("./MicrobeProfiling/output/dmg_vs_$(tit)_communities.pdf", fig);




