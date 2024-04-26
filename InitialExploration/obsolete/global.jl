#!/usr/bin/env julia

using DataFrames
using FileIO
using CSVFiles
using GZip
using StatsPlots

# load global summary
gl = load(File(format"CSV", "../data/tp-mdmg.global.weight-1.csv.gz")) |> DataFrame

# load library metadata and join to df
mt = load("../metadata/HIBO_library_metadata.csv") |> DataFrame

# join
# create a df with only data that need to be plotted
df = DataFrame([gl.sample gl.D_max gl.lambda_LR], [:sample, :D_max, :lambda_LR]);

rename!(df, 1 => "Label")
df = innerjoin(df, select(mt, [:Label, :Top_depth]), on = :Label)


# plot
gr(margins = 1.5Plots.cm)
scatter(df.D_max, df.lambda_LR,group=ismissing.(df.Top_depth), markercolor=:match, labels=["libaries" "controls"], legend=:outertopright)
xlabel!("D_max global")
ylabel!("lambda_LR global")
title!("global")
savefig("./output/global_scatter.pdf")

