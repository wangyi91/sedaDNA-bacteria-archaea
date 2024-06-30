#!/usr/bin/env julia

using Statistics, CSVFiles, DataFrames, Pipe, StatsBase
using Chron

# load model specification from OxCal
spec = @pipe load("./data/model_specification.csv") |> DataFrame |> dropmissing(_, :z)
sort!(spec,:z)

# Age interpolation at given depths
mt = load("../metadata/HIBO_library_metadata.csv") |> DataFrame

# initialise data table for age interpolation
interp = @pipe dropmissing(mt, :Top_depth_shifted) |>
                   transform(_, :Top_depth_shifted => ByRow(x->x.+0.5) => :Middle_depth) |>
                   select(_, :Label, :Middle_depth) |>
                   sort(_, :Middle_depth) |>
                   allowmissing(_)
                  
interp.yrAD = linterp1s(spec.z, spec.mu, interp.Middle_depth)
interp.sigma = linterp1s(spec.z, spec.sigma, interp.Middle_depth)

dated = @pipe transform(interp, [:yrAD, :sigma] => ByRow((a,b) -> replace([a,b], NaN => missing)) => [:yrAD, :sigma]) |>
        transform(_, :yrAD => ByRow(x -> ismissing(x) ? missing : 1950-x) => :yrBP) |>
        leftjoin(mt, _, on = :Label) |>
        select(_, Not([:Top_depth, :Index_PCR_cycle])) |>
        sort(_, :Middle_depth)

# set estimation at depths beyond dated range to missing
#transform!(interp, [:Middle_depth, :Age, :Age_025CI, :Age_975CI] => 
#                        ByRow((a,b,c,d) -> extrema(smpl.Height)[1] < -a < extrema(smpl.Height)[2] ? [a,b,c,d] : [a,missing,missing,missing]) 
#                        => [:Middle_depth, :Age, :Age_025CI, :Age_975CI])

save("./output/age_interpolation.csv", interp, quotechar=nothing)
save("../metadata/HIBO_library_metadata_dated.csv", dated, quotechar=nothing)

