#!/usr/bin/env julia
using JLD2
#using CSVFiles
using DataFrames
using Pipe
using CSV


#libs = load("./metadata/HIBO_library_metadata_dated.csv") |> DataFrame |> dropmissing
libs = CSV.File("./metadata/HIBO_library_metadata_dated.csv"; missingstring="NA") |> DataFrame |> dropmissing;
#ftax = ""; pip = "holi"
#ftax = "Viridiplantae"; pip = "amaw"
ftax = "Bacteria"; pip = "amaw"; frank = "species"; add="";add2="";add3="_holocene";add3="_top408cm";add="_ANI95";add2="_strict";
add3="_filtered"; add3="_dropped"; add3="_all";
#ftax = "Archaea"; pip = "amaw"; frank = "species"; add="";add2="";

# get kept ancient taxa
data_filtered = DataFrame();
for lib in libs.Label
    df = load_object("./DmgMixtureModel/output$(add)/obs_$(lib)_$(pip)_$(ftax)_$(frank)$(add)$(add2)_filtered.jld2") 
    data_filtered = vcat(data_filtered, df, cols=:union)
end

# or, get dropped non-ancient taxa
data_dropped = DataFrame();
for lib in libs.Label
    df = load_object("./DmgMixtureModel/output$(add)/obs_$(lib)_$(pip)_$(ftax)_$(frank)$(add)$(add2)_dropped.jld2")
    data_dropped = vcat(data_dropped, df, cols=:union)
end




# or get data for taxa above certain read count
thres_train=500
data_training = DataFrame();
for lib in libs.Label
    df = @pipe load_object("./DmgMixtureModel/data$(add)/obs_$(lib)_$(pip)_$(frank).jld2") |>
    filter(:tax_path => t -> occursin("$ftax", t), _) |>
    filter([:N_reads] => x -> x >=thres_train, _)
    data_training = vcat(data_training, df, cols=:union)
end

# if using amaw local where tax_path is not avaiable && filtering based on dmg is not necessary for direct downstream analysis (network):
data_training = x



