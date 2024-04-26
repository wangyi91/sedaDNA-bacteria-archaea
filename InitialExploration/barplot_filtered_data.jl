#!/usr/bin/env julia
# use after loading "data_filtered" from vcat_samples.jl
using StatsPlots
using DataFrames
using JLD2
using Pipe
using StatsBase
using CategoricalArrays

function tax_at_rank(path::String, rank::String)
        regex = Regex("(?<=:)([^\"]+)(?=(:\"$rank))")
        m1 = match(regex, path).match
        if occursin("__", m1)
            regex = Regex("(?<=__)([^\"]+)(?=(:\"$rank))")
        end
    return(match(regex, path).match)
end

function clean_tax(tax::String)
    if occursin("__", tax)
        out=match(r"(?<=__)(.*)", tax).match
    else out=tax
    end
    return(out)
end

ftax = "Viridiplantae"
plant = @pipe data_filtered |> filter(:tax_path => t -> occursin("$ftax", t), _)
x = plant

ftax = "Chordata"
vert = @pipe data_filtered |> filter(:tax_path => t -> occursin("$ftax", t), _)
x=vert

ftax = "Chlorophyta"
phy = @pipe data_filtered |> filter(:tax_path => t -> occursin("$ftax", t), _)
x=phy

ftax = "Cyanobacteria"
cya = @pipe data_filtered |> filter(:tax_path => t -> occursin("$ftax", t), _)
x=cya


## count number of reads for each taxon through all samples
summ = @pipe combine(groupby(x, :tax_name), nrow => "N_samples",
                     :N_reads => sum => :sum_reads,
                     :N_reads => maximum => :max_reads,
                     :N_reads => minimum => :min_reads,
                     :N_reads => median => :med_reads) |>
             sort(_, :max_reads) |>
             filter(:N_samples => n -> n>1, _)

key_tax = filter([:max_reads, :min_reads] => (m,n) -> m>=200 && n>=50, summ)[!,:tax_name] 
# or 50, 10 for bacteria
# or 40, 10 for chordata

# get info for key_tax
check = @pipe x |> filter(:tax_path => t -> tax_at_rank.(t, "genus") in clean_tax.(key_tax), _) |>
transform(_, :tax_path => ByRow(t -> tax_at_rank.(t, "genus")) => :genus)

ngr=size(unique(check.genus),1)

## reorder key_tax for plotting
ctg = CategoricalArray(check.genus)
        # for chordata
        levels!(ctg, ["Danio", "Coregonus", "Salmo",  "Perca", "Taurulus", "Esox", "Oncorhynchus", "Gadus", "Salvelinus", "Pimephales", "Cyprinus", "Barbus", "Petromyzon",
              "Rana", "Bufo",
              "Gallus",
              "Talpa", "Myodes", "Mus", "Arvicola", "Castor", "Microtus",
              "Sus", "Bos", "Cervus",
              "Homo"]) 

## color scheme
pal = pip=="amaw" ? :Accent_8 : :Set3_8

gr(margins = 5Plots.cm, size=(2000, 1000))
@pipe check |> combine(groupby(_, [:Label, :genus, :Middle_depth, :yrBP, :sigma]), 
                       nrow => "count_genera",
                      :N_reads => sum => :sum_reads) |> 
               groupedbar(Float64.(_.yrBP), sqrt.(_.sum_reads), group=ctg,
                          bar_position = :stack, legend=:outertopright, color=reshape(palette(pal, ngr)[1:ngr],(1,ngr)),
                          bar_width=180, formatter = :plain, xticks=0:1000:14000);
xlabel!("year BP");
ylabel!("√N reads");
title!("$ftax, $pip");
savefig("../InitialExploration/output/$(ftax)_$(pip)_keytaxa_barplot.pdf");


