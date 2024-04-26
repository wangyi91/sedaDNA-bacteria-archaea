#!/usr/bin/env julia

using Arrow
using DataFrames
using GZip
using CSVFiles
using Statistics
using Pipe: @pipe
using StatsPlots

include("./InitialExploration/_initial_load.jl")
include("./MicrobeProfiling/_get_annotated_df.jl") # to use the tax_at_rank function

pip="holi"
pip="amaw"

# use when remap using reads filtered by $tax
tax=".Viridiplantae"
tax=""

target=".cpnr"
target=".norwaycpnr"
target=".HB_31.dedup.bam"
target=""

alg = ".local"
alg = ".lca";

add=""; # this defult is ANI90
add="_ANI95";
add="_ANI92";

tag = "$pip$tax$target$alg$add"
tag="holi.lca";












## count number of reads for each taxon through all samples
summ = @pipe combine(groupby(x, :tax_name), nrow => "N_samples", 
                     :N_reads => sum => :sum_reads,
                     :N_reads => maximum => :max_reads,
                     :N_reads => minimum => :min_reads,
                     :N_reads => median => :med_reads) |> 
             sort(_, :max_reads) |>
             filter(:N_samples => n -> n>1, _)

## count number of taxa for each sample
summ = @pipe x|> DataFrames.transform(_, :tax_path => ByRow(p->tax_at_rank(p,"superkingdom")) => :superkingdom) |> 
                 DataFrames.transform(_, :tax_path => ByRow(p->tax_at_rank(p,"kingdom")) => :kingdom) |>
                 DataFrames.transform(_, :tax_path => ByRow(p->tax_at_rank(p,"phylum")) => :phylum) |>
                 DataFrames.transform(_, :tax_path => ByRow(p->tax_at_rank(p,"class")) => :class) |>
                 DataFrames.transform(_, :tax_path => ByRow(p->tax_at_rank(p,"family")) => :family) |>
                 DataFrames.transform(_, :tax_path => ByRow(p->tax_at_rank(p,"genus")) => :genus) |>
                 filter(:superkingdom => t -> t == "Eukaryota",_) |>
                 filter(:class => t -> !ismissing(t),_) |>
                 filter(:phylum => t -> !ismissing(t),_) |>
                 filter(:phylum => t -> t != "Streptophyta",_) |>
                 combine(groupby(_, [:Label,:phylum, :Middle_depth, :yrBP, :sigma]), nrow => "count", :N_reads => sum => :N_reads) |>
                 sort(_, :N_reads)

for ph in unique(summ.phylum) 
dp = @pipe filter(:phylum => p-> p==ph,summ) |> sort(_,:yrBP)
gr(margins = 1.5Plots.cm, size=(1000, 600))
plot(dp.yrBP, dp.count)
plot!(dp.yrBP, dp.count, seriestype=:scatter)

xlabel!("year BP")
ylabel!("count of genus in $ph")
savefig("./InitialExploration/output/summary_$(ph)_taxcount.pdf")
end

## data of controls
c = @pipe leftjoin(df, select(mt,[:Label, :Middle_depth, :yrBP, :sigma]), on = :Label) |>
          filter(:tax_path => t -> occursin("$ftax", t), _) |>
          filter([:N_reads, :tax_rank] => (x,y) -> x > 10 && y==("$frank"), _) |>
          filter(:Middle_depth => x -> ismissing(x), _)


 ### plot 				
gr(margins = 1.5Plots.cm, size=(500, 625))
bar(summ.Middle_depth, summ[!,5], 
	orientation=:h, label="", yflip=true, 
	color=:olivedrab, line=(:black, 1, 0.2)) # thistle, dodgerblue
ylabel!("Depth (cm)")
xlabel!("Count of $ftax $frank, reads>$nreads")

savefig("./output/summary_$(tag)_$(ftax)_$(frank).pdf")


      

# plot damage of samples
gr(margins = 1.5Plots.cm, size=(500, 625))

boxplot(x.Middle_depth, x.damage, label="", bar_width=20, 
        orientation=:h, yflip = true,
        line=(:black, 1, 0.5), marker=(:black, 1, 0.3, stroke(0)))
        
title!("$tag ($ftax, $frank, N_reads > $nreads)")
ylabel!("Depth (cm)")
xlabel!("damage")

savefig("./output/$(tag)_$(ftax)_$(frank)_$(nreads)_boxplot.pdf")

        
# plot damage of controls
gr(margins = 1.5Plots.cm)
boxplot(c.Label, c.damage, label="",
        line=(:black, 1, 0.5), marker=(:black, 1, 0.3, stroke(0)))

title!("$tag ($ftax, $frank, N_reads > $nreads)")
xlabel!("Control name")
ylabel!("damage")

savefig("./output/control_$(tag)_$(ftax)_$(frank)_$(nreads)_boxplot.pdf")

