#!/usr/bin/env julia
# Use after loading df, ftax and nreads from local-lca.jl

using CategoricalArrays

libs = load("../metadata/HIBO_library_metadata_dated.csv") |> DataFrame |> dropmissing
tab = DataFrame()

for frank in ["species","genus","family","order","class","phylum"]
    x = @pipe leftjoin(df, select(mt,[:Label, :Middle_depth, :yrBP, :sigma]), on = :Label) |>
          filter(:tax_path => t -> occursin("$ftax", t), _) |>
          filter([:N_reads, :tax_rank] => (x,y) -> x > nreads && y==frank, _) |>
          sort(_, :Middle_depth) |>
          dropmissing(_, [:Middle_depth, :damage])
    summ = combine(groupby(x, [:Label, :Middle_depth, :yrBP, :sigma]),
                                nrow => "taxa_count",
                                :N_reads => sum => :N_reads)
    summ[!, :rank] .= frank
    tab = vcat(tab, summ, cols=:union)
end

transform!(tab, [:Label, :N_reads] => ByRow((l,n) -> n/libs[libs.Label.==l,:].N_filtered_reads[1]) => :N_reads_proportion)

## reorder key_tax for plotting
ranks = CategoricalArray(tab.rank);
levels!(ranks, ["species","genus","family","order","class","phylum"]);

# plot for n of taxa
gr(margins = 2Plots.cm, size=(1000, 700))
plot(tab.yrBP, tab.taxa_count, group=ranks, legend=:outertopright, formatter = :plain, xticks=0:1000:14000);
xlabel!("year BP");
ylabel!("N taxa");
title!("$ftax, $pip");
savefig("../MicrobeProfiling/output/$(pip)_$(ftax)_taxa_at_rank_lineplot.pdf");

# plot for n of reads
gr(margins = 2Plots.cm, size=(1000, 700))
plot(tab.yrBP, tab.N_reads_proportion, group=ranks, legend=:outertopright, formatter = :plain, xticks=0:1000:14000);
xlabel!("year BP");
ylabel!("N reads to total reads in library");
title!("$ftax, $pip");
savefig("../MicrobeProfiling/output/$(pip)_$(ftax)_nreads_at_rank_lineplot.pdf");

