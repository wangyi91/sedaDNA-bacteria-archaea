#!/usr/bin/env julia

### Run after loading data from local-lca.jl ###


# Check only the deepest 5 samples
x = @pipe leftjoin(df, select(mt,[:Label, :Middle_depth, :yrBP, :sigma]), on = :Label) |>
                 filter(:tax_path => t -> occursin("$ftax", t), _) |>
                 filter([:N_reads, :tax_rank] => (x,y) -> x > nreads && y==("$frank"), _) |>
		 sort(_, :Middle_depth) |>
                 dropmissing(_, [:Middle_depth, :D_max])|>
                 filter(:Middle_depth => x -> x > 1000, _)

## summary for taxa that are absent in at least one of the 5 samples
summ = @pipe combine(groupby(x, :tax_name), nrow => "N_samples",
                     :N_reads => sum => :sum_reads,
                     :N_reads => maximum => :max_reads,
                     :N_reads => minimum => :min_reads,
                     :N_reads => median => :med_reads) |>
             sort(_, :max_reads)

sub = @pipe leftjoin(x, select(summ, [:tax_name, :N_samples]), on = :tax_name)

@pipe filter([:Label, :N_samples, :D_max, :lambda_LR] => (x,y,z,w) -> x=="HB_21" && y==1 && z>0.1 && w>5, sub) |> sort(_, :N_reads)
