#!/usr/bin/env julia

### Run after loading data "x" from local-lca.jl ###

# plot D vs. L vs. age
h=90;l=90;
gr(margins = 1.5Plots.cm, size=(500, 625))
@pipe dp |> plot(_.damage, _.significance, _.yrBP, marker_z=_.yrBP,
     color=:matter, seriestype=:scatter, mode="markers",
     #marker_size=log.(x.N_reads),
     marker=(log10.(_.N_reads),0.3, stroke(0)), marker_sizeref=0.3,
     label="", camera = (h,l));

xlabel!("Damage");
ylabel!("Significance");
zlabel!("Age (yr BP)");
savefig("../InitialExploration/output/age-damage_DvsL_$(tag)_$(ftax)_$(frank)_$(nreads)_camera$h$l.pdf");




# plot D vs. age, filtered data, specific taxon
dp = @pipe data_filtered |> filter(:tax_name => t -> t=="g__Fagus",_) |> select(_, [:damage, :yrBP, :N_reads])
gr(margins = 1.5Plots.cm, size=(500, 625))
@pipe dp |> plot(_.yrBP, _.damage,
     seriestype=:scatter, mode="markers",
     #marker_size=log.(x.N_reads),
     marker=(log10.(_.N_reads),0.3, stroke(0)), marker_sizeref=0.3,
     label="");

xlabel!("Age (yr BP)");
ylabel!("Damage");
savefig("../InitialExploration/output/age-damage_DvsAge_$(pip)_$(ftax)_Fagus_filtered.pdf");





# plot N vs. L vs. age

h=90;l=90;
gr(margins = 1.5Plots.cm, size=(500, 625))
@pipe dp |> plot(log1p.(_.N_x1_forward), _.significance, _.yrBP,
     marker_z=_.yrBP,
     color=:matter, seriestype=:scatter, mode="markers",
     marker=(log1p.(_.N_x1_forward),0.3, stroke(0)), marker_sizeref=0.3,
     label="", camera = (h,l));

xlabel!("log N_x1_forward");
ylabel!("Significance");
zlabel!("Age (yr BP)");
savefig("../InitialExploration/output/age-damage_NvsL_$(tag)_$(ftax)_$(frank)_$(nreads)_camera$h$l.pdf");


 ## data to plot
 dp=x
 lib="HB_21"
 dp = @pipe x |> filter(:Label => x -> x==lib,_) |> filter(:N_reads => x -> x>=nreads, _) |>
 transform(_, [:phi, :A]  => ByRow((x,y) -> (x.*y.-1)/(x.-2)) => :mode) |>
 filter(:mode => x -> x>=0, _)

# plot phi vs A vs LR
h=40;l=10;
gr(margins = 1.5Plots.cm, size=(500, 625))
@pipe dp |> plot(_.phi, _.A, _.significance, marker_z=_.significance, 
                 marker=(log1p.(_.N_x1_forward) ,0.3, stroke(0)), marker_sizeref=0.3,
                 color=:matter, seriestype=:scatter, mode="markers",
                 label="", camera = (h,l));

xlabel!("phi");
ylabel!("A");
zlabel!("significance");
savefig("../InitialExploration/output/age-damage_phiASig_$(lib)_$(tag)_$(ftax)_$(frank)_$(nreads)_camera$h$l.pdf");

#plot mode vs sig
gr(margins = 1.5Plots.cm, size=(500, 625))
@pipe dp |> plot(_.mode, _.significance, 
                 marker=(1, 0.3, stroke(0)), 
                 #marker_z=_.yrBP, color=:matter,
                 seriestype=:scatter, mode="markers",
                 label="");
xlabel!("mode");
ylabel!("significance");
savefig("../InitialExploration/output/age-damage_modesig_$(lib)_$(tag)_$(ftax)_$(frank)_$(nreads).pdf");



# plot mode against age as boxplot
gr(margins = 1.5Plots.cm, size=(500, 625))
@pipe dp |> boxplot(_.yrBP, _.mode, label="", bar_width=40,
        orientation=:h, yflip = true,
        line=(:black, 1, 0.5), marker=(:black, 1, 0.3, stroke(0)));
ylabel!("Age");
xlabel!("Mode of k(x)~(μ(x),Φ) (x=1)");

savefig("../InitialExploration/output/age-damage_modeAge_$(lib)_$(tag)_$(ftax)_$(frank)_$(nreads)_boxplot.pdf");



# histogram
gr(margins = 1.5Plots.cm, size=(500, 625))
@pipe dp |> histogram(_.mode, weights=sqrt.(max.(_.N_x1_forward, 1)));
xlabel!("mode");
savefig("../InitialExploration/output/histmode_weightedNx1_$(lib)_$(tag)_$(ftax)_$(frank)_$(nreads).pdf");


