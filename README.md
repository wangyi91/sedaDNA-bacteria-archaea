#  Prokaryote data analyses of Lake Constance _sed_ aDNA
This repository contains scripts for data analyses and visualisation for bacteria and archaea recorded in the sedimentary ancient DNA of Lake Constance.  

Before using scripts in this repository, raw sequence data (fastq files) should be first processed with scripts in the repository [shotgun-data-processing](https://github.com/wangyi91/shotgun-data-processing.git). 

Then, fastq files are mapped against reference databases using the [aMAW-pipeline] for taxonomic profiling using the [GTDB database](https://gtdb.ecogenomic.org) (this study used the [r207 release](https://data.gtdb.ecogenomic.org/releases/release207/)), and DNA damage is subsequently estimated with [metaDMG](https://github.com/metaDMG-dev/metaDMG-core).

The input file of analyses in this repository is the output file `results/taxonomic-profiling-dmg/tp-mdmg.lca.weight-1.csv.gz` of `metaDMG`. 



