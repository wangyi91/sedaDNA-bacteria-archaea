#  Prokaryote data analyses of Lake Constance _sed_ aDNA
## Background

This repository contains scripts for data analyses and visualisation for bacteria and archaea recorded in the sedimentary ancient DNA of Lake Constance.  

The manuscript of this study is currently under review in _Current Biology_.

## Prerequisite steps

Before using scripts in this repository, raw sequence data (fastq files) should be first processed with scripts in the repository [shotgun-data-processing](https://github.com/wangyi91/shotgun-data-processing.git). 

Then, fastq files are mapped against reference databases using the [aMAW-pipeline] for taxonomic profiling using the [GTDB database](https://gtdb.ecogenomic.org) (this study used the [r207 release](https://data.gtdb.ecogenomic.org/releases/release207/)), and DNA damage is subsequently estimated with [metaDMG](https://github.com/metaDMG-dev/metaDMG-core).

## Input files
The input files for analyses in this repository are: 

* The output file `results/taxonomic-profiling-dmg/tp-mdmg.lca.weight-1.csv.gz` of `metaDMG`. 
* Metadata are located in [metadata](https://github.com/wangyi91/sedaDNA-bacteria-archaea/tree/main/metadata)

## Steps of sequence analyses
[InitialExploration]
[DamageAnalysis]
[NetworkAnalysis]
[MicrobeProfilling]

## Age-depth model
[OxCal_model]
[AgeDepthModel]
