# NanobodyMiner

This git repo contains a Julia pipeline for processing Alpaca nanobody phage libraries. This code is under development, and does not yet have a license.

## Installation
Clone the repo. Go to the directory you wish to install to, and:
```
git clone git@github.com:MurrellGroup/NanobodyMiner.git
```
This will give you a directory with:
```
Demo/
Panels/
src/
README.md
```
The easiest way to use this is to duplicate (and rename) the Demo directory for for each analysis. Demo contains:
```
Fastq/
config.yaml
special.csv
NanobodyMinerPipeline.ipynb
```
Fastq/ is where the merged .fastq data files live. There should be one .fastq file for each post-panning run, and one for the "baseline", which as a special status.

config.yaml is how you control the run, by specifying which datasets get included, which primers to use, which enrichment conditions should be used to select promising candidates, what the minimum count threshold should be, etc.

special.csv contains a set of sequences that you might already have information from in this dataset. These (and their neighbours) will be excluded from the automated selection, and they will be annotated in the plots.

NanobodyMinerPipeline.ipynb is the pipeline notebook itself, in which all the processing code will run.

## Required Packages
Before this will run, the following meed to be installed, either through your Julia REPL, or at the top of the NanobodyMinerPipeline.ipynb notebook.

```
using Pkg
Pkg.add(["YAML", "PyPlot", "DataFrames", "CSV", "NearestNeighbors", "LightGraphs", "LazySets", "Distances", "SparseArrays"])
Pkg.add(PackageSpec(;name="SeqUMAP",url="https://github.com/murrellb/SeqUMAP.jl.git"))
Pkg.add(PackageSpec(name="NextGenSeqUtils", rev="1.5.3", url = "https://github.com/MurrellGroup/NextGenSeqUtils.jl.git"))
```

## Running
Once you have duplicated the Demo directory, placed the .fastq files in the Fastq/ directory, and pointed to them correctly in the config.yaml file, just open up the notebook and "Restart & Run All".

The main outputs will be in the Plots/ directory, and the automated selection will be in Selections/ - but don't just trust these!! Check them.
