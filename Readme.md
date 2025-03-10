# **A Phylogeny-aware algorithm for Quantifying Microbiome Transitions Across Scales and Habitats（MT-Tracker）**

## Introduction

Microbial transition tracker (MT-Tracker) is an easy-to-use bioinformatic package designed to analyze microbial transition patterns. ME-Tracker enables the qualitative and quantitative exploration of microbial community transition. MT-Tracker works as a plug-in tool for Parallel-META 3 and supports both 16S rRNA amplicon sequences and whole-genome sequencing (WGS) data as input.

## Software Requirement and Dependency

#### Software requirements

OpenMP library is the C/C++ parallel computing library. Most Linux releases have OpenMP already been installed in the system. In Mac OS X, to install the compiler that supports OpenMP, we recommend using the Homebrew package manager:

```
brew install gcc
```

## Installation Guide

#### Automatic installation (recommended)

**a. Download the package**

```
git clone https://github.com/qdu-bioinfo/MT-tracker.git
```

**b. Install by installer**

```
cd MT-tracker
source install.sh
```

The package should take less than 1 minute to install on a computer with the specifications recommended above.

#### Manual installation

If the automatic installer fails, MT-tracker can still be installed manually.

**a. Download the package**

```
git clone https://github.com/qdu-bioinfo/MT-tracker.git
```

**b. Configure the environment variables (the default environment variable configuration file is “~/.bashrc”)**

```
export MTTRACKER=Path to MT-tracker
export PATH="$PATH:$MTTRACKER/bin"
source ~/.bashrc
```

**c. Compile the source code**

```
cd MT-tracker
make
```

## Basic Usage

With a input 16S rRNA amplicon sequence file, e.g. sample1.fasta:

**a. Profiling by **

[Parallel-Meta Suite]: https://github.com/qdu-bioinfo/parallel-meta-suite

** **

```
PM-parallel-meta -r sample1.fasta -o sample1.out
```

The “sample1.out” folder is the profiling result.

**b. Calculate the virtual ancestor and transition direction between two samples.**

```
Mt-tracker -i sample1.out/classification1.txt sample1.out/classification2.txt
```

## Batch Processing

MT-Tracker also supports the batch input of profiling results by the following alternative two forms (compatible with Parallel-META 3):

**a. Sample list**

```
Mt-tracker -l samples.list -o samples.mtt
```

in which parameter “-l” assigns the file list of profiling results of multiple samples. The format of a sample list:

```
Sample1	/home/data/sample1.out/classification.txt
Sample2	/home/data/sample2.out/classification.txt
...	
SampleN	/home/data/sampleN.out/classification.txt
```

**b. Abundance tables**

```
Mt-tracker -T samples.OTU.Abd -o samples.mtt
```

in which parameter “-T” assigns the profiling result of OTU table of multiple samples. The format of a OTU table:

```
SampleID	Sample1	Sample2	Sample3	
OTU1	100	200	0	50
OTU2	0	300	600	100
OTU3	50	80	0	200
```

## Calculate the transition direction and probability between multiple samples

For multiple samples, after calibration and normalization, the transition probabilities among multiple samples are obtained. Input the samples.mtt file and the corresponding meta file.

```
Mt-tracker -g samples.mtt samples.meta -o output.txt
```

## Example Dataset

Here we provide a demo dataset with profiling results of 40 microbiome and corresponding meta file:

```
ls example
```

or type the following command:

```
Mt-tracker -T samples.OTU.Abd -o samples.mtt
Mt-tracker -g samples.mtt samples.meta -o output.txt
```

For other options, use Mt-tracker -h to view