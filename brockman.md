## Brockman: Brockman Representation Of Chromatin by K-mers in Mark-Associated Nucleotides
![Kenny B](/images/g4160.png)

# Overview
## What is Brockman?
Brockman is a suite of command line tools and R functions to convert genomics data into DNA k-mer words representing the regions associated with a chromatin mark, and then analyzing these k-mer sets to see how samples differ from each other.
This approach is primarily intended for single cell genomics data, and was tested most extensively on single cell ATAC-seq data. The bash scripts in particular may require some alteration for other types of genomics data.

A preprint describing the approach is available [here](https://www.biorxiv.org/content/early/2018/04/03/129247).

## What are Brockman's dependencies?
The command line tools rely on the following, and assume the shell is Bash:
* [Ruby](https://www.ruby-lang.org/en/)
* [AMUSED](https://github.com/Carldeboer/AMUSED): for counting k-mers
* [BEDTools](http://bedtools.readthedocs.io/en/latest/): For working with BED files
* [Kent Tools - twoBitToFa](http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/): For extracting genomic sequence
* [SAM Tools](http://samtools.sourceforge.net/): For working with BAM/SAM files
* [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic): For trimming sequencing reads

The R analysis tools rely on the following packages:
* [jackstraw](https://cran.r-project.org/web/packages/jackstraw/jackstraw.pdf)
* [ggplot2](https://cran.r-project.org/web/packages/ggplot2/ggplot2.pdf)
* [tsne](https://cran.r-project.org/web/packages/tsne/tsne.pdf)

## Installation

### Command line tools

At present, only anaconda installation is supported. If you haven't yet learned how to use anaconda, there's no time like the present!
#### Linux/OSX:
```bash
conda create -c bioconda -n BrockmanEnv  brockman-pipeline
```

### R library

If you don't already have `devtools`, install it:
```
install.packages("devtools")
```

Load `devtools` and install from the GitHub page:

```
library(devtools)
install_github("Carldeboer/BrockmanR")
```


## Usage

### Command line tools
See [Brockman_pipeline Example](brockman_pipe_example.md) for example data processing pipelines.


### R library
See [Brockman Analysis Example](brockman_example.md) for example analysis pipelines.

## Citation

Please cite the [Brockman preprint](https://www.biorxiv.org/content/early/2018/04/03/129247) if you find Brockman useful.

Carl de Boer, Aviv Regev. *BROCKMAN: Deciphering variance in epigenomic regulators by k-mer factorization*. bioRxiv 129247; doi: https://doi.org/10.1101/129247
