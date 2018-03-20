## Brockman: Brockman Representation Of Chromatin by K-mers in Mark-Associated Nucleotides
![Kenny B](/images/g4160.png)

# Overview
## What is Brockman?
Brockman is a suite of command line tools and R functions to convert genomics data into DNA k-mer words representing the regions associated with a chromatin mark, and then analyzing these k-mer sets to see how samples differ from each other.
This approach is primarily intended for single cell genomics data, and was tested most extensively on single cell ATAC-seq data. The bash scripts in particular may require some alteration for other types of genomics data.

## What are Broackman's dependencies?
The command line tools rely on the following, and assume the shell is Bash:
* [Ruby](https://www.ruby-lang.org/en/)
* [AMUSED](https://github.com/Carldeboer/AMUSED): for counting k-mers
* [BEDTools](http://bedtools.readthedocs.io/en/latest/): For working with BED files
* [Kent Tools - twoBitToFa](http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/): For extracting genomic sequence
* [SAM Tools](http://samtools.sourceforge.net/): For working with BAM/SAM files

The R analysis tools rely on the following packages:
* [jackstraw](https://cran.r-project.org/web/packages/jackstraw/jackstraw.pdf)
* [ggplot2](https://cran.r-project.org/web/packages/ggplot2/ggplot2.pdf)
* [tsne](https://cran.r-project.org/web/packages/tsne/tsne.pdf)

## Installation
TBD
### Command line tools


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

### R library
See [Brockman Examples](brockman_examples.md) for example analysis pipelines.
