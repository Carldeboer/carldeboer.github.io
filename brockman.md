## Brockman: Brockman Representation Of Chromatin by K-mers in Mark-Associated Nucleotides
![Kenny B](/images/g4160.png)

# Overview
## What is Brockman?
Brockman is a suite of command line tools and R functions to convert genomics data into DNA k-mer words representing the regions associated with a chromatin mark, and then analyzing these k-mer sets to see how samples differ from each other.
This approach is primarily intended for single cell genomics data, and was tested most extensively on single cell ATAC-seq data. The bash scripts in particular may require some alteration for other types of genomics data.

## What are Brockman's dependencies?
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

### Command line tools

At present, only anaconda installation is supported. If you haven't yet learned how to use anaconda, there's no time like the present!
#### Linux:
```bash
conda create -c bioconda -n BrockmanEnv  ruby samtools bedtools ucsc-twobittofa bowtie2 amused trimmomatic brockman
```
#### OSX:
OSX already has ruby installed and including it in the conda environment appears to break ruby due to some missing libraries
```bash
conda create -c bioconda -n BrockmanEnv  samtools bedtools ucsc-twobittofa bowtie2 amused trimmomatic brockman
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
The command line tool is a bash pipeline originally aiming to be used with a job scheduler, submitted as an array job. Thus, the script takes as input a table of samples (detailed below, one sample per line), the specific sample line number to do (starting from 1), 

#### Important notes:
* This is a pipeline designed for sc-ATAC-seq data.  Other data types will require modifications to the pipeline.  Please start an issue for recommended modifications and we will include modified pipelines as required/requested.
* These pipelines are designed to be robust to failure and will only try to redo a step if the previous step failed. Accordingly, if a step fails or needs to be modified, the products of that step must first be deleted before the script will attempt to re-create them.
  * If you need to redo a step, make sure to delete the files produced at this step and all subsequent files that depend on that step (which will not be recreated otherwise).
  * Once the pipeline has completed successfully, it creates a file <tempDir>/<sampleID>.done which must be deleted before the pipeline will attempt to re-run this sample.


### R library
See [Brockman Analysis Example](brockman_example.md) for example analysis pipelines.
