# Brockman Usage
![Kenny B](/images/g4160.png)

The command line tool is a bash pipeline originally aiming to be used with a job scheduler, submitted as an array job. Thus, the script takes as input a table of samples (detailed below, one sample per line), the specific sample line number to do (starting from 1), 

See [installation](brockman.md#installation) for notes on how to install Brockman.

## Usage notes

#### Important notes:
* This is a pipeline designed for sc-ATAC-seq data.  Other data types will require modifications to the pipeline.  Please start an issue for recommended modifications and we will include modified pipelines as required/requested.
* These pipelines are designed to be robust to failure and will only try to redo a step if the previous step failed. Accordingly, if a step fails or needs to be modified, the products of that step must first be deleted before the script will attempt to re-create them.
  * If you need to redo a step, make sure to delete the files produced at this step and all subsequent files that depend on that step (which will not be recreated otherwise).
  * Once the pipeline has completed successfully, it creates a file <tempDir>/<sampleID>.done which must be deleted before the pipeline will attempt to re-run this sample.
* This pipeline is designed to make use of a config file that contains specific analysis parameters.  An example config file is available in the git repo (and is described in more detail in the example below).
* The pipeline is designed to start from fastq files, one per sample.
  * It will require modification if you want it to do more/less.
  
Several resources that may be useful were included in the installation, including Nextera adapter sequences for use with read trimming, and hg19 K562 replication domains.  These can vary depending on how Brockman was installed, but you can find out where they are by doing the following:
```
$ source activate BrockmanEnv
(BrockmanEnv) $ which brockman_pipeline
~/.conda/envs/BrockmanEnv2/bin/brockman_pipeline
```
and so my resource files are here: `~/.conda/envs/BrockmanEnv2/etc/Brockman/Resources/`

## Example

Download this git repo, and `cd` into the Example directory
```bash
git clone https://github.com/Carldeboer/Brockman
cd Example
```

First we will download some example data.  Run the following script (which requires [SRAtoolkit](https://github.com/ncbi/sra-tools) to be installed):
```bash
./downloadExample.sh
```
This script contains the following:
```bash
#!/bin/bash

mkdir -p Downloaded_fastQs
cd Downloaded_fastQs
while read -r id
do
  fastq-dump --split-files --gzip  $id
done < ../SRA_entries.txt
```
and will download all the SRA entries in the file `SRA_entries.txt` and put them in a directory called `Downloaded_fastQs`.

Before we can run the pipeline on this data, we need to modify the config file, which looks like this by default:
```
range=50 # range around the Tn5 insertion site for which DNA is scanned for k-mer content
bowtieIndex="$HOME/genomes/human/hg19/Bowtie2Index" #index to use for alignment
logDir="temp_brockman_files" #where to put temp files
outDir="kmer_frequencies" #where to output results
adapterFile="../Resources/NexteraPE-PE.fa" # Nextera adapter sequences for trimming adaptor sequences from reads
maxInsertionSize=2000 # maximum insertion size for aligning to the genome
chromBED="$HOME/genomes/human/hg19/chroms.bed" # a BED file where each entry is a chromosome with the start and end marking the first and last bases of the chromosome
genome2bit="$HOME/genomes/human/hg19/allChrs.2bit" # a two-bit representation of the chromosome sequences
mitoChrom="chrM" # the mitochondrial chromosome name (excluded from analysis)
G1BED="../Resources/ReplicationDomains/hg19.K562.G1.bed" #BED file of G1 replication domains 
G2SBED="../Resources/ReplicationDomains/hg19.K562.G2S.bed" #BED file of G2 replication domains
sampleFile="sample_list.txt" # file with one sample per line, tab delimited, with first column the (unique) sample ID, the second and third columns are the paths to the corresponding fastq files
```
This file must be in bash syntax and is `source`ed from `brockman_pipeline`, so ensure that there are *no spaces* around the `=`.

Many of these parameters can be left alone, but some must be updated to allow this example to work:
* `bowtieIndex`: This is the path to the `bowtie2` index for the genome you are using.
* `genome2bit`: This is the 2bit chromosome sequence file for the genome, which are generally available from UCSC or can be created using Kent tools. For hg19, this can be found [here](http://hgdownload-test.cse.ucsc.edu/goldenPath/hg19/bigZips/)
* `chromBED`: This is a .bed file containing one entry per chromosome, ranging from 0 to the last base of the chromosome. It can easily be made from the `chrom.sizes` file located [here](http://hgdownload-test.cse.ucsc.edu/goldenPath/hg19/bigZips/). Here is an example:
```
chr1	0	249250621
chr10	0	135534747
chr11	0	135006516
...
```

If you are changing genome version, you will also want to update `G1BED` and `G2SBED` and maybe `mitoChrom` (if the identifier for the mitochondrial chromosome is different). 

Once you have updated the config file, you can run the pipeline on one example, as follows:
```
brockman_pipeline config.sh 1
```

This will run the brockman pipeline on the first (`1`) sample/line of the sampleFile (`sample_list.txt`) specified in `config.sh`

To use a job scheduler, you can run `brockman_pipeline` similar to what follows, although the below example is using GridEngine and this will change according to the job scheduler.
Here is the main job file, called `example_run_pipeine_UGER.sh`:
```bash
#!/bin/bash -l

echo Job $JOB_ID:$SGE_TASK_ID started on $HOST: `date` #print to the output file the current job id and task, and host
use Anaconda3 # or however you enable conda on your system
source activate BrockmanEnv # or whatever the conda Brockman environment was named
source config.sh #load config file (just to get $logDir and $sampleFile)

export id=`awk "NR==$SGE_TASK_ID" $sampleFile` # input the $SGE_TASK_IDth line of this file into $id
splitID=($id) #split id on whitespace into array that can be accessed like ${splitID[0]}
export id=${splitID[0]} # First column of the sample file sheet; the sample ID (must be unique)
export logPre=$logDir/$id
#redirect stdout and stderr to the log file
exec 1>>$logPre.olog
exec 2>&1
../brockman_pipeline config.sh $SGE_TASK_ID

qstat -j $JOB_ID | grep "^usage *$SGE_TASK_ID:" #display job usage stats
```

Which can be run using a command like the following:
```bash
qsub  -b y -cwd  -o run_UGER_pipe.olog -l h_vmem=20g,h_rt=99:00:00 -e run_UGER_pipe.olog -N ..example_run_pipeine_UGER.sh -t 1-10 ./example_run_pipeine_UGER.sh
```

This will run the brockman pipeline on all 10 samples in `sampleFile`.
