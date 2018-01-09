# ZIIFD - An Integrated feature detection algorithm

# Introduction
Zero-Inflated Negative Binomial Feature Detection (ZIIFD) tool is developed for
integration of multiple high-throughput sequencing data types. The tool utilises
an extended previously proposed ZINBA (Zero-inflated negative binomial algorithm)
approach to analyse multiple data tracks simultaneously. This is done by incorporating
a correlation term between the HTS data tracks, so that the algorithm can be run
for two tracks in parallel. To estimate parameters of this model, the iterative
process (EM-type algorithm) is extended from the original by including correlation
estimation based on results of logistic regression and generalised linear model
fitting steps.

The algorithm is implemented with three major steps:
1. Data pre-processing
2. Probabilistic identification of genomic enrichment regions using extended
   mixture regression setup
3. Post-processing of the enriched windows

# Required input files
Before using ZIIFD main function a user has to prepare input files. These files are
generated using code from ZINBA algorithm. In order to build window data from
aligned sequencing reads a user should first download necessary files:

1. Genome build in .2bit format. The files can be downloaded from links:
	- hg19: http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.2bit
	- hg18: http://hgdownload.cse.ucsc.edu/goldenPath/hg18/bigZips/hg18.2bit
2. Mappability files which correspond to sequencing read length. These files are
   used by ZINBA tool in order to generate the final mappability score (check
   https://code.google.com/archive/p/zinba/wikis/UsingZINBA.wiki for details):
	- hg19, map 36bp: http://www.bios.unc.edu/~nur2/map36_hg19.tgz
	- hg19, map 50bp: http://www.bios.unc.edu/~nur2/map50_hg19.tgz
	- hg18, map 36bp: http://www.bios.unc.edu/~nur2/map36.tgz
	- hg18, map 50bp: http://www.bios.unc.edu/~nur2/map50.tgz

# Required R packages
The following R packages must be installed prior to running the program:
- BSgenome
- doMC
- MASS
- matrixStats
- multicore
- R.utils
- zinba

# Implemented Functions
The following functions are provided with the software:

## buildwins
### Description
Generates input files from aligned sequencing reads that are used further by ziifd
function. The files are generated using ZINBA algorithm functions, specifically
buildwindowdata and generateAlignability. Based on the user preferences, function
either generates the mappability files from the raw alignability files or uses the
already existing ones for calculating mappability score (ZINBA in turn uses peakseq
code for generating alignability files).
### Input format
Input files should be provided in either `bed`, `tagAlign` or `bowtie` format.
### Usage
```R
buildwins(pathToZinbaLibLoc = NULL, seq = NULL, inputseq = NULL,
          generateAlign = NULL, inputDirMap = NULL, outdirAlign = NULL,
          inputDirAlign = NULL,  outdirWins =  NULL, twobitfile = NULL,
          extension = NULL, athresh = 1, winSize = 250, offset = 0,
          cnvWinSize = 100000, cnvOffset = 2500, numProc = 4)
```

### Arguments
---
pathToZinbaLibLoc:  path to ZINBA library. Default NULL.
seq:                path to mapped sample reads in bed/tagAlign/bowtie file format
inputseq:           path to mapped input reads in bed/tagAlign/bowtie file format. If
                    left blank, defaults to 'none'.
generateAlign:      boolean value specifying whether the alignability files should be
                    generated or not
inputDirMap:        input directory containing files with unzipped peakseq mappability
                    file (if generateAlign = TRUE)
outdirAlign:        output directory for writing created output alignability files (if
                    generateAlign = TRUE)
inputDirAlign:      input directory with ready alignability files (if generateAlign = FALSE)
outdirWins:         output directory for saving window data
numProc:            number of processing units. Defaults is 4.
twobitfile:         path to build of the genome the reads were mapped to, in .2bit format
extension:          average length of fragments in fragment library used
athresh:            Uniqueness threshold, number of occurrences of a given k-mer imposed
                    during alignment (1 is absolute uniqueness).
winSize:            Window size for build window data, default=250bp
offset:             Bp offset, default is no offsets (offset=0). If one's winSize is 500 and
                    would like 4 equally
would like 4 equally
                    spaced offsets, then specifying offset=125 would achieve this
cnvWinSize:         Size of windows used to calculate CNV activity in sample, default is 100000
cnvOffset:          Offset for CNV windows, typically 2500bp issuffcient, default is no offsets
---

### Output
Function generates files that contain quantified read count, gc content, mappability
score and input count (if available) values for each chromosome. The files are stored
in the user specified output folder.

### Example
This is an example of running the function with input file in BED format.
```bash
# getting mappability files
wget http://compbio.uta.fi/integrated_feature_detection_algorithm/map36_hg19.tgz
```

```R
# generating window data for ChiP-seq track profiling transcription factor occupancy
buildwins(pathToZinbaLibLoc = 'path/to/ZINBA/folder/zinba',
          seq = 'path/to/mapped/reads/in/bed/format/aligned_chip_seq_reads.bed',
          inputseq = 'path/to/mapped/reads/in/bed/format/aligned_control.bed',
          generateAlign = TRUE,
          inputDirMap = 'path/to/zipped/mappability/files/map_dir_in',
          outdirAlign = 'path/to/output/alignability/files/map_dir_out',
          outdirWins =  'path/to/folder/to/store/output/files/dir_out',
          twobitfile = 'path/to/twobitfile/hg19.2bit',
          extension = 200,
          athresh = 1,
          winSize = 250,
          numProc = 4)

# generating window data for DNase-seq/ATAC-seq/FAIRE-seq tracks
buildwins(pathToZinbaLibLoc = 'path/to/ZINBA/folder/zinba',
          seq = 'path/to/mapped/reads/in/bed/format/aligned_seq_reads.bed',
          generateAlign = TRUE,
          inputDirMap = 'path/to/zipped/mappability/files/map_dir_in',
          outdirAlign = 'path/to/output/alignability/files/map_dir_out',
          outdirWins =  'path/to/folder/to/store/output/files/dir_out',
          twobitfile = 'path/to/twobitfile/hg19.2bit',
          extension = 134,
          athresh = 4,
          winSize = 250,
          numProc = 4)
```

## ziifd
### Description
This is the main function of ZIIFD algorithm, which consists of two steps. First,
it runs an iterative extended EM-algorithm for probabilistic classification of
window read counts into one of three components: zero-inflated, background and
enriched given a set of covariates. Windows having posterior probability higher
than a user specified threshold (default parameter of 0.95) are classified as
enriched. Next, adjacent enriched windows are combined to form peaks. If specified,
the algorithm combines an output from two tracks into a consensus track, which is
preferable for analyses of two biological replicates.

Function reads previously generated input text files containing window information
from two input folders corresponding to each track subjected for the analyses. Note
that each chromosome should be summarised in a separate tab-delimited text file. If
a user provides custom input files, they should contain read count information along
with the quantified covariates for each window of fixed size. More specifically, each
file should contain columns with the following attributes: chromosome name, starting
position of window, ending position of window, read count value, input count value,
estimated GC content, mappability score and local background estimation (see example
files for clarification).

For each pair of tracks an extended EM-algorithm (“em” function) is executed, which
provides estimations of posterior probability for each window. Then function “mergepeaks”
is executed, in which the windows are selected according to user specified threshold,
adjacent windows are merged to form the features, and results are saved into the ‘output’
folder. If the consensus track is needed, the output is contained into a single BED file.
Otherwise, the output from each track is saved into a separate BED file.

### Usage
```R
ziifd(pathfilelist1, pathfilelist2, formulaZ1 = NULL, formulaB1 = NULL,
      formulaE1 = NULL, formulaZ2 = NULL, formulaB2 = NULL, formulaE2 = NULL,
      outputPath = NULL, filename1 = NULL, filename2 = NULL, threshold = 0.95,
      consensus = FALSE, numCores = 4)
```

### Arguments
> pathfilelist1/2: 
path to a folder containg a list of files for each chr for the first/second track of data
> formulaZ1/2:
formula for modelling zero-infaleded component, track 1/2
> formulaB1/2:
formula for modelling background component, track 1/2
> formulaE1/2:
formula for modelling enrichment component, track 1/2
> outputPath:
path to output folder
> filename1/2:
prefix used to denote the output files for the firts/second track that are created by ziifd
> consensus:
boolean - specifies whether the consensus track output should be created or not (default FALSE)
> threshold:
threshold of posterior probability for selecting peaks must be between 0 and 1 (default 0.95)
> numCores:
a desirable number of cores to be used for the run (default is 4)

### Output
Function generates BED files containing detected peaks with the calculated
statistics for each peak including posterior probability score characterising
enrichment probability of each peak, and false discovery rate (FDR) value. In
addition to peaks files function creates “wins.txt” files containing full
tables with statistics calculated for each window; this includes posterior
probability and FDR estimations.

### Example
```R
ziifd(pathfilelist1 = 'path/to/second/input/folder/input1',
      pathfilelist2 = 'path/to/second/input/folder/input2',
      formulaZ1 = 'exp_count ~ align_perc + exp_cnvwin_log + gcPerc + 1',
      formulaB1 = 'exp_count ~ align_perc + exp_cnvwin_log + gcPerc + 1',
      formulaE1 = 'exp_count ~ align_perc + exp_cnvwin_log + gcPerc + 1',
      formulaZ2 = 'exp_count ~ align_perc + exp_cnvwin_log + gcPerc + 1',
      formulaB2 = 'exp_count ~ align_perc + exp_cnvwin_log + gcPerc + 1',
      formulaE2 = 'exp_count ~ align_perc + exp_cnvwin_log + gcPerc + 1',
      outputPath = 'path/to/output/folder',
      filename1 = 'FAIREseqRep1',
      filename2 = 'FAIREseqRep2',
      threshold = 0.95,
      consensus = TRUE,
      numCores = 4)

# Example of analysis of FAIRE-seq data and ChiP-seq tracks of the data:
ziifd(pathfilelist1 = 'path/to/second/input/folder/input1',
      pathfilelist2 = 'path/to/second/input/folder/input2',
      formulaZ1 = 'exp_count ~ align_perc + exp_cnvwin_log + gcPerc + 1',
      formulaB1 = 'exp_count ~ align_perc + exp_cnvwin_log + gcPerc + 1',
      formulaE1 = 'exp_count ~ align_perc + exp_cnvwin_log + gcPerc + 1',
      formulaZ2 = 'exp_count ~ align_perc + exp_cnvwin_log + gcPerc + 1',
      formulaB2 = 'exp_count ~ align_perc + exp_cnvwin_log + gcPerc + 1',
      formulaE2 = 'exp_count ~ align_perc + exp_cnvwin_log + gcPerc + 1',
      outputPath = 'path/to/output/folder',
      filename1 = 'FAIREseqRep1',
      filename2 = 'FAIREseqRep2',
      threshold = 0.95,
      consensus = FALSE,
      numCores = 4)
```

