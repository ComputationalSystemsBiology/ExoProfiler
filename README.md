
ExoProfiler
========================================================================

>
> Analyzing transcription factor (TF) binding sites using ChIP-Seq and ChIP-exo data. 
>

1. Requirements
1. Pipeline presentation
1. Example analysis
1. Licence

***

Requirements
------------------------------------------------------------------------

### MatrixScanWS

Developed with:

* Python (2.7.5).

Needed packages: 

* suds (0.4.1), numpy

Other imports:

* logging, os, os.path, platform, urllib2, argparse, warnings, re


### 5PrimeCounter

Developed with:

* Python (2.7.3).

Needed packages:

* numpy, [HTSeq](https://pypi.python.org/pypi/HTSeq) (see [installation guide](http://www-huber.embl.de/users/anders/HTSeq/doc/install.html)), [pysam](https://pypi.python.org/pypi/pysam) (HTSeq dependency)
* (optional) pyfasta, FASTA package (need depend on command line options)


### ExoPlotter

Developed with:

* R (2.14.1).

Needed packages:

* gdata (2.13.3, dependency: gtools will also be installed).
* ape 

***
***


Pipeline presentation
------------------------------------------------------------------------

The full pipeline is composed of three different steps.

***

### 1. MatrixScanWS (Optional step)

MatrixScanWS takes advantage of RSAT web services in order to convert a BED file into FASTA sequences, which are then scanned using a given matrix (transfac format).

#### Steps and RSAT tools involved (as web services) : 

1. Conversion of the input bed file (option -b or --bed\_file) into a fasta file, 
using a reference genome (option -g or --genome). This step uses ```fetch_sequences```.
2. Scan the FASTA sequences for motifs, provided as a transfac matrix (option -m
or --matrix\_file). Only matches with a p-value below a threshold are reported (option
-u or --uth\_pval). This step uses ```matrix_scan```.
3. Matches are written in an output file (option -o or --output_file).

#### Creation of a validation set using permuted matrices:

If the --perm (or -p) option is included in the command, additional steps will 
be performed.
The input transfac matrix will be permuted 10 times using the tool ```convert_matrix```. Those 
permuted matrices should not be too similar to the input transfac matrix. This is checked 
using the tool ```compare_matrices```. They should not be duplicated either. This is 
checked again using ```compare_matrices```.
The permuted matrices will then be used by ```matrix_scan``` to scan the FASTA sequences 
(computed instep 1/), using the same upper threshold for the p-value as in step 2/.
Finally, permuted matrices and outputs from ```matrix_scan``` are written in the same 
folder as the output file of the input transfac matrix.  

#### Aditional information

Tool developed in Python by Celine Hernandez (version 0.4).
  
***

  
### 2. 5PrimeCounter

  
5PrimeCounter analyses a BAM file from a ChIP-Exo experiment in the context of potential transcription factor binding sites (TFBS) showing a given sequence motif of interest.  
  
Please note that in case of motif hits on both strands in the same location, only the one with highest score will be considered.  
  
As input files, 5PrimeCounter receives a BAM file (and its BAI index) and a set of sequence motifs as created by MatrixScanWS, for instance, or any output file of RSAT's tool ```matrix_scan ```. Alternatively, the sequence motifs file can be replaced by a BED file. TOCHECK  
  
  
  
#### Use cases
  
  
  
  
##### Use case 1 : basic usage

Necessary options :

 Short name | Long name | Description |
 --- | --- | --- |
 -bam \<BAM_FILE> | \--bam_file \<BAM_FILE> | Bam file of the ChIP-exo experiment. Requires the index of the bam file with name \<BAM_FILE>.bai in the same folder. |
 -i \<INPUT_SITES> | \--input_sites \<INPUT_SITES> | Input is a file containing predicted Transcription Factor Binding Sites in RSAT matrix-scan output format. It can also accept BED files: see '--input-format' option below for further details. |
 -o \<OUTPUT_PREFIX> | \--output_prefix \<OUTPUT_PREFIX> | Output file prefix. All output files will have that prefix in pathname. |

Default window size around motifs is set to 60 bases. This parameter can be modified using the -s option.


Short name | Long name | Description
--- | --- | ---
-s \<SIZE> | --size \<SIZE> | Window size around motif for which the profile will be computed.
-ob | \--output_bed | Write a BED file with the binding site regions defined by --size and/or --order_by_score.


Additionally, user can provide BED file instead of a Matrix-scan file. This makes it necessary to tell so to 5PrimeCounter by replacing by adding '-if bed' or '--input_format bed' to the command line.


Short name | Long name | Description
--- | --- | ---
-if {matrix-scan,bed} | \--input_format {matrix-scan,bed} | Format of the file containing predicted Transcription Factor Binding Sites. Defaults to 'matrix-scan'. Can also accept 'bed' as value.

  
  
  
##### Use case 2 : scan plus validation set generation


Short name | Long name | Description
--- | --- | ---
 &nbsp; | \--perm \<PERM_FOLDER> | Compute profiles for permuted matrices. ExoProfile analyses all files with the .tf extension.

  
  
  
##### Use case 3 : using input genome

Short name | Long name | Description
--- | --- | ---
-g \<GENOME_SEQ> | \--genome_seq \<GENOME_SEQ> | Reference genome sequence as FASTA file. If this optional argument is given, the consesus sequence is plotted at the bottom of profile and heatmap plots and the sequences fo all binding will be written to integer encoded matrix files. Note, while reading a fasta file the first time, an index file is build in the same directory for faster access. First execution can thus be slower.
-of | \--output_seq | Write a genomic sequences in FASTA format for the binding site regions defined by --size and --order_by_score.

  
  
  
##### Use case 4 : motif comparison (overlay mode)

A shift distance, manually computed, can be used to shift sites by a given base pair number.

Short name | Long name | Description
--- | --- | ---
-sd \<SHIFT_DIST> | \--shift_dist \<SHIFT_DIST> | Shift sites by given distance (in bp) to the right (if positive) or to the left (if negative).

  
  
  
#### Aditional information

##### Other general options

 * Change order of output motifs

Short name | Long name | Description
--- | --- | ---
 &nbsp; | &nbsp; | By default, output file motifs are sorted by occupancy level (number of total read counts).
-os | \--order_by_score  | Sort output files motif by score instead of occupancy level.

 * Change number of output motifs

Short name | Long name | Description
--- | --- | ---
-n \<NUMBER_OF_SITES> | \--number_of_sites \<NUMBER_OF_SITES> | The number of sites to be considered. For a given N take only the top N sites by occupancy level (or motif score if -os is set)
-p \<PERCENT_OF_SITES> | \--percent_of_sites \<PERCENT_OF_SITES> | The percent of sites to be considered. For a given P take only the top P percent sites by occupancy level (or motif score if -os is set)

 * Miscellaneous other options

Short name | Long name | Description
--- | --- | ---
-h | \--help | Show this help message and exit
-fs | \--flip_strand  |  Flipt the strand of motif matches from '+' to '-' and from '-' to '+' for all input sites.


##### Dependencies

5PrimeCounter depends on python packages 'numpy', 'pysam' and 'HTSeq'. If a reference genome is provided to calculate consensus sequences (see Use case 3), the program also imports the 'pyfasta' package.

Tool developed in Python by Jonas Ibn-Salem (version 0.1).

  
***

  
### 3. ExoPlotter ###

ExoPlotter generates a set of PDF plots using output files generated by 5PrimeCounter. This R script needs as input the prefix provided to 5PrimeCounter.

Tool developed in R by Jonas Ibn-Salem (version 0.1), with code from Morgane Thomas-Chollier and Samuel Collombet.


***
***


Example analysis
------------------------------------------------------------------------

***

### MatrixScanWS (Optional)

Using web services, scan a set of sequences for a motif provided as a matrix.

#### Input files

BED file : ```ExoProfiler/data/inputs_small/GR_IMR90.chip-seq_peaks.bed.center60```
Output of a ChIP-Seq experiment from a human sample, i.e. list of peaks.

Matrix (transfac) : ```ExoProfiler/data/inputs_small/MA0113.2.tf```
Matrix describing a motif as available in [JASPAR database](http://jaspar.genereg.net/cgi-bin/jaspar_db.pl?ID=MA0113.1&rm=present&collection=CORE).

#### Use cases

Command with provided example data, assuming location in folder `ExoProfiler/00_matrixScanWS`.

##### Use case 1: basic usage.

    mkdir outputs_example
    python matrixScanWS.py --bed_file ../example_data/matrixScanWS/GR_IMR90.chip-seq_peaks.bed.center60 --genome hg19 --matrix_file ../example_data/matrixScanWS/MA0113.2.tf --uth_pval 0.001 --output_file ./outputs_example/output_matrix_scan_NR3C1_0-001.txt

Output file

* `matrix_scan` output : ```./outputs_example/output_matrix_scan_NR3C1.txt```

##### Use case 2: scan plus validation set generation

Same command with --perm option.

    mkdir outputs_example
    python matrixScanWS.py --bed_file ../example_data/matrixScanWS/GR_IMR90.chip-seq_peaks.bed.center60 --genome hg19 --matrix_file ../example_data/matrixScanWS/MA0113.2.tf --uth_pval 0.001 --output_file ./outputs_example/output_matrix_scan_NR3C1_0-001.txt --perm

Output files

* `matrix_scan` output : ```./outputs_example/output_matrix_scan_NR3C1.txt```
* 10 permuted versions of input matrix, in ```./outputs_example/output_example/``` : ```indiv_motif.MA0113.2_perm1.tf```, ```indiv_motif.MA0113.2_perm2.tf```, ```indiv_motif.MA0113.2_perm3.tf```, etc.
* 10 `matrix_scan` outputs computed from these permuted versions of input matrix, in ```./outputs_example/output_example/``` : ```output_matrix_scan_NR3C1_perm1.txt```, ```output_matrix_scan_NR3C1_perm2.txt```, ```output_matrix_scan_NR3C1_perm3.txt```, etc.



***


### 5PrimeCounter

5PrimeCounter performs a profile-based analysis by calculating a 5' coverage profile, given an BAM file from a ChIP-exo experiment (accompanied by its index in BAI) and a list of potential TF binding sites (TFBS) from motif matching analysis.


#### Use cases

Command with provided example data, assuming location in folder `ExoProfiler/01_5PrimeCounter`.

##### Use case 1 : basic usage

#### Input files

* BAM file : ```ExoProfiler/data/inputs_small/IMR90.rep2.uniq_mapped.chr20_sorted.bam```
Sequence read alignment data in binary format.

* BAI file : ```ExoProfiler/data/inputs_small/IMR90.rep2.uniq_mapped.chr20_sorted.bam.bai```
Index of the BAM file, created with `samtools index`. This index must be located in the same folder as the BAM file.

* ```matrix_scan``` output : ```ExoProfiler/data/inputs_small/output_matrix_scan_NR3C1.txt```
Output from RSAT tool ```matrix_scan```, or directly generated by MatrixScanWS.


#### Command

Basic usage of 5PrimeCounter implies providing at least three parameters:

 * --input_sites for a matrix in matrix-scan format
 * --bam_file for a BAM file
 * --output_prefix to provide folder name and a prefix for output files.

Example of a simple command line :

    python ./01_5primeCounter/5PrimeCounter.py --input_sites ./data/inputs_small/output_matrix_scan_NR3C1_chr20.txt --bam_file ./data/inputs_small/IMR90.rep2.uniq_mapped.chr20_sorted.bam --output_prefix ./output_5PrimeCounter

By default, 5PrimeCounter will assume that the 'input sites' file is an output from `matrix_scan`, so there is no need to specify option `--input_format matrix-scan`. Alternatively, one can also run 5PrimeCounter with a BED file as input, but input format has to be specified. 

    python 5PrimeCounter.py --input_format bed --input_sites ../example_data/input_file.bed --bam_file ../example_data/IMR90.rep2.uniq_mapped.chr20_sorted.bam --output_prefix ./output_5PrimeCounter

Example of a case where the size of the window around the motif has been changed with option --size. Default window size is set to 60.

    python ./01_5primeCounter/5PrimeCounter.py --input_sites ./data/inputs_small/output_matrix_scan_NR3C1_chr20.txt --bam_file ./data/inputs_small/IMR90.rep2.uniq_mapped.chr20_sorted.bam --output_prefix ./output_5PrimeCounter_s100 --size 100

Example of generation of an additional BED output, in case we changed the size of our region of interest.

    python ./01_5primeCounter/5PrimeCounter.py --input_sites ./data/inputs_small/output_matrix_scan_NR3C1_chr20.txt --bam_file ./data/inputs_small/IMR90.rep2.uniq_mapped.chr20_sorted.bam --output_prefix ./output_5PrimeCounter_s100obed --size 100 --output_bed


* Output files

 File name | Content
--- | ---
\<output_prefix>.consensus.txt | Consensus motif.
\<output_prefix>.down_counts.tab | Up counts.
\<output_prefix>.up_counts.tab | Down counts.

 
##### Use case 2 : validation using permuted matrices 

Following Use Case 2 of MatrixScanWS a set of permuted motifs can be generated. One can then run 5PrimeCounter on all the generated matrix-scan outputs.

    python ./01_5primeCounter/5PrimeCounter.py --input_sites ./data/inputs_small/output_matrix_scan_NR3C1_chr20.txt --bam_file ./data/inputs_small/IMR90.rep2.uniq_mapped.chr20_sorted.bam --output_prefix ./output_5PrimeCounter_wperm --perm

* Output files

Same files as for Use Case 1, times ten, for each permuted motif.


##### Use case 3 : using input genome

Genome files in FASTA format can be quite voluminous. Default behaviour doesn't need an input genome. But if provided, additional information can be computed, and another quality check plot can be drawn with ExoPlotter. Output files in BED or FASTA can also be written.

#### Input files

* BAM file : ```ExoProfiler/data/inputs_small/IMR90.rep2.uniq_mapped.chr20_sorted.bam```
Sequence read alignment data in binary format.

* BAI file : ```ExoProfiler/data/inputs_small/IMR90.rep2.uniq_mapped.chr20_sorted.bam.bai```
Index of the BAM file, created with `samtools index`. This index must be located in the same folder as the BAM file.

* matrix\_scan output : ```ExoProfiler/data/inputs_small/output_matrix_scan_NR3C1_chr20.txt```
Output from RSAT tool ```matrix_scan```, or directly generated by MatrixScanWS.

* Reference genome (optional, see Use case 4): a FASTA file as downloaded from, for instance, UCSC's FTP site : ftp://hgdownload.cse.ucsc.edu/goldenPath/ Example of a command: `curl -o h19_chr20.fa.gz ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/chr20.fa.gz` Please note that files need to be decompressed into FASTA in order to be read by 5PrimeCounter.

#### Command

    python ./01_5primeCounter/5PrimeCounter.py --input_sites ./data/inputs_small/output_matrix_scan_NR3C1_chr20.txt --bam_file ./data/inputs_small/IMR90.rep2.uniq_mapped.chr20_sorted.bam --output_prefix ./output_5PrimeCounter_wgenomechr20 --genome_seq ./data/inputs_small/hg19_chr20.fa

NB: By default, 5PrimeCounter will assume that the 'input sites' file is an output from `matrix_scan`, so there is no need to specify option `--input_format matrix-scan`.

Output files

 File name | Content
--- | ---
\<output_prefix>.consensus.txt | Consensus motif.
\<output_prefix>.down_counts.tab | Up counts.
\<output_prefix>.up_counts.tab | Down counts.
\<output_prefix>.seq_matrix.tab | Matrix of sequences.


* Example of generation of an additional FASTA output.

    python ./01_5primeCounter/5PrimeCounter.py --input_sites ./data/inputs_small/output_matrix_scan_NR3C1_chr20.txt --bam_file ./data/inputs_small/IMR90.rep2.uniq_mapped.chr20_sorted.bam --output_prefix ./output_5PrimeCounter_outputfa --genome_seq ./data/inputs_small/hg19_chr20.fa --output_seq

Output files

 File name | Content
--- | ---
\<output_prefix>.consensus.txt | Consensus motif.
\<output_prefix>.down_counts.tab | Up counts.
\<output_prefix>.up_counts.tab | Down counts.
\<output_prefix>.seq_matrix.tab | Matrix of sequences.
\<output_prefix>.fa | Sequences of regions of interest in FASTA format.


 
##### Use case 4 : motif comparison (overlay mode)


* Overlay plots (Multiple profile plots)

TODO


##### Other options

Example of ordering output by score instead of by coverage.

    python ./01_5primeCounter/5PrimeCounter.py --input_sites ./data/inputs_small/output_matrix_scan_NR3C1_chr20.txt --bam_file ./data/inputs_small/IMR90.rep2.uniq_mapped.chr20_sorted.bam --output_prefix ./output_5PrimeCounter_orderscore --order_by_score

Example of changing number of output sites.

    python ./01_5primeCounter/5PrimeCounter.py --input_sites ./data/inputs_small/output_matrix_scan_NR3C1_chr20.txt --bam_file ./data/inputs_small/IMR90.rep2.uniq_mapped.chr20_sorted.bam --output_prefix ./output_5PrimeCounter_100sites --number_of_sites 100

    python ./01_5primeCounter/5PrimeCounter.py --input_sites ./data/inputs_small/output_matrix_scan_NR3C1_chr20.txt --bam_file ./data/inputs_small/IMR90.rep2.uniq_mapped.chr20_sorted.bam --output_prefix ./output_5PrimeCounter_20pcsites --percent_of_sites 20



***

### ExoPlotter

Plot the profile from an output of 5PrimeCounter.

#### Use cases (same as for 5PrimeCounter)

##### Use case 1 : basic usage

Prefix from 5PrimeCounter. This will be used to find output files from 5PrimeCounter, namely a consensus sequence, and up and down 5' counts.

 * Command

    Rscript exoPlotter.R <output_prefix_of_5PrimeCounter.py>

 * Output files

 File name | Content
--- | ---
output_5PrimeCounter.maxstrand_heatmap.pdf | .
output_5PrimeCounter.profile.pdf | .
output_5PrimeCounter.cluster.pdf | .

##### Use case 2 : with permutations

As for use case 1, prefix is the same as the one used for 5PrimeCounter. This will be used to find output files generated by 5PrimeCounter, namely a consensus sequence, and up and down 5' counts.
Profiles also contain a summary of the permuted profiles computed previously, as well as a Wilcoxon Rank sum test.

 * Command

    Rscript exoPlotter.R <output_prefix_of_5PrimeCounter.py> perm

 * Additional output files

Same files as for Use Case 1, plus:

 File name | Content
 --- | ---
 output_5PrimeCounter.permut_matrix_10.values.tab | .


##### Use case 3 : using input genome

###### Input and command line

As for use case 1, prefix is the same as the one used for 5PrimeCounter. This will be used to find output files generated by 5PrimeCounter, namely a consensus sequence, and up and down 5' counts.
If a label is provided saying that genome sequence was available, exoPlotter also displays a graph serving as quality check for the motif/regions alignment.

 * Command

    Rscript exoPlotter.R <output_prefix_of_5PrimeCounter.py> genome_seq

 * Additional output files

Same files as for Use Case 1, plus:

 File name | Content
 --- | ---
 output_5PrimeCounter_seq.pdf | Quality check.


##### Use case 4

TODO

***

License
-------

TODO


***


