***
***
<a id="top"></a>

ExoProfiler
========================================================================
  
  
>
> Analyzing transcription factor (TF) binding sites using ChIP-exo data.  
> 

>   
> Stephan R. Starick*, Jonas Ibn-Salem*, Marcel Jurk, Céline Hernandez, Michael I. Love, Ho-Ryun Chung, Martin Vingron, Morgane Thomas-Chollier#, Sebastiaan H. Meijsing# **ChIP-exo signal associated with DNA-binding motifs provides structural and functional insights into the genomic binding of cooperating transcription factors.**, submitted
>

Table of Content
------------------------------------------------------------------------

1. [Requirements](#requirements)
1. [Pipeline presentation](#pipeline)
1. [Example analysis](#example)
1. [License](#license)

***
***

Requirements <a id="requirements"></a>
------------------------------------------------------------------------

[Back to top](#top)

The ExoProfiler pipeline is composed of three different steps.

1. MatrixScanWS
1. 5PrimeCounter
1. ExoPlotter


### MatrixScanWS

Tool developed in Python by Celine Hernandez.

Developed with:

* Python (2.7.5).

Package dependencies: 

* suds (0.4.1), numpy



### 5PrimeCounter

Tool developed in Python by Jonas Ibn-Salem.

Developed with:

* Python (2.7.3).

Package dependencies: 

* [numpy](https://pypi.python.org/pypi/numpy), [HTSeq](https://pypi.python.org/pypi/HTSeq) (see [installation guide](http://www-huber.embl.de/users/anders/HTSeq/doc/install.html)), [pysam](https://pypi.python.org/pypi/pysam) (HTSeq dependency)
* (optional) [pyfasta](https://pypi.python.org/pypi/pyfasta). If a reference genome is provided to calculate consensus sequences (see Use case 3), 5PrimeCounter also imports the 'pyfasta' package.


### ExoPlotter

Tool developed in R by Jonas Ibn-Salem, with code from Mike Love, Morgane Thomas-Chollier, Samuel Collombet and Celine Hernandez.

Developed with:

* R (2.14.1).

Package dependencies (available on CRAN): 

* gdata (NB: gtools will also be installed).
* ape 


***
***
<a id="pipeline"></a>

Pipeline presentation
------------------------------------------------------------------------

[Back to top](#top)

The full pipeline is composed of three different steps.

1. **MatrixScanWS (Optional step)**
   * Steps and RSAT tools involved (as web services)
   * Creation of a validation set by permuting the input matrix
2. **5PrimeCounter**
   * Use cases
    * Use case 1 : basic usage
    * Use case 2 : validation using permuted matrices
    * Use case 3 : QC and FASTA creation using an input genome
  * Additional information
3. **ExoPlotter**

***

### 1. MatrixScanWS (Optional step)

MatrixScanWS takes advantage of [RSAT](http://rsat.eu) web services  in order to convert a BED file into FASTA sequences, which are then scanned using a given matrix (transfac format).

#### Steps and RSAT tools involved (as web services)

1. Conversion of the input bed file (option -b or \--bed\_file) into a fasta file, using a reference genome (option -g or --genome). This step uses [fetch_sequences](http://www.rsat.eu/fetch-sequences_form.php).
2. Scan the FASTA sequences for motifs, provided as a transfac matrix (option -m or --matrix\_file). Only matches with a p-value below a threshold are reported (option
-u or \--uth\_pval). This step uses [matrix_scan](http://www.rsat.eu/matrix-scan-quick_form.cgi).
3. Matches are written in an output file (option -o or --output_file).

#### Creation of a validation set by permuting the input matrix

If the \--perm (or -p) option is included in the command, additional steps will be performed.
The input transfac matrix will be permuted 10 times using the tool [convert_matrix](http://www.rsat.eu/convert-matrix_form.cgi). Those permuted matrices should not be too similar to the input transfac matrix. This is checked using the tool [compare_matrices](http://www.rsat.eu/compare-matrices_form.cgi). They should not be duplicated either. This is 
checked again using 'compare_matrices'.
The permuted matrices will then be used by [matrix_scan](http://www.rsat.eu/matrix-scan-quick_form.cgi) to scan the FASTA sequences (computed in step 1/), using the same upper threshold for the p-value as in step 2/.
Finally, permuted matrices and outputs from 'matrix_scan' are written in the same folder as the output file of the input transfac matrix.  

  
***

  
### 2. 5PrimeCounter

  
5PrimeCounter analyses a BAM file from a ChIP-Exo experiment in the context of potential transcription factor binding sites (TFBS) presenting a given sequence motif of interest.  

Consequently, as input files, 5PrimeCounter needs a BAM file (and its BAI index) and a set of sequence motifs as created by MatrixScanWS for instance, or any output file of RSAT's tool 'matrix_scan'.

Please note that in case of motif hits on both strands in the same location, only the one with highest score will be considered.  
  
  
#### Use cases
  
  
##### Use case 1 : basic usage

Mandatory parameters :

Short name | Long name | Description |
--- | --- | --- |
-bam \<BAM_FILE> | \--bam_file \<BAM_FILE> | Bam file of the ChIP-exo experiment. Requires the index of the bam file with name \<BAM_FILE>.bai in the same folder. |
-i \<INPUT_SITES> | \--input_sites \<INPUT_SITES> | Input is a file containing predicted Transcription Factor Binding Sites in RSAT matrix-scan output format. |
-o \<OUTPUT_PREFIX> | \--output_prefix \<OUTPUT_PREFIX> | All output file names will have that prefix. Can include a path. |

Default window size around motifs is set to 60 bases. This parameter can be modified using the -s option.

Short name | Long name | Description |
--- | --- | --- |
-s \<SIZE> | \--size \<SIZE> | Window size around motif for which the profile will be computed. |
-ob | \--output\_bed | Write a BED file with the binding site regions defined by \--size and \--order_by_score. |

  
##### Use case 2 : validation using permuted matrices

Parameter to include :

Short name | Long name | Description
--- | --- | ---
&nbsp; | \--perm | Compute profiles from matrix-scan results for permuted matrices. 5PrimeCounter searches for all files with the same name as the \<INPUT_SITES> file, plus the tag '\_perm' and a number. Files must be located in the same folder as \<INPUT_SITES>.

  
##### Use case 3 : QC and FASTA creation using an input genome

Parameter to include :

Short name | Long name | Description
--- | --- | ---
-g \<GENOME\_SEQ> | \--genome\_seq \<GENOME_SEQ> | Reference genome sequence in FASTA format. If this optional argument is given, the consensus sequence of the motif is plotted at the bottom of profile and heatmap plots. Moreover sequences for all binding regions will be written to an integer-encoded matrix file. Note, if a fasta file is read for the first time, an index is built in the same directory for faster access. First execution can thus be slower.

Dependent parameter :

Short name | Long name | Description
--- | --- | ---
-of | \--output\_seq | Write a genomic sequences in FASTA format for the binding site regions defined by \--size and \--order_by_score. Needs a genome in FASTA format to be provided using option \--genome\_seq.


  
  
#### Aditional information

List of other possible command line options.

 * Change order of output motifs

 Short name | Long name | Description
 --- | --- | ---
 &nbsp; | &nbsp; | By default, output regions are sorted by occupancy level (number of total read counts).
 -os | \--order_by_score  | Sort output regions by score instead of occupancy level.

 * Change number of output motifs

 Short name | Long name | Description
 --- | --- | ---
 -n \<NUMBER_OF_SITES> | \--number_of_sites \<NUMBER_OF_SITES> | Number of sites to be considered. For a given N take only the top N sites by occupancy level (or motif score if -os is set).
 -p \<PERCENT_OF_SITES> | \--percent_of_sites \<PERCENT_OF_SITES> | Percent of sites to be considered. For a given P take only the top P percent sites by occupancy level (or motif score if -os is set).

 * Miscellaneous

 Short name | Long name | Description
 --- | --- | ---
 -h | \--help | Show help message and exit.
 -fs | \--flip_strand  | Flip the strand of motif matches from '+' to '-' and from '-' to '+' for all input sites.

A distance, manually computed, can be provided to shift sites by a given base pair number.

 Short name | Long name | Description
 --- | --- | ---
 -sd \<SHIFT_DIST> | \--shift_dist \<SHIFT_DIST> | Shift sites by given distance (in bp) to the right (if positive) or to the left (if negative).


  
***

  
### 3. ExoPlotter ###

ExoPlotter generates a set of PDF plots using output files generated by 5PrimeCounter. This R script needs as input the prefix provided to 5PrimeCounter. Two other flags can be added to get additional plots, also written in PDF files. 

A typical command line would be :

> Rscript exoPlotter.R \<OUTPUT_PREFIX> genome_seq perm

Where :

 * \<OUTPUT_PREFIX> is the output prefix used in the 5PrimeCounter command line.
 * 'genome_seq' is a fixed label. If present, 5PrimeCounter was run with option \--genome_seq, and corresponding output files can be found using \<OUTPUT_PREFIX>.
 * 'perm' is a fixed label. If present, 5PrimeCounter was run with option \--perm, and corresponding output files can be found using \<OUTPUT_PREFIX>.

***
***

Example analysis
------------------------------------------------------------------------
<a id="example"></a>

[Back to top](#top)

1. MatrixScanWS (Optional step)
   * Input files
   * Use cases
    * Use case 1 : basic usage
    * Use case 2 : scan plus validation set generation
2. 5PrimeCounter
   * Use cases
    * Use case 1 : basic usage
    * Use case 2 : validation using permuted matrices
    * Use case 3 : QC and FASTA creation using an input genome
   * Other examples
3. ExoPlotter
   * Use cases (same as for 5PrimeCounter)
    * Use case 1 : basic usage
    * Use case 2 : validation using permuted matrices
    * Use case 3 : QC and FASTA creation using an input genome

***

### MatrixScanWS (Optional step)

Using web services, scan a set of sequences for a motif provided as a matrix.

#### Input files

BED file : ```input.bed```  
Output of a ChIP-Seq experiment from a human sample, i.e. list of peaks.

Matrix (transfac) : ```matrix.tf```  
Matrix describing a motif as available in [JASPAR database](http://jaspar.genereg.net/cgi-bin/jaspar_db.pl?ID=MA0113.1&rm=present&collection=CORE).

#### Use cases


##### Use case 1: basic usage

    python matrixScanWS.py --bed_file input.bed --genome hg19 --matrix_file matrix.tf --uth_pval 0.001 --output_file ./output_matrix_scan_0-001.txt

Output file

* RSAT 'matrix_scan' output : ```./output_matrix_scan_0-001.txt```

##### Use case 2: scan plus validation set generation

Same command with --perm option.

    python matrixScanWS.py --bed_file input.bed --genome hg19 --matrix_file matrix.tf --uth_pval 0.001 --output_file ./output_matrix_scan_0-001.txt --perm

Output files

* RSAT 'matrix_scan' output : ```./output_matrix_scan_0-001.txt```
* 10 permuted versions of the input matrix, in the same folder as the output : ```matrix_perm1.tf```, ```matrix_perm2.tf```, ```matrix_perm3.tf```, etc.
* 10 'matrix_scan' outputs computed from these permuted versions of the input matrix, in the same folder as the output : ```output_matrix_scan_perm1.txt```, ```output_matrix_scan_perm2.txt```, ```output_matrix_scan_perm3.txt```, etc.



***


### 5PrimeCounter

5PrimeCounter performs a profile-based analysis by calculating a 5' coverage profile, given an BAM file from a ChIP-exo experiment (accompanied by its index in BAI) and a list of potential TF binding sites (TFBS) from motif matching analysis.


#### Use cases


##### Use case 1 : basic usage

###### Input files

* BAM file : ```input.bam```  
Sequence read alignment data in binary format.

* BAI file : ```input.bam.bai```  
Index of the BAM file, created with `samtools index`. This index must be located in the same folder as the BAM file.

* RSAT 'matrix_scan' output : ```output_matrix_scan.txt```  
Output from RSAT 'matrix_scan', or directly generated by MatrixScanWS.


###### Command

Basic usage of 5PrimeCounter implies providing at least three parameters:

 * --input_sites for a matrix in matrix-scan format
 * --bam_file for a BAM file
 * --output_prefix to provide folder name and a prefix for output files. In following examples: './output_5PrimeCounter'.

Example of a simple command line :

    python 5PrimeCounter.py --input_sites output_matrix_scan.txt --bam_file input.bam --output_prefix ./output_5PrimeCounter

5PrimeCounter assumes that the 'input sites' file is an output from `matrix_scan`

Example of a case where the size of the window around the motif has been changed with option --size. Default window size is set to 60.

    python 5PrimeCounter.py --input_sites output_matrix_scan.txt --bam_file input.bam --output_prefix ./output_5PrimeCounter --size 100

Example of generation of an additional BED output, in case we changed the size of our region of interest.

    python 5PrimeCounter.py --input_sites output_matrix_scan.txt --bam_file input.bam --output_prefix ./output_5PrimeCounter --size 100 --output_bed


###### Output files

 File name | Content
--- | ---
\<output_prefix>.consensus.txt | Consensus motif.
\<output_prefix>.down_counts.tab | Up counts.
\<output_prefix>.up_counts.tab | Down counts.

 
##### Use case 2 : validation using permuted matrices 

Following Use Case 2 of MatrixScanWS a set of permuted motifs can be generated. One can then run 5PrimeCounter on all the generated matrix-scan outputs.

    python 5PrimeCounter.py --input_sites output_matrix_scan.txt --bam_file input.bam --output_prefix ./output_5PrimeCounter_wperm --perm

* Output files

Same files as for Use Case 1, times ten, for each permuted motif.


##### Use case 3 : QC and FASTA creation using an input genome

Genome files in FASTA format can be quite voluminous. Default behaviour of 5PrimeCounter doesn't need an input genome. But if provided, additional information can be computed, and another quality check plot can be drawn with ExoPlotter. Output file in FASTA can also be written.

###### Input files

* BAM file : ```input.bam```  
Sequence read alignment data in binary format.

* BAI file : ```input.bam.bai```  
Index of the BAM file, created with `samtools index`. This index must be located in the same folder as the BAM file.

* RSAT 'matrix\_scan' output : ```output_matrix_scan.txt```  
Output from RSAT tool ```matrix_scan```, or directly generated by MatrixScanWS.

* Reference genome : ```genome.fa```  
A FASTA file as downloaded from, for instance, UCSC's FTP site : ftp://hgdownload.cse.ucsc.edu/goldenPath/ Please note that files need to be decompressed into FASTA in order to be read by 5PrimeCounter.

###### Command

    python 5PrimeCounter.py --input_sites output_matrix_scan.txt --bam_file input.bam --output_prefix ./output_5PrimeCounter_wgenome --genome_seq ./genome.fa


###### Output files

 File name | Content
--- | ---
\<output_prefix>.consensus.txt | Consensus motif.
\<output_prefix>.down_counts.tab | Up counts.
\<output_prefix>.up_counts.tab | Down counts.
\<output_prefix>.seq_matrix.tab | Matrix of integer-encoded sequences which can be used by ExoPlotter in the next step of the pipeline.


###### Example of generation of an additional FASTA output.

    python 5PrimeCounter.py --input_sites output_matrix_scan.txt --bam_file input.bam --output_prefix ./output_5PrimeCounter_wgenome --genome_seq ./genome.fa --size 100 --output_seq

Output files

 File name | Content
--- | ---
\<output_prefix>.consensus.txt | Consensus motif.
\<output_prefix>.down_counts.tab | Up counts.
\<output_prefix>.up_counts.tab | Down counts.
\<output_prefix>.seq_matrix.tab | Matrix of sequences.
\<output_prefix>.fa | Sequences of regions of interest (100 bp) in FASTA format.




##### Other examples

Example of ordering output by score instead of by coverage.

    python 5PrimeCounter.py --input_sites output_matrix_scan.txt --bam_file input.bam --output_prefix ./output_5PrimeCounter_orderscore --order_by_score

Example of changing number of output sites.

    python 5PrimeCounter.py --input_sites output_matrix_scan.txt --bam_file input.bam --output_prefix ./output_5PrimeCounter_100sites --number_of_sites 100

    python 5PrimeCounter.py --input_sites output_matrix_scan.txt --bam_file input.bam --output_prefix ./output_5PrimeCounter_20pcsites --percent_of_sites 20



***

### ExoPlotter

Plot the profile from an output of 5PrimeCounter.

#### Use cases (same as for 5PrimeCounter)

##### Use case 1 : basic usage

Prefix from 5PrimeCounter. This will be used to find output files from 5PrimeCounter, namely a consensus sequence, and up and down 5' counts.

*Command*

    Rscript exoPlotter.R <output_prefix_of_5PrimeCounter.py>

*Output files*

 File name | Content
--- | ---
\<output_5PrimeCounter>.maxstrand_heatmap.pdf | .
\<output_5PrimeCounter>.profile.pdf | .
\<output_5PrimeCounter>.cluster.pdf | .

##### Use case 2 : validation using permuted matrices 

As for use case 1, prefix is the same as the one used for 5PrimeCounter. This will be used to find output files generated by 5PrimeCounter, namely a consensus sequence, and up and down 5' counts.  
Profiles also contain a summary of the permuted profiles computed previously, as well as a Wilcoxon Rank sum test.

*Command*

    Rscript exoPlotter.R <output_prefix_of_5PrimeCounter.py> perm

*Additional output files*

Same files as for Use Case 1, plus:

 File name | Content
 --- | ---
 \<output_5PrimeCounter>.permut_matrix_10.values.tab | Contains the resulting p-value of the Wilcoxon rank sum test.


##### Use case 3 : QC and FASTA creation using an input genome

###### Input and command line

As for use case 1, prefix is the same as the one used for 5PrimeCounter. This will be used to find output files generated by 5PrimeCounter, namely a consensus sequence, and up and down 5' counts.
If a label is provided saying that genome sequence was available, exoPlotter also displays a graph serving as quality check for the motif/regions alignment.

*Command*

    Rscript exoPlotter.R <output_prefix_of_5PrimeCounter.py> genome_seq

*Additional output files*

Same files as for Use Case 1, plus:

File name | Content
--- | ---
\<output\_5PrimeCounter\>\_seq.pdf | Quality check.



***
***

License
---------------------------------
<a id="license"></a>

[Back to top](#top)

<a rel="license" href="http://creativecommons.org/licenses/by-nc-sa/4.0/"><img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by-nc-sa/4.0/88x31.png" /></a><br /><span xmlns:dct="http://purl.org/dc/terms/" href="http://purl.org/dc/dcmitype/Text" property="dct:title" rel="dct:type">ExoProfiler</span> by <span xmlns:cc="http://creativecommons.org/ns#" property="cc:attributionName">MPIMG & IBENS</span> is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by-nc-sa/4.0/">Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License</a>.<br />Based on a work at <a xmlns:dct="http://purl.org/dc/terms/" href="https://github.com/ComputationalSystemsBiology/ExoProfiler" rel="dct:source">https://github.com/ComputationalSystemsBiology/ExoProfiler</a>.


***


