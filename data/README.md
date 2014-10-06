
Workflow to reproduce article figures
=======================

Folder structure to create
-------------------

    ExoProfiler/
      data/
        inputs/
          bam/
          bed/
          genome/
          motifs/


How to download input data
--------------

### Human genomic sequence (FASTA format) hg19 assembly

Code to be executed inside folder ExoProfiler/data/inputs, in order to download Human Genome (version 19) from UCSC.

    mkdir -p genome
    cd genome
    wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.2bit
    wget http://hgdownload.cse.ucsc.edu/admin/exe/macOSX.x86_64/twoBitToFa
    chmod u+x twoBitToFa 
    ./twoBitToFa hg19.2bit hg19.fa

Please note that other versions of twoBitToFa can be found [here](http://hgdownload.cse.ucsc.edu/admin/exe/).


### BAM files

All ChIP-seq and ChIP-exo data have been submitted to ArrayExpress and can be freely downloaded from the repository (upon publication !) using identifiers E-MTAB-2955 and E-MTAB-2956.

Cell lines:

* [IMR90](https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-2956/IMR90_GR_chip-exo.bam)
* [K562](https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-2956/K562_GR_chip-exo.bam)
* [U2OS](https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-2956/U2OS_GR_chip-exo.bam)
* CTCF, SRA ID: SRA044886 replicate 3, processed similarly to the original study (Rhee 2011).

BAM files should be downloaded into ExoProfiler/data/inputs/bam, and indexed with samtools.

    samtools index <onefile>.bam


### BED files

In the article, we used as input ChIP-seq peaks (note that the pipeline can also take as input other types of peaks within which one want to analyse ChIp-exo signal).

Cell lines:

* [IMR90](https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-2955/E-MTAB-2955.processed.1.zip/IMR90_GR_chip-seq_rep1_peaks.bed.gz)
* [K562](https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-2955/E-MTAB-2955.processed.1.zip/K562_GR_chip-seq_rep1_peaks.bed.gz)
* [U2OS](https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-2955/E-MTAB-2955.processed.1.zip/U2OS_GR_chip-seq_peaks.bed.gz)
* CTCF, GEO ID: GSM325895 converted to hg19 with UCSC liftOver.

These BED files should be downloaded into ExoProfiler/data/inputs/bed and unzipped.


### Motifs

Motifs files are already available in the GitHub repository (ExoProfiler/data/inputs/motif/), but can also be downloaded from JASPAR.


How to re-execute ExoProfiler 
--------------------------------------

Execute accompanying script `commands_article.py`.
