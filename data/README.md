
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

Code to execute inside folder ExoProfiler/data.

    wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.2bit
    wget http://hgdownload.cse.ucsc.edu/admin/exe/macOSX.x86_64/twoBitToFa
    chmod u+x twoBitToFa 
    ./twoBitToFa hg19.2bit hg19.fa
    mkdir -p inputs/genome
    mv hg19.fa inputs/genome/hg19.fa

### BAM files

All ChIP-seq and ChIP-exo data have been submitted to ArrayExpress and can be freely downloaded from the repository (upon publication !) using identifiers E-MTAB-2955 and E-MTAB-2956.

All BAM files must first be indexed with samtools.

ici mettre la ligne de commande pour indexer bam avec samtools

### BED files

In the article, we used as input ChIP-seq peaks (note that the pipeline can also take as input other types of peaks within which one want to analyse ChIp-exo signal)

### Motifs

Motifs are already available in the GitHub repository, but can also be downloaded from JASPAR.


