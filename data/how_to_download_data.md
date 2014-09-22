
How to download data
=======================

Excepted genome sequences, all data has been submitted to ArrayExpress and can be freely downloaded from the repository using identifiers E-MTAB-2731, E-MTAB-2954, E-MTAB-2955 and E-MTAB-2956.


Genome
--------------

 wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.2bit
 wget http://hgdownload.cse.ucsc.edu/admin/exe/macOSX.x86_64/twoBitToFa
 chmod u+x twoBitToFa 
 ./twoBitToFa hg19.2bit hg19.fa
 mkdir -p inputs/genome
 mv hg19.fa inputs/genome/hg19.fa

Motifs
--------------

Motifs are already available in the GitHub repository, but can also be downloaded from JASPAR.

BAM files
--------------

All BAM files must have been indexed.

BED files
--------------
