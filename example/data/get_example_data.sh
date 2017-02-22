########################################################################
# 
# This scirpt downloads and formates some example input data for the 
# ExoProfiler pipeline.
# 
# Note, this is just to document the creation of example data. You do
# not need to run this script in order to run the examples. 
# The data created here is already provided in this repository. 
#
########################################################################


#------------------------------------------------
# get .bed file with 100 rows
#------------------------------------------------
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-2955/E-MTAB-2955.processed.1.zip/IMR90_GR_chip-seq_rep1_peaks.bed.gz

gunzip IMR90_GR_chip-seq_rep1_peaks.bed.gz

tail -n +2 IMR90_GR_chip-seq_rep1_peaks.bed \
    |head -n 100 \
    > input.bed

#------------------------------------------------
# get motif from /data dir
#------------------------------------------------
    
cp ../../data/inputs/motif/MA0113.2.tf matrix.tf


#------------------------------------------------
# get .bam files
#------------------------------------------------

# extend regions in input.bed for filtering
mysql --user=genome \
    --host=genome-mysql.cse.ucsc.edu -A \
    -e "select chrom, size from hg19.chromInfo" \
     > hg19.genome

bedtools slop \
    -i input.bed \
    -g hg19.genome \
    -b 100 \
    > input.window.bed

# download full bam file from arrayexpress
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-2956/IMR90_GR_chip-exo.bam

# sort and index 
samtools sort \
    IMR90_GR_chip-exo.bam  \
    IMR90_GR_chip-exo.sorted 

samtools index IMR90_GR_chip-exo.sorted.bam

# filter .bam file to have only regions in .bed file
bedtools intersect \
    -abam IMR90_GR_chip-exo.sorted.bam \
    -b input.window.bed \
    > input.bam 

# sort and index 
samtools sort input.bam input.sorted
samtools index input.sorted.bam
   
#------------------------------------------------
# get fasta file of human chromosome 1.
#------------------------------------------------

wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/chr1.fa.gz
gunzip chr1.fa.gz

# filter fast file 

maxPos=$(cat input.window.bed |cut -f 3 |sort -n -r |head -n 1)
echo -e "chr1\t0\t${maxPos}\tchr1" > range.bed

bedtools getfasta \
	-fi chr1.fa \
	-bed range.bed \
	-name \
	-fo - \
	| fold -w 60 \
	> chr1.sub.fa

	
#------------------------------------------------
# Clean up by removing temporary files
#------------------------------------------------
rm IMR90_GR_chip-seq_rep1_peaks.bed
rm input.window.bed
rm IMR90_GR_chip-exo.bam
rm IMR90_GR_chip-exo.sorted.bam
rm IMR90_GR_chip-exo.sorted.bam.bai
rm input.bam
rm chr1.fa
rm chr1.fa.fai
rm range.bed

#------------------------------------------------
# setup latest HTSeq version
#------------------------------------------------

#~ sudo pip install Cython
#~ sudo pip install 'matplotlib>=1.4'
#~ sudo pip install https://github.com/simon-anders/htseq/archive/release_0.7.0.tar.gz


