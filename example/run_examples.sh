########################################################################
# 
# This script runs example analysis as described in the README.
#
########################################################################


#------------------------------------------------
# Step 01
#------------------------------------------------
python ../python/matrixScanWS.py \
    --bed_file data/input.bed \
    --genome hg19 \
    --matrix_file data/matrix.tf \
    --uth_pval 0.001 \
    --output_file output_matrix_scan_0-001.txt

# use case 2:
python ../python/matrixScanWS.py \
    --bed_file data/input.bed \
    --genome hg19 \
    --matrix_file data/matrix.tf \
    --uth_pval 0.001 \
    --output_file output_matrix_scan_0-001.txt \
    --perm

#------------------------------------------------
# Step 02
#------------------------------------------------

#Use case 1

python ../python/5primeCounter.py \
    --input_sites output_matrix_scan_0-001.txt \
    --bam_file data/input.sorted.bam \
    --output_prefix output_5PrimeCounter

python ../python/5primeCounter.py \
    --input_sites output_matrix_scan_0-001.txt \
    --bam_file data/input.sorted.bam \
    --output_prefix ./output_5PrimeCounter \
    --size 100
    
python ../python/5primeCounter.py \
    --input_sites output_matrix_scan_0-001.txt \
    --bam_file data/input.sorted.bam \
    --output_prefix ./output_5PrimeCounter \
    --size 100 \
    --output_bed

# Use case 2

python ../python/5primeCounter.py \
    --input_sites output_matrix_scan_0-001.txt \
    --bam_file data/input.sorted.bam \
    --output_prefix ./output_5PrimeCounter_wperm \
    --perm

# Use case 3
python ../python/5primeCounter.py \
    --input_sites output_matrix_scan_0-001.txt \
    --bam_file data/input.sorted.bam \
    --output_prefix ./output_5PrimeCounter_wgenome \
    --genome_seq data/chr1.sub.fa

python ../python/5primeCounter.py \
    --input_sites output_matrix_scan_0-001.txt \
    --bam_file data/input.sorted.bam \
    --output_prefix  ./output_5PrimeCounter_wgenome \
    --genome_seq data/chr1.sub.fa \
    --size 100 \
    --output_seq

python ../python/5primeCounter.py \
    --input_sites output_matrix_scan_0-001.txt \
    --bam_file data/input.sorted.bam \
    --output_prefix ./output_5PrimeCounter_orderscore \
    --order_by_score

python ../python/5primeCounter.py \
    --input_sites output_matrix_scan_0-001.txt \
    --bam_file data/input.sorted.bam \
    --output_prefix ./output_5PrimeCounter_100sites \
     --number_of_sites 100

python ../python/5primeCounter.py \
    --input_sites output_matrix_scan_0-001.txt \
    --bam_file data/input.sorted.bam \
    --output_prefix ./output_5PrimeCounter_20pcsites \
    --percent_of_sites 20


#------------------------------------------------
# Step 03
#------------------------------------------------


