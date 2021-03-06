#!/usr/bin/env python
"""
5PrimeCounter
=============

5PrimeCounter analyses a BAM file from a ChIP-Exo experiment in the context of potential transcription factor binding sites (TFBS) presenting a given sequence motif of interest.

Consequently, as input files, 5PrimeCounter needs a BAM file (and its BAI index) and a set of sequence motifs as created by MatrixScanWS for instance, or any output file of RSAT's tool 'matrix_scan'.

Please note that in case of motif hits on both strands in the same location, only the one with highest score will be considered.

Use cases

    Use case 1 : basic usage
    Use case 2 : validation using permuted matrices
    Use case 3 : QC and FASTA creation using an input genome


Developed with:

    Python (2.7.3).

Package dependencies:

    numpy, HTSeq (see installation guide), pysam (HTSeq dependency)
    (optional) pyfasta. If a reference genome is provided to calculate consensus sequences (see Use case 3), 5PrimeCounter also imports the 'pyfasta' package.

Tool developed in Python by Jonas Ibn-Salem.

# Last update: 26.09.2014 (CEH)

"""
epilog="""
15.10.13 Jonas Ibn-Salem <ibnsalem@molgen.mpg.de>
"""
import argparse
import os.path
# import sys
import numpy as np
import HTSeq    # For installation see: http://www-huber.embl.de/users/anders/HTSeq/doc/install.html
import subprocess


## Read and check input parameters

def commandline():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter, epilog=epilog)
    
    # Use case 1
    parser.add_argument("-bam", "--bam_file", type=str, required=True, help="Bam file of the ChIP-exo experiment. Requires the index of the bam file with name <BAM_FILE>.bai in the same folder.")
    parser.add_argument("-i", "--input_sites", type=str, required=True, help="Input is a file containing predicted Transcription Factor Binding Sites in RSAT matrix-scan output format.")
    parser.add_argument("-o", "--output_prefix", type=str, required=True, help="All output file names will have that prefix. Can include a path.")

    parser.add_argument("-s", "--size", type=int, default=60, help="Window size around motif for which the profile will be computed.")
    parser.add_argument("-ob", "--output_bed", action="store_true", help="Write a BED file with the binding site regions defined by --size and --order_by_score.")
    parser.add_argument("-if", "--input_format", type=str, choices=['matrix-scan'], default="matrix-scan", help="Input format, 'matrix-scan' output (default).")

    # Use case 2
    parser.add_argument("-pm", "--perm", action="store_true", help="Compute profiles from matrix-scan results for permuted matrices. 5PrimeCounter searches for all files with the same name as the <INPUT_SITES> file, plus the tag '_perm' and a number. Files must be located in the same folder as <INPUT_SITES>.")

    # Use case 3
    parser.add_argument("-g", "--genome_seq", type=str, help="Reference genome sequence in FASTA format. If this optional argument is given, the consensus sequence of the motif is plotted at the bottom of profile and heatmap plots. Moreover sequences for all binding regions will be written to an integer-encoded matrix file. Note, if a fasta file is read for the first time, an index is built in the same directory for faster access. First execution can thus be slower.")
    parser.add_argument("-of", "--output_seq", action="store_true", help="Write a genomic sequences in FASTA format for the binding site regions defined by --size and --order_by_score. Needs a genome in FASTA format to be provided using option --genome_seq.")

    # Other options
    parser.add_argument("-os", "--order_by_score", action="store_true", help="By default, output regions are sorted by occupancy level (number of total read counts). This option sorts output regions by score instead of occupancy level.")
    parser.add_argument("-n", "--number_of_sites", type=int, help="Number of sites to be considered. For a given N take only the top N sites by occupancy level (or motif score if -os is set).")
    parser.add_argument("-p", "--percent_of_sites", type=float, help="Percent of sites to be considered. For a given P take only the top P percent sites by occupancy level (or motif score if -os is set).")
    parser.add_argument("-d", "--down_sample_sites", type=int, help="Down sample input sites. For a given D sample D sites randomly.")
    parser.add_argument("-sd", "--shift_dist", type=int, help="Shift sites by given distance (in bp) to the right (if positive) or to the left (if negative).")
    parser.add_argument("-fs", "--flip_strand", action="store_true", help="Flip the strand of motif matches from '+' to '-' and from '-' to '+' for all input sites.")

    return parser.parse_args()


## Parse input files

def parse_matrix_scan(inFile):
    """
    Parses the RSAT matrix-scan output file and returns a list of dicts ordered by score    
    Keeps only one hit if 2 hits at the same position on both strands.
    Genomic coordinates are transfromed from 1-based (in matrix-scan output format) 
    into zero-based half open (like BED format) for internal representation and HTSeq compatibility.
    Assumes the "galaxy format" for sequence IDs in the first column.
    """
    
    # dict for unique sites
    unique_sites = {}
    
    # count number of total input sites
    n_sites = 0
    
    for line in open(inFile):
        
        # ignore comment lines
        if not line.startswith(';') and not line.startswith('#'):
            
            sp = line.strip().split('\t')
            loc = sp[0]

            chr = loc.split('_')[1]
            # the peak-coordinates are assumed now again ONE-based in matrix-scan output format! 
            peak_start = int(loc.split('_')[2]) - 1
            
            # the motif coordinates are ONE-based and relative to peak start in the matrix-scan output file format
            start = peak_start + int(sp[4]) - 1
            end = peak_start + int(sp[5])
            
            strand = sp[3].replace('DR', '+').replace('D', '+').replace('R','-')
            score = float(sp[7])
            
            # one based locus coordinates of the motif:
            motif_loc = chr + ":" + str(start+1) + "-" + str(end)
                
            # keep it only if score is greater than sits with same location
            if motif_loc not in unique_sites or score > unique_sites[motif_loc]["score"]:
            
                type = sp[1]
                ft_name = sp[2]
                seq = HTSeq.Sequence(sp[6], loc)
                
                # append region as dict with all annotations
                unique_sites[motif_loc] = {"chr":chr, "start":start, "end":end, \
                    "strand":strand, "score":score, "type":type, \
                    "ft_name":ft_name, "motif_seq":seq, "motif_loc":motif_loc, 
                    "seq_id":loc, "name":motif_loc}
            
            n_sites += 1 # increase counter for total number of sites
    
    # get list of sites sorted by motif score
    sorted_sites = sorted(unique_sites.values(), cmp=lambda x,y: cmp(x['score'], y['score']), reverse = True )
    
    print "INFO: Read {0} of {1} input regions.".format(len(sorted_sites), n_sites)
    
    return sorted_sites


def add_center(sites, size):
    """ 
    Add center of motif to sites. In case of even motifs (real center
    is between two bases) the closest base upstream from the center is 
    chosen.
    Add also extended region information as "ext_start" and "ext_end" positions.
    """
    for s in sites:
        
        # check if motif length is even
        even = (s["end"] - s["start"]) % 2 == 0

        # calculate center coordinate of binding site
        center = (s["start"]+s["end"]-1)/2
        
        # in even case on the reverse strand the center is 
        # the closesd base upstream of the real center
        if even and s["strand"] == '-':
            center += 1
            
        s["center"] = center
        s["ext_start"] = center - size/2
        s["ext_end"] = center + size/2
        # get 1-based genomic location in the format "chr:start-end"
        s["location"] = s["chr"] + ":" + str(s["ext_start"]+1) + "-" + str(s["ext_end"])

    return sites


def reads_profile(regions, bam_file, size):
    """
    Parses reads from BAM file and adds number of forward and reverse 
    5' coverage counts per position of each region.
    This function depend on the HTSeq package for fast parsing of read infromation from BAM files.
    """

    print "INFO: Begin to parse reads from BAM file for n={0} regions.".format(len(regions))
    
    # Open BAM file:
    bamHandle = HTSeq.BAM_Reader(bam_file)
    # get list of available chromosoms
    chromosomes = set([chr['SN'] for chr in bamHandle.get_header_dict()['SQ']])

    for i, reg in enumerate(regions):
        
        center = reg["center"]
        
        # initialize read-counts for all positions of this region
        up_counts = size * [0]
        down_counts = size * [0]
        
        # check if chr of region is available in BAM file:
        if reg["chr"] in chromosomes:

            # get GenomicInterval object. extend it by +-1 to for including reads on negative strand inside the interval
            iv = HTSeq.GenomicInterval( reg["chr"], max(0, reg["ext_start"]-1), reg["ext_end"]+1, reg["strand"] )

            # iterate over all reads mapping to that region (interval)
            for aln in bamHandle[ iv ]:
                            
                # consider motif on positiv stand
                if reg["strand"] == '+':
    
                    dist = aln.iv.start_d - center
                    pos = dist + size/2
                    
                    if pos >= 0 and pos < size:
                        if aln.iv.strand == '+': up_counts[pos] += 1
                        if aln.iv.strand == '-': down_counts[pos] += 1
    
                if reg["strand"] == '-':
                    
                    dist = -1 * (aln.iv.start_d - center)
                    pos = dist + size/2
    
                    if pos >= 0 and pos < size:
                        if aln.iv.strand == '+': down_counts[pos] += 1
                        if aln.iv.strand == '-': up_counts[pos] += 1
                
        # add counts to region dictionary:
        reg["up_counts"] = up_counts
        reg["down_counts"] = down_counts
    
    print "INFO: Finished parsing of BAM file."
                
    return regions

## Write functions

def write_counts(regions, output_prefix):
    """ Writes count matrix as TAB seperated file to output file. """

    # get number of rows:
    n = len(regions)
    
    upHandle = open(output_prefix + ".up_counts.tab", 'w')
    downHandle = open(output_prefix + ".down_counts.tab", 'w')

    # iterate over rows
    for reg in regions:
        
        upHandle.write('\t'.join([reg["location"]+"_"+reg["name"]] + [str(c) for c in reg["up_counts"]]) + '\n')
        downHandle.write('\t'.join([reg["location"]+"_"+reg["name"]] + [str(c) for c in reg["down_counts"]]) + '\n')
        
    upHandle.close()
    downHandle.close()


def write_region_to_bed(regions, outFile, size=60):
    """Regions to outFile in BED format"""
    
    with open(outFile, 'w') as outHandle:
    
        for reg in regions:
    
            outHandle.write("\t".join([str(c) for c in [
            reg["chr"], reg["center"]-size/2, reg["center"]+size/2, 
            reg["name"] if "name" in reg else ".", 
            reg["score"] if "score" in reg else ".", 
            reg["strand"] 
            ]]) + '\n')


def write_fasta(regions, outFile):
    """writes sequences of regions to fasta file"""
    
    with open(outFile, 'w') as outHandle:
        for reg in regions:
            outHandle.write(">" + reg["location"]+"_"+reg["strand"] + "\n")
            outHandle.write(reg["ext_seq"] + "\n")

## Read sequences

def parse_sequences(sites, size, fasta_file):
    """Adds the binding site sequences extende to 'size' per row (decoded as A=0, C=1, G=2, T=3) to each input region."""
    from pyfasta import Fasta  # Fasta package is needed to fetch sequences from genome fasta file
            
    print "INFO: Begin to fetch sequences...."
    
    f = Fasta(fasta_file, key_fn=lambda key: key.split()[0])

    for i, reg in enumerate(sites):
        
        start = reg["ext_start"]
        end = reg["ext_end"]
        
        # if motif on negativ strand, shift region by +1 to account for zero based half-open intervals
        if reg["strand"] == '-':
            start += 1
            end += 1
        
        seq = f.sequence({"chr":reg["chr"], "start":start, "stop":end}, one_based=False)

        # Note, the 'strand':reg["strand"] argument for f.sequence does not work, there seems to be a bug in the pyfasta/fasta.py code.
        seq = seq.upper()
 
        # if motif on negative strand, convert seq to reverse complement
        if reg["strand"] == '-': 
            seq = reverse_complement(seq)
        
        # add sequence to region dict
        reg["ext_seq"] = seq
        
    print "INFO: Finished sequences."
    return regions 


def reverse_complement(seq):
    """ returns the reverse complement of seq"""
    rep_dict = {"A":"T", \
                "C":"G", \
                "G":"C", \
                "T":"A"}
    revcomp = ""
    for i, base in enumerate(seq):
        if base in rep_dict:
            revcomp += rep_dict[base]
        else:
            revcomp += base

    return revcomp[::-1]


def get_consensus(sites, seq_type, m=-1):
    """return a string as consesus sequence """
    
    bases = ['A', 'C', 'G', 'T', 'N']
    n = len(sites)  # number of sites    
    if m == -1 and sites: m = len(sites[0][seq_type])
        
    # initialize count array
    # rows correspond to positions in motif sequence
    # columns correspond to bases: "A", "C", "G", "T", and "N"
    counts = np.zeros(( m, 5 ), np.int )
    consenus = ""
    for s in sites:
    
        # convert seq to HTSeq.Sequence object
        seq = HTSeq.Sequence(str(s[seq_type]))
        # count bases to counts array
        seq.add_bases_to_count_array( counts )

    base_idx =  np.argmax(counts, 1)
    
    for i in range(m):
        
        # test if at least 75% of sites have same base:
        if n>0 and counts[i, base_idx[i]]/float(n) > 0.75:
            consenus += bases[base_idx[i]]

        # test if at least 50% of sites have same base:
        elif n>0 and counts[i, base_idx[i]]/float(n) > 0.5:
            consenus += bases[base_idx[i]].lower()
        else:
            consenus += "."
    
    return consenus


## Write functions

def write_consensus(consenus, size, outFile):
    """ write extende consenus sequence to outFile"""
    if len(consensus) == size : 
        ext_consensus = consensus
    else:
        l = len(consensus)
        before =  (size - l) / 2 + 1
        after = np.ceil( (size - l)/2.0 ) - 1
        ext_consensus = before * '.' + consenus + after * '.'
        
    with open(outFile, 'w') as outHandle:
        outHandle.write(ext_consensus  + "\n")


def write_seq_matrix(seq_matrix, outFile):
    """ writes for each region the genomic sequence encoded as integers to tab seperated file"""
    base2int = {"A":"0", "C":"1", "G":"2", "T":"3", "N":"4"}

    with open(outFile, 'w') as outHandle:
        for reg in regions:
            outHandle.write('\t'.join([reg["location"]+"_"+reg["name"]] + [base2int[b] for b in reg["ext_seq"]]) + '\n')


## Ordering functions

def order_by_read_counts(regions):
    """Reorders the input list of regions by number of total read counts"""
    return sorted(regions, cmp=lambda x,y: cmp(sum(x['up_counts']+x['down_counts']), sum(y['up_counts']+y['down_counts'])), reverse = True )

def order_by_score(regions, reverse=True):
    """Reorders the list of input regions by score"""
    return sorted(regions, cmp=lambda x,y: cmp(x['score'], y['score']), reverse = reverse )


## Other functions

def flip_strands(regions):
    """ Filp the strand of all input regions """
    print "INFO: Flip strand of all input sites."
    
    for i in range(len(regions)):
        
        strand = regions[i]["strand"]
        
        if strand == "+":
            regions[i]["strand"] = "-"
        elif strand == "-":
            regions[i]["strand"] = "+"

    return regions

def shift_sites(regions, shift_dist):
    """Shift all regions by indicated distance to the right (if positive) or to the left (if negative) """

    for s in regions:
        for coord in ["start", "center", "end", "ext_start", "ext_end"]:
            
            # in case of motif on negative stand, shift in oposite direction
            if s["strand"] == '-':
                s[coord] -= shift_dist
            else:
                s[coord] += shift_dist
            
        # addjust 1-based genomic location in the format "chr:start-end"
        s["location"] = s["chr"] + ":" + str(s["ext_start"]+1) + "-" + str(s["ext_end"])    


## Main 

if __name__ == "__main__":

    # read commandline argumets
    args = commandline()
    
    # test validity of other arguments:
    if args.number_of_sites and args.percent_of_sites: 
        sys.exit("ERROR: '--number_of_sites' and '--percent_of_sites' arguments are mutually exclusive. Exit now.")

    if args.perm:
        # Split between real name and extension (matrix)
        (matrix_filename_label, matrix_filename_ext) = os.path.splitext(args.input_sites)
        nb_permutations = 10
        permuted_res = []
        permuted_suffix = []
        for index in xrange(nb_permutations):
            permuted_res.append(os.path.join(matrix_filename_label+'_perm'+str(index+1)+matrix_filename_ext))
            permuted_suffix.append('_perm'+str(index+1))

    files_to_analyze = [args.input_sites]
    suffixes = ['']
    if args.perm:
        files_to_analyze +=  permuted_res
        suffixes += permuted_suffix
    for index in xrange(len(files_to_analyze)):

        file_to_analyze = files_to_analyze[index]

        if args.input_format.lower() == "matrix-scan":
            # parse matrix-scan results
            sites = parse_matrix_scan(file_to_analyze)            
            # extend sites to region of given size:
            regions = add_center(sites, args.size)

        else:
            sys.exit("ERROR: INPUT_FORMAT shuld be one of 'matrix-scan' or 'bed'. Exit now.")
    
        # if option 'flip_strand' is given, flip strand of all input sites:
        if args.flip_strand:
            regions = flip_strands(regions)
    
        # if option shift_dist is set, shift sites by given distance:
        if args.shift_dist:
            shift_sites(regions, args.shift_dist)
        
        # down sample sites
        if args.down_sample_sites:
            import random
            regions = random.sample(regions, args.down_sample_sites)
            print "INFO: Input sites were down sampled to {0} regions.".format(len(regions))

        # parse 5' coverage counts from BAM file:
        regions = reads_profile(regions, args.bam_file, args.size)
        
        # order output regions by motif score or exo-read occupancy level
        if args.order_by_score:        
            regions = order_by_score(regions)
        else:
            regions = order_by_read_counts(regions)
        
        # take only a subset of top p percent sites if such an argument is given:
        if args.percent_of_sites:
            args.number_of_sites = int(args.percent_of_sites * len(regions)/100) 

        # take only a subset of top N sites if such an argument is given:
        if args.number_of_sites:
            regions = regions[:args.number_of_sites]
    
        # write the 5' coverage count data to ouput files
        write_counts(regions, args.output_prefix+suffixes[index])    
        
        # fetch genomic sequences, if reference seq is given and calculate consensus sequence.
        if args.genome_seq:
            regions = parse_sequences(regions, args.size, args.genome_seq)
            write_seq_matrix(regions, args.output_prefix+suffixes[index] + ".seq_matrix.tab")
        
            consensus = get_consensus(regions, "ext_seq", args.size)
        else:
            consensus = get_consensus(regions, "motif_seq")
    
        # write consensus seq to output file
        write_consensus(consensus, args.size, args.output_prefix+suffixes[index] + ".consensus.txt")
	
        # write extended regions to BED file:
        if args.output_bed:
            write_region_to_bed(regions, args.output_prefix+suffixes[index] + ".bed", args.size)
    
        # write genomic sequences to fasta file:
        if args.output_seq:
            if args.genome_seq:
                write_fasta(regions, args.output_prefix+suffixes[index] + ".fa")
            else:
                sys.exit("ERROR: Need reference genome file (--genome_seq) to write sequence of given regions to fasta file. Exit now.")
	
