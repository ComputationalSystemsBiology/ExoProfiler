#!/usr/bin/python

###########################################################'
# 
# ExoScanner
# Client to run the RSAT tool 'matrix-scan'
#
#
#
# Steps and RSAT tools involved (as web services) : 
#
# 1/ conversion of the input bed file (option -b or --bed_file) into a fasta file, 
# using a reference genome (option -g or --genome). This step uses 'fetch_sequences'.
# 2/ scan the FASTA sequences for motifs, provided as a transfac matrix (option -m
# or --matrix_file). Only matches with a p-value below a threshold are reported (option
# -u or --uth_pval). This step uses 'matrix_scan'.
# 3/ Matches are written in an output file (option -o or --output_file).
#
#
# Validation :
#
# If the --perm (or -p) option is included in the command, additional steps will 
# be performed.
# The transfac matrix will be permuted 10 times using the tool 'convert_matrix'. Those 
# permuted matrices should not be too similar to the transfac matrix. This is checked 
# using the tool 'compare_matrices'. They should not be duplicated either. This is 
# checked again using 'compare_matrices'.
# The permuted matrices will then be used by 'matrix_scan' to scan the FASTA sequences 
# (computed instep 1/), using the same upper threshold for the p-value as in step 2/.
# Finally, permuted matrices and outputs from 'matrix_scan' are written in the same 
# folder as the output file of the input transfac matrix.  
#
#
# usage: matrixScanWS.py [-h] -b <BED_FILE> -g <GENOME_NAME> -m
#                             <MATRIX_FILE> -o <OUTPUT_FILE> -u <P-VALUE> [-p]
# 
# 
# optional arguments:
#   -h, --help            show this help message and exit
#
#   -b <BED_FILE>, --bed_file <BED_FILE>
#                         File in BED format.
#   -g <GENOME_NAME>, --genome <GENOME_NAME>
#                         Genome on which BED file is to be matched. Example: mm9, hg19.
#
#   -m <MATRIX_FILE>, --matrix_file <MATRIX_FILE>
#                         Matrix in transfac format.
#   -u <P-VALUE>, --uth_pval <P-VALUE>
#                         Upper threshold for p-value. Should be in [0.0-1.0].
#
#   -o <OUTPUT_FILE>, --output_file <OUTPUT_FILE>
#                         Output file name for matrix_scan result on the input matrix. 
#
#   -p, --perm            (Optional) Flag indicating that permuted matrices
#                         should also be computed and used for scanning. Output file names from permuted 
#                         matrices will contain the tag '_perm' and an index. The files will be created in
#                         the same folder as the output file.
# 
# Last update: 26.09.2014 (CEH)
# 
#
###########################################################'


matrixscanClientVersion = '0.4 - 30/07/2014'


## Tested with python 2.7.5 and suds version 0.4.1
## Other imports:
## - logging
## - os
## - os.path
## - platform
## - urllib2
## - argparse
## - warnings
## - re
## - numpy


# Help on arguments
# python matrixScanWS.py -h
# python matrixScanWS.py --help


# Command line with argument
# python matrixScanWS.py -b <BED_FILE> -g <GENOME_NAME> -m <MATRIX_FILE> -o <OUTPUT_FILE> -u <P-VALUE> [-p]
# python matrixScanWS.py --bed_file <BED_FILE> --genome <GENOME_NAME> --matrix_file <MATRIX_FILE> --output_file <OUTPUT_FILE> --uth_pval <P-VALUE> [--perm]

# python matrixScanWS.py -b fetch-sequences_Schmidt_2011_mm9_CEBPA_SWEMBL_R0.12_702peaks.bed -g mm9 -m param_matrix.txt -u 0.001 -o ./output.txt -p

#####################################################################################################'

## SET UP



###########################################################'
## Define log options for suds

# Import log package
import logging
# Configure default log package to INFO
logging.basicConfig(level=logging.INFO)
# Configure log of suds clients to DEBUG for verbose output concerning Client request
logging.getLogger('suds.client').setLevel(logging.ERROR)




###########################################################'
## Create the SOAP client to request RSAT services


# Load Client class from suds
from suds.client import Client
# Define URL for RSAT services 
#wsdlUrl =  'http://rsat.ulb.ac.be/rsat/web_services/RSATWS.wsdl'
rsatURL = 'http://pedagogix-tagc.univ-mrs.fr/rsat/'
wsdlUrl =  rsatURL + '/web_services/RSATWS.wsdl'
# Create the client
client = Client(wsdlUrl)
# Need the service interface to perform requests
rsat_service = client.service

# Define client header (optional but convenient)
import os, platform
userAgent = 'RSAT-Client/v%s (%s; Python %s; %s)' % (
    matrixscanClientVersion, 
    os.path.basename( __file__ ),
    platform.python_version(), 
    platform.system()
)
httpHeaders = {'User-agent': userAgent}
client.set_options(headers=httpHeaders)

# Define time out limit
client.set_options(timeout=300)

###########################################################'
## Define functions to make a service perform the desired request using provided arguments

# Transform bed files into fasta with Galaxy header
def call_fetch_sequences(service, bed_content, genome_4bed):

	# Wrap all arguments into a named list
	#http://rsat.ulb.ac.be/web_services/RSATWS_documentation.xml
	arguments_fetch = {
		'input' : bed_content,
		'genome' : genome_4bed,
		'header_format' : 'galaxy',
		'output' : 'client'
	}

	# Perform SOAP request on RSAT server
	result = service.fetch_sequences(arguments_fetch)
	
	# Return only the desired result
	return result.client


# Create permuted matrices if user provided option --perm (or -p) in the command line
def call_convert_matrix(service, matrix_tf):

	# Wrap all arguments into a named list (dictionary), including some hard-coded characters
	#http://rsat.ulb.ac.be/web_services/RSATWS_documentation.xml
	arguments_conv = {
		'matrix' : matrix_tf,
		'from' : 'transfac',
		'to' : 'transfac',
		'perm' : 1,
		'output' : 'client'
	}

	## Perform SOAP request on RSAT server
	result = service.convert_matrix(arguments_conv)
	
	return result.client


# Call compare_matrices Web Service
def perform_compare_matrices(service, args):

	# Start the job on the RSAT server...
	result = service.compare_matrices(args)

	# But the server starts it in background. We have to access the output file 
	# (extracted from the executed command) and download/read the file ourselves
	command = result.command
	url_to_download = rsatURL + command[command.rindex('-o \'$RSAT/public_html/')+21:-1]

	# Open URL of the output file
	import urllib2, numpy
	# 	response = urllib2.urlopen(url_to_download)
	#  	print response.read()

	# Read content of output
	b = numpy.genfromtxt(urllib2.urlopen(url_to_download), delimiter='\t', names=True, dtype=None)
	# Values in b['m2']
	#print b['m2']

	return b

# Compare_matrices : Case where we have two sets of matrices
def call_compare_matrices_2sets(service, matrix1, matrix2, Ncor_max):

	# Wrap all arguments into a named list
	#http://rsat.ulb.ac.be/web_services/RSATWS_documentation.xml
	arguments_comp = {
		'matrix_1' : matrix1,
		'format1' : 'transfac',
		'matrix_2' : matrix2,
		'format2' : 'transfac',
		'strand' : 'DR',
		'lth' : ['Ncor '+str(Ncor_max)],
		'return' : 'matrix_number,Ncor',
		'output' : 'both'
	}
	return perform_compare_matrices(service, arguments_comp)

# Compare_matrices : Case where we have one set of matrices
def call_compare_matrices_1set(service, matrices, Ncor_max):
	
	# Wrap all arguments into a named list
	#http://rsat.ulb.ac.be/web_services/RSATWS_documentation.xml
	arguments_comp = {
		'matrix' : matrices,
		'format' : 'transfac',
		'distinct' : 1,
		'lth' : ['Ncor '+str(Ncor_max)],
		'return' : 'matrix_number,Ncor',
		'output' : 'both'
	}
	return perform_compare_matrices(service, arguments_comp)


# Scan sequences for a given matrix
def call_matrix_scan(service, fasta_content_str, matrix_str, uth_val):

	# Wrap all arguments into a named list
	#http://rsat.ulb.ac.be/web_services/RSATWS_documentation.xml
	arguments = {
		'sequence' : fasta_content_str,
		'matrix' : matrix_str,
		'uth' : ['pval '+str(uth_val)],
 		'quick' : 1,
		'str' : 2,
		'origin' : 'start',
	
		'background_input' : 1, # this option requires 'markov'
		'markov' : 1,
	
		'pseudo' : 1,
		'n_treatment' : 'score',
		'matrix_format' : 'transfac',
		'background_pseudo' : 0.01,
# 		'verbosity' : 1
	}

	## Perform SOAP request on RSAT server
	result = service.matrix_scan(arguments)
	return result.client


#####################################################################################################'

## MAIN EXECUTION STEPS



#####################################################################################################'
## Read and process command-line


# See https://docs.python.org/2/library/argparse.html
import argparse

# Create the parser
parser = argparse.ArgumentParser(description='Client to run the RSAT tool \'matrix-scan\'', epilog='Version '+matrixscanClientVersion)

# List possible arguments
parser.add_argument('-b', '--bed_file', metavar='<BED_FILE>', type=argparse.FileType('r'), nargs=1, help='File in BED format.', required=True)
parser.add_argument('-g', '--genome', metavar='<GENOME_NAME>', type=str, nargs=1, help='Genome on which BED file is to be matched.', required=True)
parser.add_argument('-m', '--matrix_file', metavar='<MATRIX_FILE>', type=argparse.FileType('r'), nargs=1, help='Matrix in transfac format.', required=True)
# to force uth_pvalue to be between 0 and 1 (floating values)
def restricted_float(x):
    x = float(x)
    if x < 0.0 or x > 1.0:
        raise argparse.ArgumentTypeError("%r not in range [0.0, 1.0]" % (x,))
    return x
parser.add_argument('-u', '--uth_pval', metavar='<P-VALUE>', type=restricted_float, nargs=1, help='Upper threshold for p-value. Should be in [0.0-1.0].', required=True)
parser.add_argument('-o', '--output_file', metavar='<OUTPUT_FILE>', type=argparse.FileType('w', 0), nargs=1, help='Output file name for matrix_scan result on the input matrix.', required=True)
parser.add_argument('-p', '--perm', action='store_true', help='(Optional) Flag indicating that permuted matrices should also be computed and used for scanning. Output file names from permuted matrices will contain the tag \'_perm\' and an index. The files will be created in the same folder as the output file.', required=False)

# Parse command line
args = parser.parse_args()


#########################################################'
## Prepare sequences

# Read content of the bed file
bed_filecontent = args.bed_file[0].read()

# Transform BED to FASTA using service: fetch_sequence
fasta_content = call_fetch_sequences(rsat_service, bed_filecontent, args.genome[0])



#########################################################'
## Prepare matrix

# Input matrix (file content)
matrix = args.matrix_file[0].read()

# Input matrix (file name)
matrix_filename = args.matrix_file[0].name


#########################################################'
## Perform SOAP request on RSAT server with the matrix and the FASTA sequences

# Call web service matrix_scan
result = call_matrix_scan(rsat_service, fasta_content, matrix, args.uth_pval[0])

# Output file for this matrix, from arguments in the command line
output_ffile = args.output_file[0]
# Output file name
output_filename = args.output_file[0].name

# Write result in output file
output_ffile.write(result)
output_ffile.close()



#####################################################################################################'

## OPTIONAL STEPS



#########################################################'
## If asked by the user, perform same pipeline on a certain number of permuted matrices

# --perm option (or -p) is present?
if(args.perm):

	####################
	# First, we permute X times the input matrix using convert_matrix,
	# checking each time with compare_matrices that the obtained permutation
	# is not too similar to the original matrix

	# container to store the matrices	
	nb_permutations = 10
	permuted_and_NcorOK = []
	
	# correlation threshold below which matrices are considered different enough
	Ncor_threshold = 0.4
	
	# Loop until we get a certain amount of permuted matrices corresponding to our criterium
	current_index = 0
	while current_index < nb_permutations:
	
		# Call convert_matrix to get a permuted matrix
		one_perm_matrix = call_convert_matrix(rsat_service, matrix)

		# Compare with the original transfac matrix
		# We have to test if permuted matrices are not too similar to the original matrix
		# by calling compare_matrices
		# Returns a result with size >0 only if the permuted matrix doesn't fit our criterium
		result = call_compare_matrices_2sets(rsat_service, matrix, one_perm_matrix, Ncor_threshold)

		# Keep permutation only if not too similar (result will not be empty if matrix too similar)
		if result.size == 0:
			permuted_and_NcorOK.append(one_perm_matrix)
			current_index += 1
		
	####################
	
	# Now we have 10 permuted matrices
	# Just check that we didn't get doubles during this process
	result_comp = call_compare_matrices_1set(rsat_service, ''.join(permuted_and_NcorOK), 0.99)

	if result_comp.size > 0:
		import warnings; warnings.warn("Biased control: duplication among permuted matrices!")

	####################

	## Call matrix_scan on each of the permuted matices
	
	# First rename matrices with sequential index

	import re
	# Extract AC and Id of original matrix
	orig_AC = re.search(r'^AC.*$', matrix, re.M).group()
	orig_ID = re.search(r'^ID.*$', matrix, re.M).group()
	# Trim them
	orig_AC = re.sub(r'\s*$', '', orig_AC)
	orig_ID = re.sub(r'\s*$', '', orig_ID)
	# Rename AC and ID of all permuted matrices (sequential indices)
	for index in xrange(nb_permutations):	
		permuted_and_NcorOK[index] = re.sub(r'AC.*\n', orig_AC+'_perm'+str(index+1)+'\n', permuted_and_NcorOK[index], 1)
		permuted_and_NcorOK[index] = re.sub(r'ID.*\n', orig_ID+'_perm'+str(index+1)+'\n', permuted_and_NcorOK[index], 1)
		
	
	# Create/write files containing each permuted matrix, in folder where output has to be written

	import os.path
	# Get output folder provided by user
	output_folder = os.path.dirname(output_filename)
	# Split file name and extension (needed to create filenames for permuted matrices outputs)
	(output_filename_label, output_filename_ext) = os.path.splitext(output_filename)

	# Extract name of matrix file
	matrix_filename_name = os.path.split(matrix_filename)[1] # discard file path
	# Split between real name and extension (matrix)
	(matrix_filename_label, matrix_filename_ext) = os.path.splitext(matrix_filename_name)

	for index in xrange(nb_permutations):
		# Write each permuted matrix
		file_perm = open(os.path.join(output_folder, matrix_filename_label+'_perm'+str(index+1)+matrix_filename_ext), 'w')
		file_perm.write(permuted_and_NcorOK[index])
		file_perm.close()

		# Call matrix_scan (as RSAT Web Service)
		result_perm = call_matrix_scan(rsat_service, fasta_content, permuted_and_NcorOK[index], args.uth_pval[0])
	
		# Write outputs for each permuted matrix
		fileo_perm = open(output_filename_label+'_perm'+str(index+1)+output_filename_ext, 'w')
		fileo_perm.write(result_perm)
		fileo_perm.close()



#import sys; sys.exit(0)
