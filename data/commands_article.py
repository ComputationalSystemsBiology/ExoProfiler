#!/usr/bin/env python

bed_IMR90 = './inputs/bed/IMR90_GR_chip-seq_rep1_peaks.bed'
bed_K562 = './inputs/bed/K562_GR_chip-seq_rep1_peaks.bed'
bed_U2OS = './inputs/bed/U2OS_GR_chip-seq_peaks.bed'
# bed_CTCF = './inputs/bed/GSM325897_HeLa-CTCF-bsites.txt.hg19'

bam_IMR90 = './inputs/bam/IMR90_GR_chip-exo.bam'
bam_K562 = './inputs/bam/K562_GR_chip-exo.bam'
bam_U2OS = './inputs/bam/U2OS_GR_chip-exo.bam'
# bam_CTCF = './inputs/bam/SRR346403.m1_v3.bam.merged.bam.sorted.bam'

motif_GR = './inputs/motif/MA0113.2.tf'
motif_STAT = './inputs/motif/MA0137.3.tf'
motif_combi = './inputs/motif/oligos_6nt_mkv4_m5.tf'
motif_FOX = './inputs/motif/MA0148.3.tf'
motif_GATA1 = './inputs/motif/MA0035.3.tf'
# motif_CTCF = './inputs/motif/MA0139.1.tf'

genome_hg19 = './inputs/genome/hg19.fa'

def_logsdir = './logs'
def_outputdir = './outputs'
def_figuresdir = './figures'

tool_matrixScanWS = '../python/matrixScanWS.py'
tool_5PrimeCounter = '../python/5PrimeCounter.py'
tool_exoPlotter = '../R/exoPlotter.R'

import os
import shutil
import subprocess


if os.path.exists(def_logsdir):
    shutil.rmtree(def_logsdir, ignore_errors=True)
os.makedirs(def_logsdir)

if os.path.exists(def_outputdir):
    shutil.rmtree(def_outputdir, ignore_errors=True)

if os.path.exists(def_figuresdir):
    shutil.rmtree(def_figuresdir, ignore_errors=True)


##########################################################################################
# Processing

## BAM : IMR90, motif: GR (JASPAR MA0113.2)

# Preprocessing : call matrixScanWS

output_prepdir_IMR90_GR = os.path.join(def_outputdir, 'preprocessing_IMR90_GR')
os.makedirs(output_prepdir_IMR90_GR)

output_preprocessing_IMR90_GR = os.path.join(output_prepdir_IMR90_GR, 'output_matrix-scan_IMR90_GR-MA0113-2_0-0001.txt')
log_stdout = open(os.path.join(def_logsdir, "log_prep_IMR90_GR.txt"), "wb")
subprocess.check_call(['python', tool_matrixScanWS, 
                       '--bed_file', bed_IMR90,
                       '--genome', 'hg19',
                       '--matrix_file', motif_GR,
                       '--uth_pval', '0.0001',
                       '--output_file', output_preprocessing_IMR90_GR,
                       '--perm'
                       ], stdout=log_stdout)
log_stdout.close()

# Compute profiles

output_profdir_IMR90_GR = os.path.join(def_outputdir, 'profile_IMR90_GR')
os.makedirs(output_profdir_IMR90_GR)

output_prof_IMR90_GR = os.path.join(output_profdir_IMR90_GR, 'output_5PrimeCounter_IMR90_GR-MA0113-2_0-0001')
log_stdout = open(os.path.join(def_logsdir, 'log_prof_IMR90_GR.txt'), "wb")
subprocess.check_call(['python', tool_5PrimeCounter,
                       '--input_sites', output_preprocessing_IMR90_GR,
                       '--bam_file', bam_IMR90,
                       '--output_prefix', output_prof_IMR90_GR,
                       '--genome_seq', genome_hg19,
                       '--perm'
                       ], stdout=log_stdout)
log_stdout.close()

# Plot results
log_stdout = open(os.path.join(def_logsdir, 'log_plot_IMR90_GR.txt'), "wb")
subprocess.check_call(["Rscript", tool_exoPlotter, 
                       output_prof_IMR90_GR, 
                       'genome_seq', 
                       'perm'
                       ], stdout=log_stdout)
log_stdout.close()

##########################################################################################
# Processing

## BAM : IMR90, motif: STAT (JASPAR MA0137.3)

# Preprocessing : call matrixScanWS

output_prepdir_IMR90_STAT = os.path.join(def_outputdir, 'preprocessing_IMR90_STAT')
os.makedirs(output_prepdir_IMR90_STAT)

output_preprocessing_IMR90_STAT = os.path.join(output_prepdir_IMR90_STAT, 'output_matrix-scan_IMR90_STAT-MA0137-3_0-0001.txt')
log_stdout = open(os.path.join(def_logsdir, 'log_prep_IMR90_STAT.txt'), "wb")
subprocess.check_call(['python', tool_matrixScanWS, 
                       '--bed_file', bed_IMR90,
                       '--genome', 'hg19',
                       '--matrix_file', motif_STAT,
                       '--uth_pval', '0.0001',
                       '--output_file', output_preprocessing_IMR90_STAT,
                       '--perm'
                       ], stdout=log_stdout)
log_stdout.close()

# Compute profiles

output_profdir_IMR90_STAT = os.path.join(def_outputdir, 'profile_IMR90_STAT')
os.makedirs(output_profdir_IMR90_STAT)

output_prof_IMR90_STAT = os.path.join(output_profdir_IMR90_STAT, 'output_5PrimeCounter_IMR90_STAT-MA0137-3_0-0001')
log_stdout = open(os.path.join(def_logsdir, 'log_prof_IMR90_STAT.txt'), "wb")
subprocess.check_call(['python', tool_5PrimeCounter,
                       '--input_sites', output_preprocessing_IMR90_STAT,
                       '--bam_file', bam_IMR90,
                       '--output_prefix', output_prof_IMR90_STAT,
                       '--genome_seq', genome_hg19,
                       '--perm'
                       ], stdout=log_stdout)
log_stdout.close()

# Plot results

log_stdout = open(os.path.join(def_logsdir, 'log_plot_IMR90_STAT.txt'), "wb")
subprocess.check_call(["Rscript", tool_exoPlotter, 
                       output_prof_IMR90_STAT, 
                       'genome_seq', 
                       'perm'
                       ], stdout=log_stdout)
log_stdout.close()

##########################################################################################
# Processing

## BAM : IMR90, motif: Combi (oligos_6nt_mkv4_m5)

# Preprocessing : call matrixScanWS

output_prepdir_IMR90_Combi = os.path.join(def_outputdir, 'preprocessing_IMR90_Combi')
os.makedirs(output_prepdir_IMR90_Combi)

output_preprocessing_IMR90_Combi = os.path.join(output_prepdir_IMR90_Combi, 'output_matrix-scan_IMR90_Combi-oligos-6nt-mkv4-m5_0-0001.txt')
log_stdout = open(os.path.join(def_logsdir, 'log_prep_IMR90_Combi.txt'), "wb")
subprocess.check_call(['python', tool_matrixScanWS, 
                       '--bed_file', bed_IMR90,
                       '--genome', 'hg19',
                       '--matrix_file', motif_combi,
                       '--uth_pval', '0.0001',
                       '--output_file', output_preprocessing_IMR90_Combi,
                       '--perm'
                       ], stdout=log_stdout)
log_stdout.close()

# Compute profiles

output_profdir_IMR90_Combi = os.path.join(def_outputdir, 'profile_IMR90_Combi')
os.makedirs(output_profdir_IMR90_Combi)

output_prof_IMR90_Combi = os.path.join(output_profdir_IMR90_Combi, 'output_5PrimeCounter_IMR90_Combi-oligos-6nt-mkv4-m5_0-0001')
log_stdout = open(os.path.join(def_logsdir, 'log_prof_IMR90_Combi.txt'), "wb")
subprocess.check_call(['python', tool_5PrimeCounter,
                       '--input_sites', output_preprocessing_IMR90_Combi,
                       '--bam_file', bam_IMR90,
                       '--output_prefix', output_prof_IMR90_Combi,
                       '--genome_seq', genome_hg19,
                       '--perm'
                       ], stdout=log_stdout)
log_stdout.close()

# Plot results
log_stdout = open(os.path.join(def_logsdir, 'log_plot_IMR90_Combi.txt'), "wb")
subprocess.check_call(["Rscript", tool_exoPlotter, 
                       output_prof_IMR90_Combi, 
                       'genome_seq', 
                       'perm'
                       ], stdout=log_stdout)
log_stdout.close()

##########################################################################################
# Processing

## BAM : IMR90, motif: FOX (JASPAR MA0148.3)

# Preprocessing : call matrixScanWS

output_prepdir_IMR90_FOX = os.path.join(def_outputdir, 'preprocessing_IMR90_FOX')
os.makedirs(output_prepdir_IMR90_FOX)

output_preprocessing_IMR90_FOX = os.path.join(output_prepdir_IMR90_FOX, 'output_matrix-scan_IMR90_FOX-MA0148-3_0-0001.txt')
log_stdout = open(os.path.join(def_logsdir, 'log_prep_IMR90_FOX.txt'), "wb")
subprocess.check_call(['python', tool_matrixScanWS, 
                       '--bed_file', bed_IMR90,
                       '--genome', 'hg19',
                       '--matrix_file', motif_FOX,
                       '--uth_pval', '0.0001',
                       '--output_file', output_preprocessing_IMR90_FOX,
                       '--perm'
                       ], stdout=log_stdout)
log_stdout.close()

# Compute profiles

output_profdir_IMR90_FOX = os.path.join(def_outputdir, 'profile_IMR90_FOX')
os.makedirs(output_profdir_IMR90_FOX)

output_prof_IMR90_FOX = os.path.join(output_profdir_IMR90_FOX, 'output_5PrimeCounter_IMR90_FOX-MA0148-3_0-0001')
log_stdout = open(os.path.join(def_logsdir, 'log_prof_IMR90_FOX.txt'), "wb")
subprocess.check_call(['python', tool_5PrimeCounter,
                       '--input_sites', output_preprocessing_IMR90_FOX,
                       '--bam_file', bam_IMR90,
                       '--output_prefix', output_prof_IMR90_FOX,
                       '--genome_seq', genome_hg19,
                       '--perm'
                       ], stdout=log_stdout)
log_stdout.close()

# Plot results
log_stdout = open(os.path.join(def_logsdir, 'log_plot_IMR90_FOX.txt'), "wb")
subprocess.check_call(["Rscript", tool_exoPlotter, 
                       output_prof_IMR90_FOX, 
                       'genome_seq', 
                       'perm'
                       ], stdout=log_stdout)
log_stdout.close()


##########################################################################################
# Processing


## BAM : K562, motif: GR (JASPAR MA0113.2)

# Preprocessing : call matrixScanWS

output_prepdir_K562_GR = os.path.join(def_outputdir, 'preprocessing_K562_GR')
os.makedirs(output_prepdir_K562_GR)

output_preprocessing_K562_GR = os.path.join(output_prepdir_K562_GR, 'output_matrix-scan_K562_GR-MA0113-2_0-0001.txt')
log_stdout = open(os.path.join(def_logsdir, "log_prep_K562_GR.txt"), "wb")
subprocess.check_call(['python', tool_matrixScanWS, 
                       '--bed_file', bed_K562,
                       '--genome', 'hg19',
                       '--matrix_file', motif_GR,
                       '--uth_pval', '0.0001',
                       '--output_file', output_preprocessing_K562_GR,
                       '--perm'
                       ], stdout=log_stdout)
log_stdout.close()

# Compute profiles

output_profdir_K562_GR = os.path.join(def_outputdir, 'profile_K562_GR')
os.makedirs(output_profdir_K562_GR)

output_prof_K562_GR = os.path.join(output_profdir_K562_GR, 'output_5PrimeCounter_K562_GR-MA0113-2_0-0001')
log_stdout = open(os.path.join(def_logsdir, 'log_prof_K562_GR.txt'), "wb")
subprocess.check_call(['python', tool_5PrimeCounter,
                       '--input_sites', output_preprocessing_K562_GR,
                       '--bam_file', bam_K562,
                       '--output_prefix', output_prof_K562_GR,
                       '--genome_seq', genome_hg19,
                       '--perm'
                       ], stdout=log_stdout)
log_stdout.close()

# Plot results
log_stdout = open(os.path.join(def_logsdir, 'log_plot_K562_GR.txt'), "wb")
subprocess.check_call(["Rscript", tool_exoPlotter, 
                       output_prof_K562_GR, 
                       'genome_seq', 
                       'perm'
                       ], stdout=log_stdout)
log_stdout.close()

##########################################################################################
# Processing


## BAM : U2OS, motif: GR (JASPAR MA0113.2)

# Preprocessing : call matrixScanWS

output_prepdir_U2OS_GR = os.path.join(def_outputdir, 'preprocessing_U2OS_GR')
os.makedirs(output_prepdir_U2OS_GR)

output_preprocessing_U2OS_GR = os.path.join(output_prepdir_U2OS_GR, 'output_matrix-scan_GR-MA0113-2_0-0001.txt')
log_stdout = open(os.path.join(def_logsdir, "log_prep_U2OS_GR.txt"), "wb")
subprocess.check_call(['python', tool_matrixScanWS, 
                       '--bed_file', bed_U2OS,
                       '--genome', 'hg19',
                       '--matrix_file', motif_GR,
                       '--uth_pval', '0.0001',
                       '--output_file', output_preprocessing_U2OS_GR,
                       '--perm'
                       ], stdout=log_stdout)
log_stdout.close()

# Compute profiles

output_profdir_U2OS_GR = os.path.join(def_outputdir, 'profile_U2OS_GR')
os.makedirs(output_profdir_U2OS_GR)

output_prof_U2OS_GR = os.path.join(output_profdir_U2OS_GR, 'output_5PrimeCounter_U2OS_GR-MA0113-2_0-0001')
log_stdout = open(os.path.join(def_logsdir, 'log_prof_U2OS_GR.txt'), "wb")
subprocess.check_call(['python', tool_5PrimeCounter,
                       '--input_sites', output_preprocessing_U2OS_GR,
                       '--bam_file', bam_U2OS,
                       '--output_prefix', output_prof_U2OS_GR,
                       '--genome_seq', genome_hg19,
                       '--perm'
                       ], stdout=log_stdout)
log_stdout.close()

# Plot results
log_stdout = open(os.path.join(def_logsdir, 'log_plot_U2OS_GR.txt'), "wb")
subprocess.check_call(["Rscript", tool_exoPlotter, 
                       output_prof_U2OS_GR, 
                       'genome_seq', 
                       'perm'
                       ], stdout=log_stdout)
log_stdout.close()

# ##########################################################################################
# # Processing
# 
# 
# ## BAM : CTCF, motif: GR (JASPAR MA0113.2)
# 
# # Preprocessing : call matrixScanWS
# 
# output_prepdir_CTCF_GR = os.path.join(def_outputdir, 'preprocessing_CTCF_GR')
# os.makedirs(output_prepdir_CTCF_GR)
# 
# output_preprocessing_CTCF_GR = os.path.join(output_prepdir_CTCF_GR, 'output_matrix-scan_CTCF_GR-MA0113-2_0-0001.txt')
# log_stdout = open(os.path.join(def_logsdir, "log_prep_CTCF_GR.txt"), "wb")
# subprocess.check_call(['python', tool_matrixScanWS, 
#                        '--bed_file', bed_CTCF,
#                        '--genome', 'hg19',
#                        '--matrix_file', motif_GR,
#                        '--uth_pval', '0.0001',
#                        '--output_file', output_preprocessing_CTCF_GR,
#                        '--perm'
#                        ], stdout=log_stdout)
# log_stdout.close()
# 
# # Compute profiles
# 
# output_profdir_CTCF_GR = os.path.join(def_outputdir, 'profile_CTCF_GR')
# os.makedirs(output_profdir_CTCF_GR)
# 
# output_prof_CTCF_GR = os.path.join(output_profdir_CTCF_GR, 'output_5PrimeCounter_CTCF_GR-MA0113-2_0-0001')
# log_stdout = open(os.path.join(def_logsdir, 'log_prof_CTCF_GR.txt'), "wb")
# subprocess.check_call(['python', tool_5PrimeCounter,
#                        '--input_sites', output_preprocessing_CTCF_GR,
#                        '--bam_file', bam_CTCF,
#                        '--output_prefix', output_prof_CTCF_GR,
#                        '--genome_seq', genome_hg19,
#                        '--perm'
#                        ], stdout=log_stdout)
# log_stdout.close()
# 
# # Plot results
# log_stdout = open(os.path.join(def_logsdir, 'log_plot_CTCF_GR.txt'), "wb")
# subprocess.check_call(["Rscript", tool_exoPlotter, 
#                        output_prof_CTCF_GR, 
#                        'genome_seq', 
#                        'perm'
#                        ], stdout=log_stdout)
# log_stdout.close()

##########################################################################################
# Processing


## BAM : K562, motif: GATA1 (JASPAR MA0035.3)

# Preprocessing : call matrixScanWS

output_prepdir_K562_GATA1 = os.path.join(def_outputdir, 'preprocessing_K562_GATA1')
os.makedirs(output_prepdir_K562_GATA1)

output_preprocessing_K562_GATA1 = os.path.join(output_prepdir_K562_GATA1, 'output_matrix-scan_K562_GATA1-MA0035-3_0-0001.txt')
log_stdout = open(os.path.join(def_logsdir, 'log_prep_K562_GATA1.txt'), "wb")
subprocess.check_call(['python', tool_matrixScanWS, 
                       '--bed_file', bed_K562,
                       '--genome', 'hg19',
                       '--matrix_file', motif_GATA1,
                       '--uth_pval', '0.0001',
                       '--output_file', output_preprocessing_K562_GATA1,
                       '--perm'
                       ], stdout=log_stdout)
log_stdout.close()

# Compute profiles

output_profdir_K562_GATA1 = os.path.join(def_outputdir, 'profile_K562_GATA1')
os.makedirs(output_profdir_K562_GATA1)

output_prof_K562_GATA1 = os.path.join(output_profdir_K562_GATA1, 'output_5PrimeCounter_K562_GATA1-MA0035-3_0-0001')
log_stdout = open(os.path.join(def_logsdir, 'log_prof_K562_GATA1.txt'), "wb")
subprocess.check_call(['python', tool_5PrimeCounter,
                       '--input_sites', output_preprocessing_K562_GATA1,
                       '--bam_file', bam_K562,
                       '--output_prefix', output_prof_K562_GATA1,
                       '--genome_seq', genome_hg19,
                       '--perm'
                       ], stdout=log_stdout)
log_stdout.close()

# Plot results

log_stdout = open(os.path.join(def_logsdir, 'log_plot_K562_GATA1.txt'), "wb")
subprocess.check_call(["Rscript", tool_exoPlotter, 
                       output_prof_K562_GATA1, 
                       'genome_seq', 
                       'perm'
                       ], stdout=log_stdout)
log_stdout.close()


# ##########################################################################################
# # Processing
# 
# 
# ## BAM : CTCF, motif: CTCF (JASPAR MA0139.1)
# 
# # Preprocessing : call matrixScanWS
# 
# output_prepdir_CTCF_CTCF = os.path.join(def_outputdir, 'preprocessing_CTCF_CTCF')
# os.makedirs(output_prepdir_CTCF_CTCF)
# 
# output_preprocessing_CTCF_CTCF = os.path.join(output_prepdir_CTCF_CTCF, 'output_matrix-scan_CTCF_CTCF-MA0139-1_0-0001.txt')
# log_stdout = open(os.path.join(def_logsdir, "log_prep_CTCF_CTCF.txt"), "wb")
# subprocess.check_call(['python', tool_matrixScanWS, 
#                        '--bed_file', bed_CTCF,
#                        '--genome', 'hg19',
#                        '--matrix_file', motif_CTCF,
#                        '--uth_pval', '0.0001',
#                        '--output_file', output_preprocessing_CTCF_CTCF,
#                        '--perm'
#                        ], stdout=log_stdout)
# log_stdout.close()
# 
# # Compute profiles
# 
# output_profdir_CTCF_CTCF = os.path.join(def_outputdir, 'profile_CTCF_CTCF')
# os.makedirs(output_profdir_CTCF_CTCF)
# 
# output_prof_CTCF_CTCF = os.path.join(output_profdir_CTCF_CTCF, 'output_5PrimeCounter_CTCF_CTCF-MA0139-1_0-0001')
# log_stdout = open(os.path.join(def_logsdir, 'log_prof_CTCF_CTCF.txt'), "wb")
# subprocess.check_call(['python', tool_5PrimeCounter,
#                        '--input_sites', output_preprocessing_CTCF_CTCF,
#                        '--bam_file', bam_CTCF,
#                        '--output_prefix', output_prof_CTCF_CTCF,
#                        '--genome_seq', genome_hg19,
#                        '--size', '100',
#                        '--perm'
#                        ], stdout=log_stdout)
# log_stdout.close()
# 
# # Plot results
# 
# log_stdout = open(os.path.join(def_logsdir, 'log_plot_CTCF_CTCF.txt'), "wb")
# subprocess.check_call(["Rscript", tool_exoPlotter, 
#                        output_prof_CTCF_CTCF, 
#                        'genome_seq', 
#                        'perm'
#                        ], stdout=log_stdout)
# log_stdout.close()


##########################################################################################
# Figure 1
            

figures_dir = os.path.join(def_figuresdir, 'figure_1')
os.makedirs(figures_dir)

#3 images : QC , heatmap, profile w. perm
shutil.copy(output_prof_IMR90_GR+'_seq.pdf', figures_dir)
shutil.copy(output_prof_IMR90_GR+'.maxstrand_heatmap.pdf', figures_dir)
shutil.copy(output_prof_IMR90_GR+'.profile-perm.pdf', figures_dir)



##########################################################################################
# Figure 2


figures_dir = os.path.join(def_figuresdir, 'figure_2')
os.makedirs(figures_dir)

shutil.copy(output_prof_IMR90_GR+'.profile.pdf', figures_dir)



##########################################################################################
# Figure 4


### Motif STAT (A)

figures_dir = os.path.join(def_figuresdir, 'figure_4A')
os.makedirs(figures_dir)

shutil.copy(output_prof_IMR90_STAT+'.profile-perm.pdf', figures_dir)

### Motif Combi (B)

figures_dir = os.path.join(def_figuresdir, 'figure_4B')
os.makedirs(figures_dir)

shutil.copy(output_prof_IMR90_Combi+'.profile-perm.pdf', figures_dir)

### Motif FOX (C)

figures_dir = os.path.join(def_figuresdir, 'figure_4C')
os.makedirs(figures_dir)

shutil.copy(output_prof_IMR90_FOX+'.profile-perm.pdf', figures_dir)

##########################################################################################
# Figure 6


figures_dir = os.path.join(def_figuresdir, 'figure_6')
os.makedirs(figures_dir)

shutil.copy(output_prof_IMR90_FOX+'.cluster.pdf', figures_dir)

##########################################################################################
# Figure S3


##########################################################################################
# Figure S2

figures_dir = os.path.join(def_figuresdir, 'figure_S2')
os.makedirs(figures_dir)

shutil.copy(output_prof_IMR90_GR+'.maxstrand_heatmap.pdf', figures_dir)
shutil.copy(output_prof_IMR90_GR+'.profile-perm.pdf', figures_dir)
shutil.copy(output_prof_IMR90_GR+'.permut_matrix_10.values.tab', figures_dir)

shutil.copy(output_prof_K562_GR+'.maxstrand_heatmap.pdf', figures_dir)
shutil.copy(output_prof_K562_GR+'.profile-perm.pdf', figures_dir)
shutil.copy(output_prof_K562_GR+'.permut_matrix_10.values.tab', figures_dir)

shutil.copy(output_prof_U2OS_GR+'.maxstrand_heatmap.pdf', figures_dir)
shutil.copy(output_prof_U2OS_GR+'.profile-perm.pdf', figures_dir)
shutil.copy(output_prof_U2OS_GR+'.permut_matrix_10.values.tab', figures_dir)

# shutil.copy(output_prof_CTCF_GR+'.maxstrand_heatmap.pdf', figures_dir)
# shutil.copy(output_prof_CTCF_GR+'.profile-perm.pdf', figures_dir)
# shutil.copy(output_prof_CTCF_GR+'.permut_matrix_10.values.tab', figures_dir)

##########################################################################################
# Figure S4

figures_dir = os.path.join(def_figuresdir, 'figure_S4')
os.makedirs(figures_dir)

shutil.copy(output_prof_K562_GATA1+'.profile-perm.pdf', figures_dir)
shutil.copy(output_prof_K562_GATA1+'.permut_matrix_10.values.tab', figures_dir)

# ##########################################################################################
# # Figure S7
# 
# figures_dir = os.path.join(def_figuresdir, 'figure_S7')
# os.makedirs(figures_dir)
# 
# shutil.copy(output_prof_CTCF_CTCF+'.profile-perm.pdf', figures_dir)
# shutil.copy(output_prof_CTCF_CTCF+'.permut_matrix_10.values.tab', figures_dir)

