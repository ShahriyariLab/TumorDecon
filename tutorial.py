# tutorial.py

# The next 2 lines only needed if running tutorial from inside the package directory
import sys
sys.path.insert(0,"..")

import TumorDecon as td
td_install_loc = td.get_td_Home()

# read in colon cancer from TCGA Pancancer:
# rna = td.read_pancancer_rna_file(td_install_loc+'data/coadred_data_RNA_Seq_v2_expression_median.txt')

# OR: download from cbioportal:
rna = td.download_from_cbio(url="http://download.cbioportal.org/coadread_tcga_pan_can_atlas_2018.tar.gz")

# Drop any rows (genes) that contain a NaN value from rna expression data
rna.dropna(axis=0, inplace=True)

#################################################################################
######################### Frequency Methods #####################################
#################################################################################
# Frequency methods (cibersort and DeconRNASeq) require a signature matrix:
# read in LM22 file as signature matrix
sig = td.read_lm22_file(td_install_loc+'data/LM22.txt')

# Run cibersort on ALL patients:
# optional argments include:
#    'scaling': how to scale the data (rna mixture & signature matrix) before applying DeconRNASeq (must be either 'None', 'zscore', or 'minmax')
#    'scaling_axis': whether to scale data by normalizing each cell types / patient individually (0) or each gene individually (1). default is 0
#    'nu','C','kernel': all arguments to be passed to sklearn's implementation of NuSVR. See sklearn documentation for details
#       additionally, can pass in 'nu':'best' to choose the nu from {0.25, 0.5, 0.75 that minimizes the root mean square error}
#       defaults are nu='best', C=1.0, kernel=linear
ciber_freqs = td.tumor_deconvolve(rna, 'cibersort',  patient_IDs='ALL', cell_signatures=sig, args={'nu':'best', 'scaling':'minmax', 'scaling_axis':1})
print(ciber_freqs)

# Can visualize results with boxplots:
# (Currently only works for frequencies generated for the 22 cell types in the LM22 signature matrix):
td.cell_frequency_boxplot(ciber_freqs)

# Can also visualize with barcharts:
td.cell_frequency_barchart(ciber_freqs)

# Clustermaps:
td.hierarchical_clustering(ciber_freqs)

# Pair plots:
td.pair_plot(ciber_freqs)

# Run DeconRNASeq:
# optional argments include:
#    'scaling': Same as in cibersort
#    'scaling_axis': Same as in cibersort
#    'check_sig': whether to check the condition number of the signature matrix before solving (default: False)
#    'print_results': whether to print the results of each fit while function is running (default: True)
decon_freqs = td.tumor_deconvolve(rna, 'DeconRNASeq',  patient_IDs='ALL', cell_signatures=sig, args={'scaling':'minmax', 'scaling_axis':0, 'print_results':False})
print(decon_freqs)

# #################################################################################
# ######################### Rank-Based Methods ####################################
# #################################################################################
# Rank-based methods (ssGSEA and SingScore) require a list of up-regulated genes:
gene_set = td.read_ssGSEA_up_genes()

####
# When running tumor_deconvolve(), the default is to apply the method on all patients
# in the rna data frame (patient_IDs='ALL'). To instead run the method on only a
# select list of patients, pass in a list of patient identifiers (column names) instead
####
patient_subset = ['TCGA-3L-AA1B-01', 'TCGA-4N-A93T-01', 'TCGA-4T-AA8H-01', 'TCGA-5M-AAT4-01',
                    'TCGA-5M-AAT5-01', 'TCGA-5M-AAT6-01', 'TCGA-5M-AATA-01', 'TCGA-5M-AATE-01',
                    'TCGA-A6-2675-01', 'TCGA-A6-2682-01', 'TCGA-A6-2684-01', 'TCGA-A6-2685-01']

# Run ssGSEA on only a few patients:
# optional arguments include:
#   'alpha': Weight used in ssGSEA method (default: 1.0)
#   'print_progress': whether to print patient ID as ssGSEA iterates through (default: True)
ssgsea_scores = td.tumor_deconvolve(rna, 'ssGSEA',  patient_IDs=patient_subset, up_genes=gene_set, args={'alpha':0.5, 'print_progress':False})
print(ssgsea_scores)

# SingScore can be run with just an up-regulated gene set (unidirectional):
singscore_unidirectional = td.tumor_deconvolve(rna, 'singscore',  patient_IDs='ALL', up_genes=gene_set)
print(singscore_unidirectional)

# or with both an up-regulated and down-regulated gene set (bidirectional):
# (Can also get up-regulated and down-regulated genes from a signature matrix, such as the LM22 file):
up_set, down_set = td.find_up_down_genes_from_sig(sig, down_cutoff=0.4, up_cutoff=4.0)
singscore_bidirectional = td.tumor_deconvolve(rna, 'singscore',  patient_IDs=patient_subset, up_genes=up_set, down_genes=down_set)
print(singscore_bidirectional)
