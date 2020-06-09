# tutorial.py

import TumorDecon as td

"""
To perform digital cytometry on a mixture data, given cell signatures and the deconvolution method of your choice, run:

 td.tumor_deconvolve(mixture_data, method, patient_IDs='ALL', cell_signatures=None, up_genes=None, down_genes=None)

Explanation of Inputs:
    - mixture_data: pandas df of rna gene expression data.
        Rows are genes (indexed by 'Hugo_Symbol') and columns are patients
    - method: Method for tumor deconvolution. Must be either 'CIBERSORT', 'DeconRNASeq', 'ssGSEA', or 'SingScore'
    - patient_IDs: list of patient IDs to run DeconRNASeq for.
        Alternatively, can use the string 'ALL' to run for all patients
    - cell_signatures: pandas df of Signature gene expression values for given cell types.
        Rows are genes (indexed by 'Hugo_Symbol') and columns are cell types
        *** Required for 'CIBERSORT' and 'DeconRNASeq' methods ***
            (ignored for 'ssGSEA' and 'SingScore' methods)
    - up_genes: dictionary. Keys are cell types, vaues are a list of up-regulated genes (hugo symbols) for
        the given cell type
        *** Required for 'ssGSEA' and 'SingScore' methods ***
            (ignored for 'CIBERSORT' and 'DeconRNASeq' methods)
    - down_genes: dictionary. Keys are cell types, vaues are a list of down-regulated genes (hugo symbols) for
        the given cell type
        *** Optional for 'SingScore' method ***
            (ignored for 'CIBERSORT', 'DeconRNASeq', and 'ssGSEA' methods)
    - args: Dictionary of arguments to pass on to the given method
Output:
    - pandas dataframe of fractions/ranks of each cell type in the cell signatures, for each specified patient in 'patient_IDs' list.
        Rows are indexed by cell type, columns are patient IDs
"""

# Sample data location
data_loc = td.get_td_Home()+"data/"

## read in sample data (download colon cancer gene expressions directly from cbioportal):
rna = td.download_from_cbio(url="http://download.cbioportal.org/coadread_tcga_pan_can_atlas_2018.tar.gz", fetch_missing_hugo=False)

## Drop any rows (genes) that contain a NaN value from rna expression data
rna.dropna(axis=0, inplace=True)

################################################################################
############################# Linear Models ####################################
################################################################################
## Linear models (cibersort and DeconRNASeq) require a signature matrix, and output cell fractions:

## read in LM22 file as signature matrix
sig = td.read_lm22_file(data_loc+'LM22.txt')

## Run cibersort on ALL patients:
## optional argments include:
##    'scaling': how to scale the data (rna mixture & signature matrix) before applying DeconRNASeq (must be either 'None', 'zscore', or 'minmax')
##    'scaling_axis': whether to scale data by normalizing each cell types / patient individually (0) or each gene individually (1). default is 0
##    'nu','C','kernel': all arguments to be passed to sklearn's implementation of NuSVR. See sklearn documentation for details
##       additionally, can pass in 'nu':'best' to choose the nu from {0.25, 0.5, 0.75} that minimizes the root mean square error
##       defaults are nu='best', C=1.0, kernel=linear
##   'print_progress': whether to print patient ID as ssGSEA iterates through (default: False)
ciber_freqs = td.tumor_deconvolve(rna, 'cibersort',  patient_IDs='ALL', cell_signatures=sig, args={'nu':'best', 'scaling':'minmax'})
print(ciber_freqs)

## Run DeconRNASeq:
## optional argments include:
##    'scaling': Same as in cibersort
##    'scaling_axis': Same as in cibersort
##    'check_sig': whether to check the condition number of the signature matrix before solving (default: False)
##    'print_results': whether to print the results of each fit while function is running (default: False)
##   'print_progress': whether to print patient ID as ssGSEA iterates through (default: False)
decon_freqs = td.tumor_deconvolve(rna, 'DeconRNASeq',  patient_IDs='ALL', cell_signatures=sig, args={'scaling':'minmax', 'print_results':False})
print(decon_freqs)

# ##############################################################################
# ######################### Rank-Based Methods #################################
# ##############################################################################
# Rank-based methods (ssGSEA and SingScore) require a list of up-regulated genes

## Read in up-regulated gene set derived from LM22:
up_geneset = td.read_geneset(data_loc+"LM22_up_genes.csv")

####
# When running tumor_deconvolve(), the default is to apply the method on all patients
# in the rna data frame (patient_IDs='ALL'). To instead run the method on only a
# select list of patients, pass in a list of patient identifiers (column names) instead
####
patient_subset = ['TCGA-3L-AA1B-01', 'TCGA-4N-A93T-01', 'TCGA-4T-AA8H-01', 'TCGA-5M-AAT4-01',
                    'TCGA-5M-AAT5-01', 'TCGA-5M-AAT6-01', 'TCGA-5M-AATA-01', 'TCGA-5M-AATE-01',
                    'TCGA-A6-2675-01', 'TCGA-A6-2682-01', 'TCGA-A6-2684-01', 'TCGA-A6-2685-01']

## Run ssGSEA on only a few patients:
## optional arguments include:
##   'alpha': Weight used in ssGSEA method (default: 1.0)
##   'norm': whether to normalize enrichment scores, as done by Barbie et al., 2009, online methods, pg. 2
##   'print_progress': whether to print patient ID as ssGSEA iterates through (default: False)
ssgsea_scores = td.tumor_deconvolve(rna, 'ssGSEA',  patient_IDs=patient_subset, up_genes=up_geneset, args={'alpha':0.5, 'norm':True, 'print_progress':False})
print(ssgsea_scores)

## SingScore can be run with just an up-regulated gene set (unidirectional):
singscore_unidirectional = td.tumor_deconvolve(rna, 'singscore',  patient_IDs='ALL', up_genes=up_geneset)
print(singscore_unidirectional)

## or with both an up-regulated and down-regulated gene set (bidirectional):
## (Can also derive up-regulated and down-regulated genes from a signature matrix, such as LM6 file):
LM6 = td.read_lm22_file(data_loc+'LM6.txt')
up_geneset_LM6, down_geneset_LM6 = td.find_up_down_genes_from_sig(LM6, down_cutoff=0.4, up_cutoff=4.0)

singscore_bidirectional = td.tumor_deconvolve(rna, 'singscore',  patient_IDs=patient_subset, up_genes=up_geneset_LM6, down_genes=down_geneset_LM6)
print(singscore_bidirectional)

################################################################################
######################### Visualize Results ####################################
################################################################################
## WARNING: These functions currently only work for fractions/ranks generated with LM22 cell signatures

# Select a subset of the patients to visualize:
results = ciber_freqs.loc[patient_subset]

## Can visualize results with boxplots:
td.cell_frequency_boxplot(results)

## Can also visualize with barcharts:
td.cell_frequency_barchart(results)

# Clustermaps:
td.hierarchical_clustering(results)

# and Pair plots:
td.pair_plot(results)
