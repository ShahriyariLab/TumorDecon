# tutorial.py

"""
If you use this package (or some parts of these codes), please cite:
    T. Le, R. Aronow, A. Kirshtein, L. Shahriyari,
    "A review of digital cytometry methods: estimating the relative abundance of cell types in a bulk of cells",
    Briefings in Bioinformatics, 2020, https://doi.org/10.1093/bib/bbaa219.
"""

"""
## To install TumorDecon with pip, use:

pip install git+https://github.com/kristyhoran/singscore
pip install TumorDecon
"""
import TumorDecon as td

"""
To perform digital cytometry on a mixture data, given cell signatures and the deconvolution method of your choice, run:

 td.tumor_deconvolve(mixture_data, method, patient_IDs='ALL', sig_matrix=None, up_genes=None, down_genes=None)

Explanation of Inputs:
    - mixture_data: pandas df of rna gene expression data.
        Rows are genes (indexed by 'Hugo_Symbol') and columns are patients
    - method: Method for tumor deconvolution. Must be either 'CIBERSORT', 'DeconRNASeq', 'ssGSEA', or 'SingScore'
    - patient_IDs: list of patient IDs to run DeconRNASeq for.
        Alternatively, can use the string 'ALL' to run for all patients
    - sig_matrix: pandas df of Signature gene expression values for given cell types.
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

## read in sample data (download colon cancer gene expressions directly from cbioportal)
## and fetch any missing Hugo Symbols from NCBI website:
rna = td.download_by_name('cbio', 'Colorectal Adenocarcinoma', fetch_missing_hugo=True)

## Can alternatively read in a data file already downloaded:
# rna = td.read_rna_file(data_loc+'coadred_data_RNA_Seq_v2_expression_median.txt')

################################################################################
############################# Linear Models ####################################
################################################################################
## Linear models (cibersort and DeconRNASeq) require a signature matrix, and output cell fractions:

## read in LM22 file as signature matrix
sig = td.read_sig_file(data_loc+'LM22.txt')
print("LM22 Signature Matrix:")
print(sig)

## We can also use the custom signature matrix created in 'sig_matrix_tutorial.py'.
    ## Note that this signature matrix uses Ensembl Gene IDs instead of Hugo Symbols.
    ## The read_sig_file() function can convert these to Hugo Symbols:
path_to_custom_sig = '~/TumorDecon/kmeans_signature_matrix_qval.txt' # Change to reflect your path
sig2 = td.read_sig_file(path_to_custom_sig, geneID='Ensembl_Gene_ID')
print("Custom Signature Matrix created in 'sig_matrix_tutorial.py':")
print(sig2)

## Run cibersort on ALL patients:
## optional argments include:
##    'scaling': how to scale the data (rna mixture & signature matrix) before applying DeconRNASeq (must be either 'None', 'zscore', or 'minmax')
##    'scaling_axis': whether to scale data by normalizing each cell types / patient individually (0) or each gene individually (1). default is 0
##    'nu','C','kernel': all arguments to be passed to sklearn's implementation of NuSVR. See sklearn documentation for details
##       additionally, can pass in 'nu':'best' to choose the nu from {0.25, 0.5, 0.75} that minimizes the root mean square error
##       defaults are nu='best', C=1.0, kernel=linear
##   'print_progress': whether to print patient ID as cibersort works through each patient (default: False)
ciber_freqs = td.tumor_deconvolve(rna, 'cibersort',  patient_IDs='ALL', sig_matrix=sig, args={'nu':'best', 'scaling':'minmax'})
print("CIBERSORT Results:")
print(ciber_freqs)

## Save results to a file:
ciber_freqs.to_csv("cibersort_results.csv", index_label='Patient_ID')

## Run DeconRNASeq:
## optional argments include:
##    'scaling': Same as in cibersort
##    'scaling_axis': Same as in cibersort
##    'check_sig': whether to check the condition number of the signature matrix before solving (default: False)
##    'print_progress': whether to print the results of each fit while function is running (default: False)
decon_freqs = td.tumor_deconvolve(rna, 'DeconRNASeq',  patient_IDs='ALL', sig_matrix=sig, args={'scaling':'minmax', 'print_progress':False})
print("DeconRNASeq Results:")
print(decon_freqs)

# ##############################################################################
# ######################### Rank-Based Methods #################################
# ##############################################################################
# Rank-based methods (ssGSEA and SingScore) require a list of up-regulated genes

## Read in up-regulated gene set derived from LM22:
# up_geneset = td.read_geneset(data_loc+"LM22_up_genes.csv")
up_geneset = td.read_geneset()
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
print("ssGSEA Results:")
print(ssgsea_scores)

## SingScore can be run with just an up-regulated gene set (unidirectional):
singscore_unidirectional = td.tumor_deconvolve(rna, 'singscore',  patient_IDs=patient_subset, up_genes=up_geneset)
print("SingScore Results:")
print(singscore_unidirectional)

## or with both an up-regulated and down-regulated gene set (bidirectional):
## (Can also derive up-regulated and down-regulated genes from a signature matrix, such as LM6 file):
LM6 = td.read_sig_file(data_loc+'LM6.txt')
up_geneset_LM6, down_geneset_LM6 = td.find_up_down_genes_from_sig(LM6, down_cutoff=0.4, up_cutoff=4.0)

singscore_bidirectional = td.tumor_deconvolve(rna, 'singscore',  patient_IDs=patient_subset, up_genes=up_geneset_LM6, down_genes=down_geneset_LM6)
print("Bidirectional SingScore Results:")
print(singscore_bidirectional)

################################################################################
######################### Visualize Results ####################################
################################################################################

# Select a subset of the patients to visualize:
results = ciber_freqs.loc[patient_subset]

## WARNING: The following 4 functions currently only work for fractions/ranks
## generated with LM22 cell signatures. If you are working with a different group
## of cell types, you will need to define your own dictionary for combing cell types
## An example of this is provided later on in this file.

## Can visualize results with boxplots. To save a plot, add the "save_as" argument:
td.cell_frequency_boxplot(results, save_as="boxplots.png")

## Can also visualize with barcharts:
td.cell_frequency_barchart(results, save_as="barcharts.png"))

# Pair Plots:
td.pair_plot(results, save_as="pairplots.png"))

# and Clustermaps.
td.hierarchical_clustering(results, save_as="clustermaps.png")

## Each of the plots above use matplotlib and seaborn, and users can customize
## plot parameters as in these packages. For example:
td.cell_frequency_boxplot(results, font_scale=0.75, rcParams={'figure.figsize':(5,3)})

## Note that these visualization functions simplify the results by summing the frequencies/scores of related cell types!
## You can use the "combine_celltypes" function to create such a simplified output dataframe.

## For example, LM22 contains 3 types of macrophages (M0,M1,M2). To combine these into one general
## column called "Macrophages", we define a dictionary, and pass it to the combine_celltypes() function:
dict = {'Macrophages':['Macrophages M0', 'Macrophages M1', 'Macrophages M2']}
print(td.combine_celltypes(results, cols_to_combine=dict))

## If no argument is passed in for "cols_to_combine", the function will attempt to use a sensible categorization
## based on the cell types present in the LM22 signature matrix:
print("\nOriginal Columns:")
print(results.columns)
results_simplified = td.combine_celltypes(results)
print("\nSimplified Columns:")
print(results_simplified.columns)
