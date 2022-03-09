# Reproducible Example, demonstrating all the capabilities of the TumorDecon package

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

## read in sample data (download colon cancer gene expressions directly from cbioportal)
## and fetch any missing Hugo Symbols from NCBI website:
rna = td.download_by_name('cbio', 'Uveal Melanoma', fetch_missing_hugo=False)

## Can alternatively read in a data file already downloaded:
## The location "td.get_td_Home()" refers to the install location of the TumorDecon package
data_loc = td.get_td_Home()+"data/"
# rna = td.read_rna_file(data_loc+'coadred_data_RNA_Seq_v2_expression_median.txt')


"""
The following section walks through an example of generating a signature matrix from raw
data files. The data used in this tutorial are single cell profiles from the Gene
Expression Omnibus (GEO) repository. They are available to download on the NCBI website,
as well as from our GitHub page, at:

https://github.com/ShahriyariLab/TumorDecon/tree/master/TumorDecon/data/sig_matrix_tutorial
"""

################################################################################
########  BEGIN: Building a signature matrix from raw data files ###############
################################################################################
# Define location where data is downloaded:
single_cell_data_loc = "sig_matrix_tutorial/" # CHANGE THIS to where you have downloaded the data

# Define dictionary of cell types and all raw files that contain data for that cell type:
cell_file_dict = {'CD4': ['GSE107011_Processed_data_TPM2.txt', 'ensembl_version_GSE97861.txt', 'ensembl_version_GSE97862.txt',
                          'ensembl_version_GSE113891.txt', 'ensembl_version_GSE114407.txt', 'ensembl_version_GSE115978.txt'],
                  'CD8': ['ensembl_version_GSE98638.txt', 'ensembl_version_GSE114407.txt', 'GSE107011_Processed_data_TPM2.txt'],
                  'B_': ['GSE107011_Processed_data_TPM2.txt', 'ensembl_version_GSE114407.txt', 'ensembl_version_GSE115978.txt'],
                  'mono': ['ensembl_version_GSE114407.txt', 'GSE107011_Processed_data_TPM2.txt'],
                  'NK': ['GSE107011_Processed_data_TPM2.txt', 'ensembl_version_GSE115978.txt'],
                  'Endo': ['ensembl_version_GSE102767.txt', 'ensembl_version_GSE113839.txt', 'ensembl_version_GSE115978.txt'],
                  'Fibro': ['ensembl_version_GSE113839.txt', 'GSE109448_rnaseq_gene_tpm.tsv', 'GSE109449_singlecell_rnaseq_gene_tpm.txt'],
                  'Neutro': ['GSE107011_Processed_data_TPM2.txt']}

## Run the batch correction function, and save the corrected data under the name "batch_corrected_sample_data.txt"
td.batch_correct_datasets(single_cell_data_loc, cell_file_dict, outfile="batch_corrected_sample_data.txt")

## Call the function to create the signature matrix from our batch corrected dataset:
## This will output 2 text files:
#   - "batch_corrected_sample_data_clustered.txt"
#   - "kmeans_signature_matrix_qval.txt"
cell_types_to_include = ["CD8", "CD4", "B_cell", "NK", "mono", "Endothelial", "Fibroblast", "Neutrophils"]
td.create_signature_matrix("batch_corrected_sample_data.txt", cell_types_to_include)


## Note that the "create_signature_matrix()" function expects a non-clustered data file as input.
   # However, if your datafile is large, running these two processes (clustering & generation of sig matrix)
   # in sequence can take a very long time, or you may run into memory issues.
   # To avoid this, you can also input an already-clustered data file (saved as an intermediate step in a previous run),
   # and include the argument clustered=True, to skip the clustering step in following runs:
td.create_signature_matrix("batch_corrected_sample_data_clustered.txt", cell_types_to_include, clustered=True, outfile="kmeans_signature_matrix_qval_2.txt")

## The two signature matrices created above will be identical

## We can then use this signature matrix 'kmeans_signature_matrix_qval.txt' to run TumorDecon! (next step!)

################################################################################
##########  END: Building a signature matrix from raw data files ###############
################################################################################

################################################################################
####################### BEGIN: Linear Models ###################################
################################################################################
## Linear models (cibersort and DeconRNASeq) require a signature matrix, and output cell fractions:

# Read in default signature matrix file (LM22):
sig = td.read_sig_file()
# (Can also pass in a filename of an alternative signature)
print("LM22 Signature Matrix:")
print(sig)

## We can also use the custom signature matrix we created above:
    ## Note that this signature matrix uses Ensembl Gene IDs instead of Hugo Symbols.
    ## The read_sig_file() function can convert these to Hugo Symbols, by specifying the geneID:
path_to_custom_sig = 'kmeans_signature_matrix_qval.txt' # CHANGE this to match where you saved the file above
sig2 = td.read_sig_file(path_to_custom_sig, geneID='Ensembl_Gene_ID')
print("Custom Signature Matrix created:")
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
decon_freqs.to_csv("DeconRNASeq_results.csv", index_label='Patient_ID')

################################################################################
######################### END: Linear Models ###################################
################################################################################

################################################################################
####################### BEGIN: Rank-Based Methods ##############################
################################################################################
# Rank-based methods (ssGSEA and SingScore) require a list of up-regulated genes

# Read in default up-regulated gene set (derived from single-cell reference profiles of LM22, in Le et al. (2020)):
up_geneset = td.read_geneset()
# (Can also pass in a filename of an alternative gene set)

####
# When running tumor_deconvolve(), the default is to apply the method on all patients
# in the rna data frame (patient_IDs='ALL'). To instead run the method on only a
# select list of patients, pass in a list of patient identifiers (column names) instead
####
patient_subset = ['TCGA-RZ-AB0B-01', 'TCGA-V3-A9ZX-01', 'TCGA-V3-A9ZY-01', 'TCGA-V4-A9E5-01',
                  'TCGA-V4-A9E7-01', 'TCGA-V4-A9E8-01', 'TCGA-V4-A9E9-01', 'TCGA-V4-A9EA-01',
                  'TCGA-V4-A9EC-01', 'TCGA-V4-A9ED-01',	'TCGA-V4-A9EE-01', 'TCGA-V4-A9EF-01']

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
######################### END: Rank-Based Methods ##############################
################################################################################


################################################################################
###################### BEGIN: Visualize Results ################################
################################################################################

## We will visualize the cibersort results:
results = ciber_freqs

## Can also select just a subset of the patients to visualize:
# results = ciber_freqs.loc[patient_subset]

## WARNING: The following 4 functions assume the 22 cell types present in LM22.
## If you are working with a signature matrix / different group of cell types, you
## will need to define your own dictionary for combing cell types - an example of
## this is provided later on in this tutorial.

## Each of the following plots use matplotlib and seaborn, and users can pass in
## plot parameters for these packages such as font_scale, rcParams, etc
## To save a plot, add the "save_as" argument:

## Can visualize results with boxplots:
td.cell_frequency_boxplot(results, title="Distribution of Cell Types in all Samples", font_scale=0.75, title_fontsize=15, save_as="example_boxplot.png")
# CHANGE PATH
## Barcharts:
td.cell_frequency_barchart(results, title="Distribution of Cell Types within Individual Samples", save_as="example_barchart.png")
# CHANGE PATH
# Cluster Maps:
td.hierarchical_clustering(results, title="Cluster Map", figsize=(15,15), save_as="example_clutermap.png")
# CHANGE PATH
# and Pair Plots:
td.pair_plot(results, save_as="example_pairplot.png")
# CHANGE PATH

## As mentioned above, these visualization functions simplify the results by summing the frequencies/scores of related cell types!
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

################################################################################
######################## END: Visualize Results ################################
################################################################################
