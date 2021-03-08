# sig_matrix_tutorial.py - authored by Rachel Aronow & Shaya Akbarinejad

"""
This tutorial walks through an example of generating a signature matrix from raw
data files. The data used in this tutorial are single cell profiles from the Gene
Expression Omnibus (GEO) repository. They are available to download on the NCBI
website, as well as on from our GitHub page, at:
https://github.com/ShahriyariLab/TumorDecon/tree/master/TumorDecon/data/sig_matrix_tutorial
"""

import TumorDecon as td

################################################################################
######### Example 1: Build a signature matrix from raw data files ##############
################################################################################

# Sample data location - change to reflect where you have downloaded the data:
data_loc = td.get_td_Home()+"data/sig_matrix_tutorial/"

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

# Run the batch correction function, and save the corrected data under the name "batch_corrected_sample_data.txt"
td.batch_correct_datasets(data_loc, cell_file_dict, outfile="batch_corrected_sample_data.txt")


# Call the function to create the signature matrix from our batch corrected dataset:
# This will output 2 text files:
#   - "batch_corrected_sample_data_clustered.txt"
#   - "kmeans_signature_matrix_qval.txt"
cell_types_to_include = ["CD8", "CD4", "B_cell", "NK", "mono", "Endothelial", "Fibroblast", "Neutrophils"]
td.create_signature_matrix("batch_corrected_sample_data.txt", cell_types_to_include)


# Note that the "create_signature_matrix()" function expects a non-clustered data file as input.
   # However, if your datafile is large, running these two processes (clustering & generation of sig matrix)
   # in sequence can take a very long time, or you may run into memory issues.
   # To avoid this, you can also input an already-clustered data file (saved as an intermediate step in a previous run),
   # and include the argument clustered=True, to skip the clustering step in following runs:
td.create_signature_matrix("batch_corrected_sample_data_clustered.txt", cell_types_to_include, clustered=True, outfile="kmeans_signature_matrix_qval_2.txt")

# The two signature matrices created above should be identical

# We can then use this signature matrix 'kmeans_signature_matrix_qval.txt' to run TumorDecon! (tutorial.py)

################################################################################
######### Example 2: Build signature matrix from LM22 source profiles ##########
################################################################################

# Coming Soon!
