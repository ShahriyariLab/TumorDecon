# data_utils.py - functions for reading in pancancer datasets and manipulating data

def get_td_Home():
    import os
    # Returns the path where this file is stored
    return os.path.realpath(__file__).strip("data_utils.py")

def read_rna_file(rna_file_path, identifier='hugo', fetch_missing_hugo=False):
    """
    Read in a cbioportal pancancer or TCGA Xena Hug txt file containing mixture gene expression data, and return a pandas dataframe.
    Columns are patients, Rows are genes
    Inputs:
        - rna_file_path: string. Relative or full path to the RNA seq file
            This file is tab seperated and includes two columns, 'Hugo_Symbol' and 'Entrez_Gene_Id', for each
            gene preceding the expression values for each patient
        - identifier: string. Determines which gene identifier to use to index the rna data.
            Must be set to either:
                - 'hugo' for Hugo Symbols
                - 'entrez' for Entrez Gene ID
        - fetch_missing_hugo: boolean. Whether to fetch missing Hugo Symbols (by Entrez Gene ID) from ncbi website
    Output:
        - pandas data frame (num_genes x num_patients) with either Hugo_Symbol or Entrez_Gene_Id as index
    """

    import pandas as pd
    from TumorDecon.hugoify import hugo

    rna = pd.read_csv(rna_file_path,  sep='\t')

    # Xena Hub:
    if "sample" in rna.columns:
        rna['Hugo_Symbol'] = rna['sample'].astype(str)
        rna = rna.set_index(['Hugo_Symbol']).drop(['sample'], axis=1)
    # CbioPortal TCGA Pan Cancer Atlas:
    elif "Hugo_Symbol" in rna.columns:
        rna[['Hugo_Symbol', 'Entrez_Gene_Id']] = rna[['Hugo_Symbol', 'Entrez_Gene_Id']].astype(str)
        if identifier == 'hugo':
            if fetch_missing_hugo:
                print("Fetching missing Hugo Symbols for genes by Entrez ID...")
                rna = hugo(rna)
            rna = rna.set_index(['Hugo_Symbol']).drop(['Entrez_Gene_Id'], axis=1)
            # Remove rows/genes that don't have a Hugo Symbol name:
            rna = rna.drop('nan',errors='ignore') # ignore error that arrises if no 'nan's
        elif identifier == 'entrez':
            rna = rna.set_index(['Entrez_Gene_Id']).drop(['Hugo_Symbol'], axis=1)
            # Remove rows/genes that don't have a Entrez ID:
            rna = rna.drop('nan')
        else:
            raise ValueError("gene identifier must be set to 'hugo' or 'entrez'")
    else:
        raise ReadError("Data Format not recognized")

    # Drop duplicate genes:
    # remove duplicates:
    # rna = rna.loc[~rna.index.duplicated(keep='first')] # argument keep = first, last, max, mean

    # Drop NA (if any):
    rna.dropna(axis=0, how='any', inplace=True)
    rna.dropna(axis=1, how='any', inplace=True)

    return rna


def download_by_name(source, type, fetch_missing_hugo=True):
    """
    Function to download TCGA data from either UCSC Xena Hub or cBioPortal given
    the desired source and cancer type, instead of a URL.
        Inputs:
            - source: either "xena" or "cbio"
            - type: specific strings. Each will be presented as options in a dropdown menu in the GUI
            - fetch_missing_hugo: boolean. Whether to fetch missing Hugo Symbols (by Entrez Gene ID) from ncbi website
        Outputs:
            - pandas df of mixture data
    """
    if source in ["xena", "Xena", "UCSC", "ucsc"]:
        type_dict = {'Acute Myeloid Leukemia':'LAML','Adrenocortical Cancer':'ACC', 'Bile Duct Cancer':'CHOL', 'Bladder Cancer':'BLCA', 'Breast Cancer':'BRCA',
                     'Cervical Cancer':'CESC', 'Colon and Rectal Cancer':'COADREAD', 'Colon Cancer':'COAD', 'Endometrioid Cancer':'UCEC', 'Esophageal Cancer':'ESCA',
                     'Glioblastoma':'GBM', 'Head and Neck Cancer':'HNSC', 'Kidney Chromophobe':'KICH', 'Kidney Clear Cell Carcinoma':'KIRC', 'Kidney Papillary Cell Carcinoma':'KIRP',
                     'Large B-cell Lymphoma':'DLBC', 'Liver Cancer':'LIHC', 'Lower Grade Glioma':'LGG', 'Lower Grade Glioma and Glioblastoma':'GBMLGG', 'Lung Adenocarcinoma':'LUAD',
                     'Lung Cancer':'LUNG', 'Lung Squamous Cell Carcinoma':'LUSC', 'Melanoma':'SKCM', 'Mesothelioma':'MESO', 'Ocular Melanomas':'UVM', 'Ovarian Cancer':'OV',
                     'Pancreatic Cancer':'PAAD', 'Pheochromocytoma and Paraganglioma':'PCPG', 'Prostate Cancer':'PRAD', 'Rectal Cancer':'READ', 'Sarcoma':'SARC', 'Stomach Cancer':'STAD',
                     'Testicular Cancer':'TGCT', 'Thymoma':'THYM', 'Thyroid Cancer':'THCA', 'Uterine Carcinosarcoma':'UCS'}
        if type in type_dict.keys():
            urlcode = type_dict[type]
        else:
            raise ValueError("type {!r} not available from UCSC Xena".format(type))
        # Download and return:
        data = download_from_xena(url="https://tcga-xena-hub.s3.us-east-1.amazonaws.com/latest/TCGA."+urlcode+".sampleMap/HiSeqV2.gz", fetch_missing_hugo=fetch_missing_hugo)
        return data
    elif source in ["cbio", "CbioPortal", "cbioportal", "cBioPortal"]:
        type_dict = {'Acute Myeloid Leukemia':'laml','Adrenocortical Carcinoma':'acc', 'Bladder Urothelial Carcinoma':'blca', 'Brain Lower Grade Glioma':'lgg', 'Breast Invasive Carcinoma':'brca',
                     'Cervical Squamous Cell Carcinoma':'cesc', 'Cholangiocarcinoma':'chol', 'Colorectal Adenocarcinoma':'coadread', 'Diffuse Large B-Cell Lymphoma':'dlbc', 'Esophageal Adenocarcinoma':'esca',
                     'Glioblastoma Multiforme':'gbm', 'Head and Neck Squamous Cell Carcinoma':'hnsc', 'Kidney Chromophobe':'kich', 'Kidney Renal Clear Cell Carcinoma':'kirc', 'Kidney Renal Papillary Cell Carcinoma':'kirp',
                     'Liver Hepatocellular Carcinoma':'lihc', 'Lung Adenocarcinoma':'luad', 'Lung Squamous Cell Carcinoma':'lusc', 'Mesothelioma':'meso', 'Ovarian Serous Cystadenocarcinoma':'ov', 'Pancreatic Adenocarcinoma':'paad',
                     'Pheochromocytoma and Paraganglioma':'pcpg', 'Prostate Adenocarcinoma':'prad', 'Sarcoma':'sarc', 'Skin Cutaneous Melanoma':'skcm', 'Stomach Adenocarcinoma':'stad', 'Testicular Germ Cell Tumors':'tgct',
                     'Thymoma':'thym', 'Thyroid Carcinoma':'thca', 'Uterine Carcinosarcoma':'ucs', 'Uterine Corpus Endometrial Carcinoma':'ucec', 'Uveal Melanoma':'uvm'}
        if type in type_dict.keys():
            urlcode = type_dict[type]
        else:
            raise ValueError("type {!r} not available from CbioPortal".format(type))
        # Download and return:
        data = download_from_cbio(url="https://cbioportal-datahub.s3.amazonaws.com/"+urlcode+"_tcga_pan_can_atlas_2018.tar.gz", fetch_missing_hugo=fetch_missing_hugo)
        return data
    else:
        raise ValueError("source ({!r}) must be 'xena' or 'cbioportal'".format(source))



def download_from_cbio(url="https://cbioportal-datahub.s3.amazonaws.com/uvm_tcga_pan_can_atlas_2018.tar.gz", save_location=get_td_Home()+"data/downloaded/", delete_after=False, fetch_missing_hugo=True):
    """
    Function to download data directly from cbioportal and read it in as a pandas df.
        Inputs:
            - url: url to data
            - save_loc: where to save the data. Default is in TumorDecon directory, under data/downloaded
            - delete_after: not implemented yet. Whether to delete data after reading
            - fetch_missing_hugo: boolean. Whether to fetch missing Hugo Symbols (by Entrez Gene ID) from ncbi website
        Outputs:
            - pandas df of mixture data
    """
    import wget
    import tarfile
    import os

    # Check if file already downloaded, if not, download it:
    file = save_location+url.replace('https://cbioportal-datahub.s3.amazonaws.com/','')
    if not os.path.exists(file):
        # Download file and save it locally:
        print("Downloading data from cbioportal...")
        wget.download(url, save_location)
    # Unzip if applicable
    folder = file.replace(".tar.gz","")
    if file.endswith("tar.gz") and not os.path.exists(folder):
        tar = tarfile.open(file, "r:gz")
        tar.extract(folder.split("/")[-1]+"/data_RNA_Seq_v2_expression_median.txt", path=save_location)
        # tar.extractall(folder.strip(".tar.gz"))
        tar.close()
    # Read in data:
    return read_rna_file(folder+"/data_RNA_Seq_v2_expression_median.txt", fetch_missing_hugo=fetch_missing_hugo)


def download_from_xena(url="https://tcga-xena-hub.s3.us-east-1.amazonaws.com/latest/TCGA.UCS.sampleMap/HiSeqV2.gz", save_location=get_td_Home()+"data/downloaded/", delete_after=False, fetch_missing_hugo=True):
    """
    Function to download data directly from TCGA Xena Hub and read it in as a pandas df.
        Inputs:
            - url: url to data
            - save_loc: where to save the data. Default is in TumorDecon directory, under data/downloaded
            - delete_after: not implemented yet. Whether to delete data after reading
            - fetch_missing_hugo: boolean. Whether to fetch missing Hugo Symbols (by Entrez Gene ID) from ncbi website
        Outputs:
            - pandas df of mixture data
    """
    import wget
    import gzip
    import os

    # Since all cancer types have the same filename, re-download every time:
    zipfile = save_location+url.split('/')[-1]
    if os.path.exists(zipfile):
        os.remove(zipfile)
    # Download file and save it locally:
    print("Downloading data from Xena Hub...")
    wget.download(url, save_location)
    # pandas read_csv can handle the zipped file - no need to unzip
    # Read in data:
    logdf = read_rna_file(zipfile, fetch_missing_hugo=fetch_missing_hugo)
    # Transform from log2(x+1):
    df = 2.**logdf - 1
    return df


# # OBSOLETE! Use read_sig_file() with no arguments to read in LM22.txt instead
# def read_lm22_file(file_path=get_td_Home()+"data/LM22.txt"):
#     """
#     Read in the LM22 (or LM6) signature matrix file (containing signature gene expression data for 22 different cell types)
#     and return a pandas dataframe
#     Inputs:
#         - file_path: string. Relative or full path to the LM22/LM6 signature matrix file
#             This file is tab seperated. Columns are cell types, Rows are gene Hugo Symbols.
#     Output:
#         - pandas data frame. Columns are cell types, Rows are genes, indexed by Hugo Symbol.
#     """
#
#     import pandas as pd
#
#     print("WARNING: 'read_lm22_file()'' depreciated. Use 'read_sig_file()' with no arguments to read in 'LM22.txt' signatures")
#
#     lm22 = pd.read_csv(file_path, sep='\t')
#     # lm22['Hugo_Symbol'] = lm22['Gene symbol']
#     # lm22 = lm22.drop(['Gene symbol'], axis = 1)
#     new_column_names = lm22.columns.tolist()
#     new_column_names[0] = 'Hugo_Symbol'
#     lm22.columns = new_column_names
#     lm22 = lm22.set_index(['Hugo_Symbol'])
#
#     return lm22


def read_sig_file(file_path=get_td_Home()+"data/LM22.txt", geneID="Hugo_Symbol"):
    """
    Reads in a signature matrix file (containing signature gene expression data for
    a number of different cell types), converts the gene identifiers to Hugo Symbols,
    and returns a pandas dataframe.
        ** IF NO FILE PATH GIVEN, LM22 IS ASSUMED **
    Inputs:
        - file_path: string. Relative or full path to signature matrix file.
            This file is tab seperated. Columns are cell types, Rows are genes
            If no file_path given, LM22 is assumed.
        - geneID: string in ["Hugo_Symbol", "Ensembl_Gene_ID", "Entrez_Gene_ID"].
            Describes how genes are labeled in the signature matrix file
    Output:
        - pandas data frame. Columns are cell types, Rows are genes, indexed by Hugo Symbol.
    """
    import pandas as pd

    # Read file:
    sig = pd.read_csv(file_path, sep='\t')
    # Get name of column with gene identifiers:
    gene_col = list(sig)[0]

    to_be_dropped = []
    if geneID == "Hugo_Symbol":
        sig['Hugo_Symbol'] = sig[gene_col]

    elif geneID == 'Ensembl_Gene_ID':
        sig['Hugo_Symbol'], to_be_dropped = convert_emsembl(sig[gene_col])

    elif geneID == 'Entrez_Gene_ID':
        raise ValueError("Conversion from Entrez Gene ID to Hugo Symbol not yet supported in TumorDecon")

    else:
        raise ValueError("geneID ({!r}) must be set to 'Hugo_Symbol', 'Ensembl_Gene_ID', or 'Entrez_Gene_ID'".format(geneID))

    # Drop old gene ID column, rename columns, and set Index to Hugo Symbol:
    sig = sig.drop([gene_col], axis = 1)
    new_column_names = sig.columns.tolist()
    sig.columns = new_column_names
    sig = sig.set_index(['Hugo_Symbol'])

    # Drop columns without a Hugo Symbol:
    sig = sig.drop(to_be_dropped)
    return sig

def convert_emsembl(ensembl):
    """
    Input:
        - ordered array of Ensembl Gene IDs to be converted to Hugo Symbols
    Outputs:
        - hugo: array. Ordered array of Hugo Symbols correponsing to the Ensembl Gene IDs
        - to_be_dropped: array. List of genes without a correponsing hugo symbol (to be dropped)
    """
    import pandas as pd
    import numpy as np

    hugo = np.array(ensembl)
    # Will drop any genes that don't have a corresponding Hugo Symbol
    to_be_dropped = []

    # Create dictionary for conversion:
    ensembl_hugo = pd.read_csv(get_td_Home()+"data/gene_ID_conversion/hugo_enseml_synonym.txt", header=0, sep="\t")
    ensembl_hugo2 = pd.read_csv(get_td_Home()+"data/gene_ID_conversion/GSE87692_Primary_HsESF_TPMs_073015.txt", header=0, sep="\t")
    ensembl_hugo2 = ensembl_hugo2[['ensemble_geneid', 'Associate Gene Name']]
    ensembl_hugo2.set_index('Associate Gene Name', inplace=True)
    ensembl_hugo.dropna(subset=['Ensembl ID(supplied by Ensembl)', 'Ensembl gene ID'], how='all', inplace=True)
    ensembl_hugo = ensembl_hugo.set_index('Approved symbol')

    for i, gene in enumerate(ensembl):
        # first search in the ensembl_hugo2 then if not found in ensembl_hugo
        is_found = ensembl_hugo2.loc[ensembl_hugo2['ensemble_geneid'] == gene]
        if len(is_found) > 0 and is_found.index[0] != "#N/A":
            hugo[i] = is_found.index[0]

        else:
            t = ensembl_hugo.loc[ensembl_hugo['Ensembl ID(supplied by Ensembl)'] == gene]
            # first search in 'Ensembl ID(supplied by Ensembl)' column, if not found, search the other col
            if len(t) > 0:
                if type(t.index[0]) == str:
                    hugo[i] = t.index[0]

            else:
                t2 = ensembl_hugo.loc[ensembl_hugo['Ensembl gene ID'] == gene]
                if len(t2) > 0:
                    if type(t2.index[0]) == str:
                        hugo[i] = t2.index[0]
                    else:
                        to_be_dropped.append(gene)
                else:
                    to_be_dropped.append(gene)
    if len(to_be_dropped) > 0:
        print("Missing Hugo Symbols for the following genes (will be dropped from signature matrix):")
        print(to_be_dropped)
    return hugo, to_be_dropped


def df_normalization(df, scaling, axis=0):
    """
    Inputs:
        - df: pandas data frame
        - scaling: string, must be either 'zscore' or 'minmax'.
        - axis: 0 or 1. Whether to apply normalization across individual columns (0) or individual rows (1)
    Output:
        - df_norm: normalized (minmax) data
    """
    import pandas as pd
    import numpy as np
    from sklearn import preprocessing
    if scaling == "zscore":
        df_norm = preprocessing.scale(df, axis=axis)
        df_norm = pd.DataFrame(df_norm)
    elif scaling == "r-zscore":
        # uses delta degrees of freedom = 1 (as R std does), instead of 0 (as sklearn does)
        mean = np.mean(df, axis=axis)
        std = np.std(df, axis=axis, ddof=1)
        df_norm = (df - mean) / std
        df_norm = pd.DataFrame(df_norm)
    elif scaling == "minmax":
        min_max_scaler = preprocessing.MinMaxScaler()
        if axis == 1:
            df_norm = df.transpose()
        else:
            df_norm = df.copy()
        df_norm = min_max_scaler.fit_transform(df_norm)
        df_norm = pd.DataFrame(df_norm)
        if axis == 1:
            df_norm = df_norm.transpose()
    # add column names back to the dataframe
    df_norm.columns = df.columns
    # add row index back to the dataframe
    df_norm.index = df.index
    return df_norm


def keep_common_genes(rna_df, sig_df):
    """
    Given two dataframes, eliminate all genes that are not present in both dataframes.
    Both dataframes must be indexed with the same gene identifier
    Inputs:
        - rna_df: pandas df. Rows are genes (indexed by 'Hugo_Symbol') and columns are patients
        - sig_df: pandas df. Rows are genes (indexed by 'Hugo_Symbol') and columns are cell types
    Output:
        - rna_df_red: pandas data frame. rna_df reduced to only contain genes also present in sig_df
        - sig_df_red: pandas data frame. sig_df reduced to only contain genes also present in rna_df
    """
    import pandas as pd

    # Eliminate genes that are not present in both rna df and signature matrix
    shared_index = pd.merge(rna_df, sig_df, how='inner', left_index=True, right_index=True)
    rna_df_red = rna_df.loc[shared_index.index]
    sig_df_red = sig_df.loc[shared_index.index]

    # Make sure genes listed in same order for both dataframes lines up:
    rna_df_red = rna_df_red.sort_index()
    sig_df_red = sig_df_red.sort_index()

    return rna_df_red, sig_df_red


def find_up_down_genes_from_sig(sig_df, down_cutoff=0.4, up_cutoff=4.0, show_plots=False):
    """
    Function divides each gene expression value in Signature Matrix (sig_df) by the median value of the
        given gene across all cell types.
    Then all genes with ratios below 'down_cutoff' are considered "down regulated genes" while all genes with
    ratios above 'up_cutoff' are considered "up regulated genes".
    Inputs:
        - sig_df: Pandas df. Signature matrix (rows are genes, columns are cell types) indexed by Hugo Symbol
        - down_cutoff: value to use for the cutoff point as to what should be considered a "down regulated gene"
        - up_cutoff: value to use for the cutoff point as to what should be considered a "up regulated gene"
        - show_plots: boolean. Whether to show a plot of the sorted signature matrix ratios for each cell type.
            Can be used to help choose the cutoff points
    Outputs:
        - up: dictionary. Keys are cell types, values are a list of up-regulated genes for that cell type
        - down: dictionary. Keys are cell types, values are a list of down-regulated genes for that cell type
    """
    import pandas as pd
    import matplotlib.pyplot as plt

    ref = sig_df.copy()
    ref['reference'] = ref.median(axis=1)
    sig_ratios = sig_df.divide(ref['reference'], axis=0)

    if show_plots:
        t = range(sig_ratios.shape[0])
        for cell in sig_ratios.columns:
            plt.figure()
            vals = sig_ratios[cell]
            vals = vals.sort_values()
            axes = plt.gca()
            axes.set_ylim([0,10])
            axes.plot(t, vals)
            plt.title(cell)
            plt.show()

    # Generate list of up-regulated and down-regulated gene sets
    up = {}
    down = {}
    for cell in sig_ratios.columns:
        down[cell] = sig_ratios[sig_ratios[cell] < down_cutoff].index.tolist()
        up[cell] = sig_ratios[sig_ratios[cell] > up_cutoff].index.tolist()

    return up, down


def get_top_ranked_genes_from_sig(sig_df, sig_size=50):
    """
    Function returns the top 'sig_size' most expressed genes in the signature matrix, for each cell type
    Inputs:
        - sig_df: Pandas df. Signature matrix (rows are genes, columns are cell types) indexed by Hugo Symbol
        - sig_size: int. Number of genes to select
    Outputs:
        - dictionary. Keys are cell types, values are a list of the "sig_size" top most expressed genes in the
        signature matrix for that cell type
    """
    import pandas as pd

    top_genes = pd.DataFrame()
    for cell in sig_df.columns:
        top_genes[cell] = sig_df[cell].sort_values(ascending=False).index[:sig_size]

    return top_genes.to_dict('list')



def variance_threshold_selector(data, threshold=0.5):
    """
    Function to fit sklearn's variance threshold feature selection method on a df where the features (genes) are
    the rows in a pandas df (instead of columns)
    Inputs:
        - data: pandas df with genes as rows and either patients/celltype as column
        - threshold: The cut-off value to use for variance threshold selection
    Outputs:
        - reduced data frame (eliminates genes with variance < threshold) across patients/celltypes
    """
    from sklearn.feature_selection import VarianceThreshold

     # Transpose to work with sklearn's implementation, which assumes features are columns:
    data = data.T
    # Fit selector:
    selector = VarianceThreshold(threshold)
    selector.fit(data)
    # Apply transformation to df:
    data_red = data[data.columns[selector.get_support(indices=True)]]
    # Transpose back so genes are rows again:
    return data_red.T


def read_ssGSEA_up_genes(file_path=get_td_Home()+'data/Gene_sets.csv'):
    """
    Function to read in csv file containing the gene sets from the ssGSEA paper
    Inputs:
        none (hard-coded path to correct csv file)
    Outputs:
        - dictionary. Keys are the 26 cell types used in ssGSEA paper, values are a list of up regulated genes
            for that cell type
    """
    import pandas as pd

    gene_sets = pd.read_csv(file_path)
    up_genes = {}
    for cell in gene_sets['Cell type'].unique():
        up_genes[cell] = gene_sets[gene_sets['Cell type'] == cell]['Symbol'].tolist()

    return up_genes

def read_geneset(file_path=get_td_Home()+"data/LM22_up_genes.csv"):
    """
        Function to read in a custom csv file containing the up or down regulated
            genes for each cell type
        Inputs:
            - file_path: string. Relative or full path to signature matrix file.
            File should have columns named for each cell type, and rows should
            contain Hugo Symbols for genes to be considered as up (or down)
            regulated for that cell type. If not all cell types have the same number of up(down)
            regulated genes, excess rows in each column should be coded as "NaN"
                If no file_path given, the up-regulated gene set derived in
                Le et al. (2020) is assumed.
        Outputs:
            - dictionary. Keys are the column names of the given csv file (cell types),
                values are a list of up (down) regulated genes for that cell type
    """
    import pandas as pd
    import numpy as np
    gene_sets = pd.read_csv(file_path)
    chosen_genes = {} # up or down genes
    for cell in gene_sets.columns.unique():
        chosen_genes[cell] = list(gene_sets[cell].dropna())
    return chosen_genes


def combine_celltypes(df, cols_to_combine=None):
    """
    Function to sum related cell types into a single column
    Inputs:
        - df: pandas dataframe. Output of td.tumor_deconvolve()
        - cols_to_combine: dictionary. Keys are the desired names of any new cell type columns,
            values are an arary of current column names to combine under the key name
                - Default = dictionary for combining common cell types from LM22
    Outputs:
        - Pandas dataframe with columns combined as specified by the dictionary.
            All unmentioned column names are left as they are
    """
    import pandas as pd

    if cols_to_combine is not None:
        from collections.abc import Mapping
        if isinstance(cols_to_combine, Mapping):
            pass
        else:
            raise TypeError("cols_to_combine must be a dictionary")
    else:
        # Use LM22 as default
        print("WARNING: No dictionary defined for combining columns... Attempting to use default dict for LM22 signatures")
        cols_to_combine = {'B cells':['B cells naive', 'B cells memory'],
                           'CD4 T cells':['T cells CD4 naive', 'T cells CD4 memory resting', 'T cells CD4 memory activated','T cells follicular helper', 'T cells regulatory (Tregs)'],
                           'CD8 T cells':['T cells CD8'],
                           'NK cells':['NK cells resting', 'NK cells activated'],
                           'Macrophages':['Macrophages M0', 'Macrophages M1', 'Macrophages M2'],
                           'Mast cells':['Mast cells resting','Mast cells activated'],
                           'Dendritic cells':['Dendritic cells resting', 'Dendritic cells activated']
                          }
    df2 = df.copy()
    for cell_type in cols_to_combine.keys():
        try:
            df2[cell_type] = df2[cols_to_combine[cell_type]].sum(axis=1)
            df2.drop(cols_to_combine[cell_type],axis=1,inplace=True)
        except KeyError as e:
            print("WARNING: Failed to combine some columns: ")
            print("KeyError: "+str(e))
            pass

    return df2



def corr_table(methods, results, cell_types, true_freqs):
    import pandas as pd
    import numpy as np
    from scipy.stats.stats import pearsonr
    from scipy.stats import spearmanr

    p_corr_per_cell = pd.DataFrame(index=cell_types, columns=methods)
    p_corr_per_sample = pd.DataFrame(index=true_freqs.index, columns=methods)
    s_corr_per_cell = pd.DataFrame(index=cell_types, columns=methods)
    s_corr_per_sample = pd.DataFrame(index=true_freqs.index, columns=methods)
    for method in methods:
        for cell in cell_types:
            p_corr_per_cell.loc[cell, method] = pearsonr(results[method][cell], true_freqs[cell])[0]
            s_corr_per_cell.loc[cell, method] = spearmanr(results[method][cell], true_freqs[cell])[0]
        for i in range(len(true_freqs)):
            p_corr_per_sample.loc[true_freqs.index[i], method] = pearsonr(results[method][cell], true_freqs[cell])[0]
            s_corr_per_sample.loc[true_freqs.index[i], method] = spearmanr(results[method][cell], true_freqs[cell])[0]
    return np.abs(p_corr_per_cell), np.abs(p_corr_per_sample), np.abs(s_corr_per_cell), np.abs(s_corr_per_sample)


def corr_mean_std(corr_per_cell, corr_per_sample):
    import pandas as pd
    import numpy as np

    corr = pd.DataFrame(data=np.zeros((corr_per_cell.shape[1], 4)), index=corr_per_cell.columns, columns=['Mean_corr_per_sample', 'Std_corr_per_sample', 'Mean_corr_per_cell', 'Std_corr_per_cell'])
    corr['Mean_corr_per_sample'] = np.mean(corr_per_sample, axis=0)
    corr['Std_corr_per_sample'] = np.std(corr_per_sample, axis=0)
    corr['Mean_corr_per_cell'] = np.mean(corr_per_cell, axis=0)
    corr['Std_corr_per_cell'] = np.std(corr_per_cell, axis=0)
    return corr



def flatten_corr_per_cell(corr_per_cell):
    import pandas as pd
    import numpy as np

    corr_per_cell2 = pd.DataFrame(data=np.zeros((corr_per_cell.shape[0]*corr_per_cell.shape[1], 3)), columns=['Method', 'Cell_type', 'Correlation'])
    methods = []
    for method in corr_per_cell.columns:
        methods.extend([method]*corr_per_cell.shape[0])
    corr_per_cell2['Method'] = methods
    corr_per_cell2['Cell_type'] = list(corr_per_cell.index)*corr_per_cell.shape[1]
    corr = []
    for method in corr_per_cell.columns:
        corr.extend(list(corr_per_cell[method]))
    corr_per_cell2['Correlation'] = corr

    return corr_per_cell2


def predicted_truth_bycell(method, method_freqs, exp_freqs, cell_types):
    import pandas as pd

    df = pd.concat([method_freqs[cell_types], exp_freqs[cell_types]])
    df['Method'] = [method]*exp_freqs.shape[0] + ['Ground truth']*exp_freqs.shape[0]
    df = pd.DataFrame(data=np.zeros((exp_freqs.shape[0]*len(cell_types), 3)),
                              columns=[method, 'Ground truth', 'Cell type'])
    method_fractions, true_fractions, cell_list = [], [], []
    for cell in cell_types:
        method_fractions.extend(list(method_freqs[cell]))
        true_fractions.extend(list(exp_freqs[cell]))
        cell_list.extend([cell]*exp_freqs.shape[0])
    df[method] = method_fractions
    df['Ground truth'] = true_fractions
    df['Cell type'] = cell_list

    return df
