# data_utils.py - functions for reading in pancancer datasets and manipulating data

def read_pancancer_rna_file(rna_file_path, identifier='hugo'):
    """
    Read in a pancancer csv file containing mixture gene expression data, and return a pandas dataframe.
    Columns are patients, Rows are genes
    Inputs:
        - rna_file_path: string. Relative or full path to the RNA seq file
            This file is tab seperated and includes two columns, 'Hugo_Symbol' and 'Entrez_Gene_Id', for each
            gene preceding the expression values for each patient
        - identifier: string. Determines which gene identifier to use to index the rna data.
            Must be set to either:
                - 'hugo' for Hugo Symbols
                - 'entrez' for Entrez Gene ID
    Output:
        - pandas data frame (num_genes x num_patients) with either Hugo_Symbol or Entrez_Gene_Id as index
    """

    import pandas as pd

    rna = pd.read_csv(rna_file_path,  sep='\t')
    rna[['Hugo_Symbol', 'Entrez_Gene_Id']] = rna[['Hugo_Symbol', 'Entrez_Gene_Id']].astype(str)

    if identifier == 'hugo':
        rna = rna.set_index(['Hugo_Symbol']).drop(['Entrez_Gene_Id'], axis=1)
    elif identifier == 'entrez':
        rna = rna.set_index(['Entrez_Gene_Id']).drop(['Hugo_Symbol'], axis=1)
    else:
        raise ValueError("gene identifier must be set to 'hugo' or 'entrez'")

    return rna



def read_lm22_file(lm22_file_path):
    """
    Read in the LM22 signature matrix file (containing signature gene expression data for 22 different cell types)
    and return a pandas dataframe
    Inputs:
        - lm22_file_path: string. Relative or full path to the LM22 signature matrix file
            This file is tab seperated. Columns are cell types, Rows are gene Hugo Symbols.
    Output:
        - pandas data frame. Columns are cell types, Rows are genes, indexed by Hugo Symbol.
    """

    import pandas as pd

    lm22 = pd.read_csv(lm22_file_path, sep='\t')
    # lm22['Hugo_Symbol'] = lm22['Gene symbol']
    # lm22 = lm22.drop(['Gene symbol'], axis = 1)
    new_column_names = lm22.columns.tolist()
    new_column_names[0] = 'Hugo_Symbol'
    lm22.columns = new_column_names
    lm22 = lm22.set_index(['Hugo_Symbol'])

    return lm22

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
    from sklearn import preprocessing
    if scaling == "zscore":
        df_norm = preprocessing.scale(df, axis=axis)
        df_norm = pd.DataFrame(df_norm)
    elif scaling == "minmax":
        min_max_scaler = preprocessing.MinMaxScaler()
        if axis == 1:
            df_norm = df.transpose()
        df_norm = min_max_scaler.fit_transform(df_norm)
        df_norm = pd.DataFrame(df_norm)
        if axis == 1:
            df_norm = df_norm.transpose()
    # add column names back to the dataframe
    df_norm.columns = df.columns
    # add row index back to the dataframe
    df_norm.index = df.index
    return df_norm

# def read_custom_geneset():
# columns are cell types, rows are lists of gene identifiers (hugo symbols)

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


def read_ssGSEA_up_genes(filepath='data/Gene_sets.csv'):
    """
    Function to read in csv file containing the gene sets from the ssGSEA paper
    Inputs:
        none (hard-coded path to correct csv file)
    Outputs:
        - dictionary. Keys are the 26 cell types used in ssGSEA paper, values are a list of up regulated genes
            for that cell type
    """
    import pandas as pd

    gene_sets = pd.read_csv(filepath)
    up_genes = {}
    for cell in gene_sets['Cell type'].unique():
        up_genes[cell] = gene_sets[gene_sets['Cell type'] == cell]['Symbol'].tolist()

    return up_genes
