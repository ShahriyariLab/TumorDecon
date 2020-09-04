# signature_matrix.py - functions originally authored by Shaya Akbarinejad

def create_signature_matrix(infile, cell_types, clustered=False, max_clusters=10, intermfile=None, outfile="kmeans_signature_matrix_qval.txt"):
    """
    Function:
        -takes in batch corrected datasets as input, does clustering on each cell-type
        (with the silhoette method), does differential expression analysis to get the genes that
        are significantly expressed among one cluster (using adjusted p value) and outputs signature matrix
    Inputs:
        - infile: string. Path to text file of batch corrected datasets (output of "batch_correct_datasets()" function)
        - cell_types: List of cell types to include in signature matrix
        - clustered: boolean. Whether the infile is pre-clustered or not.
        - max_clusters: int. Maximum value of K to ty for k-means clustering.
            Ignored if clustered=True
        - intermfile: string. Where to save intermediate step of clustered batch corrected datasets.
            Ignored if clustered=True
            Default is infile+"_clustered.txt"
        - outfile: string. Where to save output files to.
            Default is "kmeans_signature_matrix_qval.txt"
    Outputs:
        - Signature Matrix (as named by "outfile")
        - (if clustered=False & intermfile != None) text file of Clustered (via silhoette method) Batch Corrected Datasets
    """

    import pandas as pd

    print("Reading batch-corrected dataset file %s..." % infile)

    # Cluster each cell if not already clustered
    if not clustered:
        all_cells = pd.read_csv(infile, sep="\t")
        all_cells.set_index(list(all_cells)[0], inplace=True)
        print("Clustering batch-corrected datasets...")
        clustered_cells = cluster_each_cell(cell_types, all_cells, max_clusters)
        # Save clustered results to text file (intermediate step):
        if intermfile is None:
            intermfile = infile.replace(".txt","_clustered.txt")
        print("Saving clustered data sets to %s..." % intermfile)
        clustered_cells.to_csv(intermfile, sep="\t")

    # If already clustered, simply read in clustered data file:
    if clustered:
        clustered_cells = pd.read_csv(infile, sep="\t")
        clustered_cells.set_index(list(clustered_cells)[0], inplace=True)

    # Compute signature matrix and save to 'outfile':
    clustered_cells[clustered_cells < 0] = 0
    print("Generating signature matrix from clustered batch-corrected datasets...")
    mean_of_cluster = get_mean_of_each_cluster(cell_types, clustered_cells)
    final_sign_mat = get_differentially_expr_genes(cell_types, clustered_cells, mean_of_cluster)
    print("Saving signature matrix to %s..." % outfile)
    final_sign_mat.to_csv(outfile, sep="\t")



def find_optimal_k(cell_expr, max_clusters):
    """
    Args:
        cell_expr: (n_samples, n_genes) pd DataFrame - Transposed DataFrame of cell expression
        max_clusters: integer - Specifies maximum number of of clusters to be considered
    Returns:
        k_opt: integer - optimal number of clusters based on the Silhouette method
    """
    import numpy as np
    from sklearn.cluster import KMeans
    from sklearn.metrics import silhouette_score

    sample_num = cell_expr.shape[0]
    k_vals = np.zeros(len(range(2, min(max_clusters, sample_num))))
    range_n_clusters = list(range(2, min(max_clusters, sample_num)))
    for i, n_clusters in enumerate(range_n_clusters):
        clusterer = KMeans(n_clusters=n_clusters, random_state=10, n_init=10, n_jobs=-1, tol=1e-3)
        cluster_labels = clusterer.fit_predict(cell_expr)
        silhouette_avg = silhouette_score(cell_expr, cluster_labels)
        k_vals[i] = silhouette_avg

    k_opt = np.argmax(k_vals) + 2
    print("Optimal K = %d" % k_opt)
    return k_opt


def cluster_each_cell(cells, cell_expr, max_clusters):
    """
    Args:
        cells: list of strings - Specifies what type of cells we have
        cell_expr: (n_genes, n_samples) pd DataFrame - DataFrame of cell expression
    Returns:
        out_df: (n_genes, n_samples) pd DataFrame - clustered cell expression. each column shows the cell type and the
                cluster number that sample belongs to, e.g. CD4_subtype_1.
    """
    import pandas as pd
    from sklearn.cluster import KMeans

    out_df = pd.DataFrame(index=cell_expr.index)
    for cell in cells:
        one_cell_expr = cell_expr.loc[:, cell_expr.columns.str.contains(cell)]

        assert(one_cell_expr.shape[1] != 0), cell + " is not present in the file"

        print("%s: Finding optimal K for K-means (max # of clusters = %d)" % (cell, max_clusters))
        k_opt = find_optimal_k(one_cell_expr.T, max_clusters)
        clusterer = KMeans(n_clusters=k_opt, random_state=10, n_init=10, n_jobs=-1, tol=1e-3)
        cluster_labels = clusterer.fit_predict(one_cell_expr.T)

        one_cell_expr.columns = [cell + "_subtype_" + str(cluster_labels[i] + 1) + "." + str(i) for i in range(len(cluster_labels))]
        out_df[list(one_cell_expr)] = one_cell_expr
    out_df = out_df.rename_axis("gene")
    return out_df


def get_cluster_names(cells, cell_expr):
    """
    If the clustered dataset is imported from an external file, this function gets the number of clusters of
    each cell-type in the "cell" list, and assign name to each cluster in the format of cellName_subtype_clusterNum

    Args:
        cells: list of strings - Specifies what type of cells we have
        cell_expr: (n_genes, n_samples) pd DataFrame - Transposed DataFrame of cell expression
    Returns:
        cluster_name_list: list of strings - list of all cells and clusters in the file in the format of CD4_subtype_1.
    """
    import pandas as pd
    import numpy as np

    cluster_name_list = []
    for cell in cells:
        cell_type_list = list(cell_expr.loc[:, cell_expr.columns.str.contains(cell)])

        assert(len(cell_type_list) != 0), cell + " is not present in the file"

        cell_type_list = [cell_name[:cell_name.find(".")] for cell_name in cell_type_list]
        num_cell_cluster = len(np.unique(cell_type_list))
        for i in range(num_cell_cluster):
            cluster_name_list.append(cell + "_subtype_" + str(i+1))

    return cluster_name_list


def get_mean_of_each_cluster(cells, cell_expr):
    """
    Args:
        cells: list of strings - Specifies what type of cells we have
        cell_expr: (n_genes, n_samples) pd DataFrame - Transposed DataFrame of cell expression
    Returns:
        mean_df: (n_genes, n_clusters) pd Dataframe - mean value of each cluster across all samples of that cluster.
    """
    import pandas as pd

    cluster_name_list = get_cluster_names(cells, cell_expr)
    mean_df = pd.DataFrame(index=cell_expr.index)
    for cluster_name in cluster_name_list:
        one_cell_mean = cell_expr.loc[:, cell_expr.columns.str.contains(cluster_name)].mean(axis=1)
        mean_df[cluster_name] = one_cell_mean
    return mean_df


def detect_highly_expr_genes(one_cluster, other_clusters):
    """
    Args:
        one_cluster: (n_genes, n_samples) pd DataFrame- A single cluster with its samples
        other_clusters: (n_genes, n_samples) pd DataFrame - Other clusters!
    Returns:
        list of strings - returns the first 100 genes that are highly expressed in one_cluster, but lowly expressed in
        other_clusters
    """
    import pandas as pd
    import numpy as np

    one_cluster = np.log2(one_cluster + 1)
    other_clusters = np.log2(other_clusters + 1)
    one_cluster_min = one_cluster.min(axis=1)
    other_clusters_mean = other_clusters.mean(axis=1)
    diff = one_cluster_min - other_clusters_mean
    sub_df = diff[diff > 1]

    sub_df = sub_df.sort_values(ascending=False)
    if len(sub_df) > 100:
        sub_df = sub_df.iloc[:100]
    return list(sub_df.index)


def get_differentially_expr_genes(cells, cell_expr, mean_df, p_val_thresh=0.05):
    """
    Args:
        cells: list of strings - Specifies what type of cells we have
        cell_expr: (n_genes, n_samples) pd DataFrame - Transposed DataFrame of cell expression
        mean_df: (n_genes, n_clusters) pd Dataframe - mean value of each cluster across all samples of that cluster.
        p_val_thresh: float - Threshold for adjuster p value
    Returns:
        signature matrix - pd DataFrame
    """
    import pandas as pd
    import numpy as np
    import scipy.stats as stats
    import statsmodels.stats.multitest as multi

    cluster_name_list = get_cluster_names(cells, cell_expr)
    final_selected_genes = []
    for cluster_name in cluster_name_list:
        one_cell_expr = cell_expr.loc[:, cell_expr.columns.str.contains(cluster_name)]
        other_cell_exp = cell_expr.loc[:, ~cell_expr.columns.str.contains(cluster_name)]

        _, p_val = stats.ttest_ind(one_cell_expr, other_cell_exp, axis=1, equal_var=False)
        _, q_val, _, _ = multi.multipletests(p_val)
        one_cell_expr = one_cell_expr[q_val < p_val_thresh]
        other_cell_exp = other_cell_exp[q_val < p_val_thresh]
        selected_genes = detect_highly_expr_genes(one_cell_expr, other_cell_exp)
        final_selected_genes += selected_genes
    final_selected_genes = np.unique(final_selected_genes)

    return mean_df.loc[final_selected_genes]
