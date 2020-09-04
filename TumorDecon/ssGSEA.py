# ssGSEA.py

def ssGSEA(mix_data, up_genes, alpha=1.0, ties_method="average"):
    """
    Performs ssGSEA on a single sample for a single cell type and returns the enrichment score
    Inputs:
        - mix_data: Series of mixture RNA gene expression values, indexed on gene Hugo_Symbol, for a single patient
        - up_genes: a list of hugo identifiers referencing the up-regulated gene set to use to calculate enrichment score
        - alpha: Weight used in ssGSEA method
        - ties_method: see pandas.DataFrame.rank(). How to treat ties
    Output:
        - ssGSEA enrichment Score (float)
    """

    import pandas as pd
    import numpy as np

    mix_orig = mix_data.copy()#.reset_index(inplace=False)

    # Sort mixture data in descending order:
    mix_data = mix_data.rank(axis=0, method=ties_method).astype('int64') #, method='first', ascending=False)
    mix_data = mix_data.sort_values(ascending = False)

    # Sum all genes present in up_genes (signature):
    sm = 0
    # Count number of genes encountered that are not present in up_genes (signature)
    not_count = 0

    # Iterate through and compute cumulative sums:
    P_G = []
    P_NG = []
    for i, val in enumerate(mix_data):
        gene = mix_data.index[i]
        if gene in up_genes:
            sm += pow(np.abs(val), alpha)
        else:
            not_count += 1.

        P_G.append(sm)
        P_NG.append(not_count)

    # Calculate denominator for P_G using only genes in BOTH mix data and up_genes
    mix_sig_sum = np.sum(np.power(np.abs(mix_data.loc[mix_data.index.intersection(up_genes)]), alpha))

    P_G = np.array(P_G) / mix_sig_sum
    P_NG = np.array(P_NG) / not_count # denominator N_N_G_diff

    rank_diff = np.subtract(P_G, P_NG)
    rank = np.sum(rank_diff)
    return rank


def ssGSEA_main(rna_df, up_genes=None, patient_IDs='ALL', args={}):
    """
    This function does the following:
        - parses the dictionary 'args' for the arguments to pass on to the ssGSEA method.
        - If up_genes not passed in to function, uses the list of up-regulated genes provided in the ssGSEA paper
        - Runs ssGSEA() for each patient specified in 'patient_IDs' argument
        - Combines the resulting scores/ranks into a pandas dataframe (num_celltypes x num_patients)
    Inputs:
        - rna_df: pandas df of rna gene expression data.
            Rows are genes (indexed by 'Hugo_Symbol') and columns are patients
        - up_genes: dictionary. Keys are cell types, vaues are a list of up-regulated genes (hugo symbols) for
            the given cell type
        - patient_IDs: list of patient IDs to run ssGSEA on.
            Alternatively, can use the string 'ALL' to run for all patients
        - args: dictionary containing any of the following:
            - alpha: see ssGSEA()
            - ties_method: see ssGSEA()
            - print_progress: whether to print patient ID as ssGSEA iterates through
            - norm: whether or not to normalize enrichment scores by using the entire data set, as indicated
                by Barbie et al., 2009, online methods, pg. 2
    Outputs:
        - scores: pandas df. Contains ssGSEA enrichment scores for each cell type declared in 'up_genes'
            dictionary, for each patient in 'patient_IDs' list.
            Rows are indexed by cell type, columns are patient IDs
    """
    import pandas as pd
    import numpy as np
    from collections import Mapping
    from .data_utils import read_ssGSEA_up_genes

    # Read in optional arguments, or set them as default
    # Assert values are of the right data type when passed to ssGSEA() function.

    # Patient_ID must be 'ALL' or an array of specific patient IDs.
    if patient_IDs == 'ALL':
        patient_list = rna_df.columns
    elif not isinstance(patient_IDs, type([])):
        raise ValueError("patient_IDs should be either 'ALL', or an array of IDs (not a single ID)")
    else:
        patient_list = patient_IDs

    if 'alpha' in args.keys():
        alpha = args['alpha']
    else:
        alpha = 1.0
    if 'ties_method' in args.keys():
        ties_method = args['ties_method']
    else:
        ties_method = 'average'
    if 'print_progress' in args.keys():
        print_progress = args['print_progress']
    else:
        print_progress = False
    if 'norm' in args.keys():
        norm = args['norm']
    else:
        norm = False


    # Check if up_genes argument was passed in:
    if up_genes is not None:
        if not isinstance(up_genes, Mapping):
            raise ValueError("up_genes argument must be a dictionary")
    else:
        # use the up-regulated genes list that was in the ssGSEA paper:
        up_genes = read_ssGSEA_up_genes()

    # Get ranks for each cell type, for each patient:
    print("Running ssGSEA...")
    scores = pd.DataFrame()
    scores['Patient_ID'] = list(up_genes.keys())
    scores = scores.set_index(['Patient_ID'])
    for patient in patient_list:
        if patient in rna_df.columns:
            if print_progress:
                print(patient)
            patient_scores = [].copy()
            for cell_type in up_genes:
                score = ssGSEA(rna_df[patient], up_genes[cell_type], alpha=alpha, ties_method=ties_method)
                patient_scores.append(score)
            scores[patient] = patient_scores
        else:
            raise ValueError("patient_ID ({!r}) not present in rna dataframe".format(patient))

    if norm:
        range = np.max(np.array(scores)) - np.min(np.array(scores))
        scores = scores / range

    scores = scores.transpose()
    return scores
