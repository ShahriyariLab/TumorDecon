# cibersort.py

def cibersort(rna_sample, sig_df, nu=0.5, C=1.0, kernel='linear', shrinking=True):
    """
    Uses NuSVR from sklearn to solve for the cell type frequencies
    Inputs:
        - rna_sample: pandas series. A single sample of rna expression values
        - sig_df: pandas df of Signature gene expression values for given cell types.
            Rows are genes (indexed by 'Hugo_Symbol') and columns are cell types
        - nu: see sklearn's NuSVR
        - C: see sklearn's NuSVR
        - kernel: see sklearn's NuSVR
        - shrinking: see sklearn's NuSVR
    Outputs:
        - weights: NuSVR solution vector for the given sample.
            (Negative values in the solution vector are set to 0 for interpretation as cell type frequencies)
    """
    from sklearn.svm import NuSVR
    import numpy as np
    from sklearn.model_selection import GridSearchCV

    # If a numerical of nu not explicitly specified, use gridsearch to find the best nu:
    if nu == 'best':
        gridsearch = GridSearchCV(NuSVR(C=C, kernel=kernel, max_iter=-1, shrinking=shrinking), cv=5, param_grid={"nu": [0.25, 0.5, 0.75]}, scoring='neg_mean_squared_error', refit=True)
        gridsearch.fit(sig_df, rna_sample)
        nu = gridsearch.best_params_['nu']

    # Fit nuSVR with best (or specified) value of nu:
    clf = NuSVR(nu=nu, C=C, kernel=kernel, max_iter=-1, shrinking=shrinking, tol=1e-3)
    clf.fit(sig_df, rna_sample)

    # Replace negative "frequencies" with 0:
    weights = np.array(clf.coef_)[0] # equivalent to np.matmul(np.array(clf.dual_coef_), np.array(clf.support_vectors_))[0]
    weights[weights<0] = 0
    # Sum to 1 contraint:
    weights = weights / np.sum(weights)

    return weights


def cibersort_main(rna_df, sig_df, patient_IDs='ALL', args={}):
    """
    This function does the following:
        - parses the dictionary 'args' for the arguments to pass on to the CIBERSORT method.
        - eliminates genes from rna_df and sig_df that are not present in both data sets
        - Runs cibersort() for each patient specified in patient_IDs' argument
        - Combines the resulting frequencies into a pandas dataframe (num_celltypes x num_patients)
    Inputs:
        - rna_df: pandas df of rna gene expression data.
            Rows are genes (indexed by 'Hugo_Symbol') and columns are patients
        - sig_df: pandas df of Signature gene expression values for given cell types.
            Rows are genes (indexed by 'Hugo_Symbol') and columns are cell types
        - patient_IDs: list of patient IDs to run cibersort for.
            Alternatively, can use the string 'ALL' to run for all patients
        - args: dictionary containing any of the following:
            - print_progress: whether to print patient ID as cibersort iterates through
            - scaling: string, must be either 'None', 'zscore', or 'minmax'. Determines how to scale the mixture data and signature matrix before applying CIBERSORT
            - scaling_axis: 0 or 1. Whether to scale mixture data and signature matrix by normalizing each column (patient/celltype) individually (scaling_axis=0) or each row (gene) individually (scaling_axis=1).
            - nu: see sklearn's NuSVR. If nu='best', will use the nu in {0.25, 0.5, 0.75} that minimizes root mean square b/w data and prediction. If nu a float, will compute cibersort for that specific nu
            - C: see sklearn's NuSVR
            - kernel:: see sklearn's NuSVR
            - shrinking: see sklearn's NuSVR
    Outputs:
        - cell_freqs: pandas df. Contains cell type frequencies for each patient in 'patient_IDs' list.
            Rows are indexed by cell type, columns are patient IDs
    """

    import pandas as pd
    import numpy as np
    from .data_utils import keep_common_genes
    from .data_utils import df_normalization

    # Patient_ID must be 'ALL' or an array of specific patient IDs.
    if patient_IDs == 'ALL':
        patient_list = rna_df.columns
    elif not isinstance(patient_IDs, type([])):
        raise ValueError("patient_IDs should be either 'ALL', or an array of IDs (not a single ID)")
    else:
        patient_list = patient_IDs

    if 'print_progress' in args.keys():
        print_progress = args['print_progress']
    else:
        print_progress = False

    if 'scaling' in args.keys():
        scaling = args['scaling']
        if scaling not in ['None', 'none', 'zscore', 'minmax', 'r-zscore']:
            raise ValueError("scaling ({!r}) must be set to 'none', 'zscore' or 'minmax'".format(scaling))
        else:
            scaling = args['scaling']
    else: # default
        scaling = "minmax"

    if 'scaling_axis' in args.keys():
        scaling_axis = args['scaling_axis']
        if scaling_axis not in [0, 1]:
            raise ValueError("scaling_axis ({!r}) must be 0 or 1".format(scaling_axis))
    else:
        scaling_axis = 0

    if 'nu' in args.keys():
        nu = args['nu']
        if isinstance(nu, str):
            if nu != 'best':
                raise ValueError("nu ({!r}) must be either a float, or 'best'".format(nu))
    else:
        nu = 'best'

    if 'C' in args.keys():
        C = args['C']
    else:
        C = 1.0
    if 'kernel' in args.keys():
        kernel = args['kernel']
    else:
        kernel = 'linear'
    if 'shrinking' in args.keys():
        shrinking = args['shrinking']
    else:
        shrinking = True

    # eliminate genes not present in both rna and sig dfs, and ensure they list genes in the same order:
    rna_df, sig_df = keep_common_genes(rna_df, sig_df)

    # Scale data:
    if scaling in ['zscore', 'minmax', 'r-zscore']:
        sig_df = df_normalization(sig_df, scaling=scaling, axis=scaling_axis)
        rna_df = df_normalization(rna_df, scaling=scaling, axis=scaling_axis)

    # For each patient, run cibersort to get cell type frequencies, and save results to pandas df:
    cell_freqs_df = pd.DataFrame(columns = patient_list)
    cell_freqs_df['Patient_ID'] = sig_df.columns
    cell_freqs_df = cell_freqs_df.set_index(['Patient_ID'])

    print("Running CiberSort...")
    for patient in patient_list:
        if patient in rna_df.columns:
            if print_progress:
                print(patient)
            cell_freqs_df[patient] = cibersort(rna_df[patient], sig_df, nu=nu, C=C, kernel=kernel, shrinking=shrinking)
        else:
            raise ValueError("patient_ID ({!r}) not present in rna dataframe".format(patient))
            return
    cell_freqs_df = cell_freqs_df.transpose()
    return cell_freqs_df
