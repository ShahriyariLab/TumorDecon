# DeconRNASeq.py

def compute_least_squares(x, A, b):
    '''
    Inputs:
        x = vector of length n (frequency vector of length n_celltypes)
    Parameters:
        A = matrix of size p x n (signature matrix of size n_genes x n_celltypes)
        b = vector of size p (mixture data vector of length n_genes)
    Output:
        value of ||Ax-b||^2
    '''
    import numpy as np
    return np.linalg.norm(np.matmul(A,x)-b)


def compute_ridge(x, A, b, reg=1.0):
    '''
    Inputs:
        x = vector of length n (frequency vector of length n_celltypes)
    Parameters:
        A = matrix of size p x n (signature matrix of size n_genes x n_celltypes)
        b = vector of size p (mixture data vector of length n_genes)
        reg = regulation constant (lamba)
    Output:
        value of ||Ax-b||^2 + ||x||^2
    '''
    import numpy as np
    return np.linalg.norm(np.matmul(A,x)-b) + reg*np.linalg.norm(x)


def compute_lasso(x, A, b, reg=1.0):
    '''
    Inputs:
        x = vector of length n (frequency vector of length n_celltypes)
    Parameters:
        A = matrix of size p x n (signature matrix of size n_genes x n_celltypes)
        b = vector of size p (mixture data vector of length n_genes)
        reg = regulation constant (lamba)
    Output:
        value of ||Ax-b||^2 + Sum_i(|x_i|)
    '''
    import numpy as np
    return np.linalg.norm(np.matmul(A,x)-b) + reg*np.sum(np.abs(x))



def DeconRNASeq(S, m, formulation="qp", reg_constant=1.0, print_results=False, label=''):
    """
    Solves for the vector x that minimizes ||Sx-m||^2 under the additional constraints:
        Sum(x_i) = 1
        x_i >=0 for all i
    Depends on the function 'compute_least_squares(x, A, b)'
    Inputs:
        - S: numpy array (n_genes x n_celltypes), Signature matrix
        - m: numpy array (n_genes x 1), Mixture data vector for a given patient
        - formulation: string, must be either 'qp', 'ridge', or 'lasso'. Determines how to formulate the optimization problem:
            'qp': solve a QP (quadratic programming) problem min_x||Ax-b||^2 with strict condition x.sum()=1, xi>=0
            'ridge': solve a ridge regression problem min_x{||Ax-b||^2 + ||x||^2} with condition xi>=0
            'lasso': solve a lasso regression problem min_x{||Ax-b||^2 + Sum_i{|x_i|} with condition xi>=0
        - reg_constant: float, regularization constant for lasso/ridge regression
        - print_result: boolean, whether or not to print the solver (scipy.optimize) results
        - label: patient ID label. Only used if print_results = True
    Outputs:
        - xopt: optimal frequency vector found by scipy.optimize
        - (if print_results) prints the scipy.optimize output to standard output
    """
    import scipy.optimize as opt
    import numpy as np
    import warnings

    # Initialize frequencies with random numbers between 0 and 1:
    x0 = np.random.rand(S.shape[1])

    if formulation == 'qp':
        # Set bounds on x_i values to be between 0 and 1
        bnds = ((0., 1.),) * S.shape[1]

        # Define equality constraint that Sum(x_i) = 1
        cons = ({'type': 'eq',
                 'fun' : lambda x: np.array([np.sum(x)-1])})

        # Find the x that minimizes this least squares problem:
        soln = opt.minimize(compute_least_squares, x0, bounds=bnds, constraints=cons, args=(S,m))

    else:
        # Set bounds on x_i to be positive (choose a "large enough" value for the upper bound)
        bnds = ((0., 1000.),) * S.shape[1]

        # Find the x that minimizes this unconstrained ridge/lasso Regression problem:
        if formulation == 'ridge':
            soln = opt.minimize(compute_ridge, x0, bounds=bnds, constraints=None, args=(S,m,reg_constant))
        elif formulation == 'lasso':
            soln = opt.minimize(compute_lasso, x0, bounds=bnds, constraints=None, args=(S,m,reg_constant))

    if print_results:
        print()
        print(label)
        print(soln)
        print()
    else:
        if soln['success'] == False:
            warnings.warn("DeconRNASeq optimization did not terminate successfully. For more details: re-run with optional argument 'print_results:True'", category=UserWarning)

    xopt = soln['x']

    return xopt



def DeconRNASeq_main(rna_df, sig_df, patient_IDs='ALL', args={}):
    """
    This function does the following:
        - parses the dictionary 'args' for the arguments to pass on to the DeconRNASeq method.
        - eliminates genes from rna_df and sig_df that are not present in both data sets
        - Runs DeconRNASeq() for each patient specified in patient_IDs' argument
        - Combines the resulting frequencies into a pandas dataframe (num_celltypes x num_patients)
    Inputs:
        - rna_df: pandas df of rna gene expression data.
            Rows are genes (indexed by 'Hugo_Symbol') and columns are patients
        - sig_df: pandas df of Signature gene expression values for given cell types.
            Rows are genes (indexed by 'Hugo_Symbol') and columns are cell types
        - patient_IDs: list of patient IDs to run DeconRNASeq for.
                Alternatively, can use the string 'ALL' to run for all patients
        - args: dictionary containing any of the following:
            - check_sig: boolean, whether or not to check the condition number of the signature matrix
            - scaling: string, must be either 'None', 'zscore', or 'minmax'. Determines how to scale the signature matrix and mixture data before solving for optimal x
            - scaling_axis: 0 or 1. Whether to scale mixture data and signature matrix by normalizing each column (celltype/patient) separately (scaling_axis=0) or each row (gene) separately (scaling_axis=1).
            - formulation: see DeconRNASeq()
            - reg_constant: see DeconRNASeq()
            - print_result: see DeconRNASeq()
    Outputs:
        - cell_freqs: pandas df. Contains cell type frequencies for each patient in 'patient_IDs' list.
            Rows are indexed by cell type, columns are patient IDs
    """
    import pandas as pd
    import numpy as np
    from .data_utils import keep_common_genes
    from .data_utils import df_normalization

    # Read in optional arguments, or set them as default
    # Assert values are of the right data type when passed to DeconRNASeq() function.

    # formulation must be 'qp', 'ridge', or 'lasso'
    if 'formulation' in args.keys():
        formulation = args['formulation']
        if formulation not in ['qp','ridge','lasso']:
            raise ValueError("Formulation ({!r}) must be set to 'qp', 'ridge', or 'lasso'".format(formulation))
    else:
        formulation = 'qp'
    # reg_constant must be a double
    if 'reg_constant' in args.keys():
        reg_constant = args['reg_constant']
    else:
        reg_constant = 1.0
    if 'check_sig' in args.keys():
        check_sig = args['check_sig']
        if not isinstance(check_sig, bool):
            raise ValueError("check_sig ({!r}) must be a boolean variable".format(check_sig))
    else:
        check_sig = False
    if 'scaling' in args.keys():
        scaling = args['scaling']
        if scaling not in ['None', 'none', 'zscore', 'minmax', 'r-zscore']:
            raise ValueError("scaling ({!r}) must be set to 'none', 'zscore' or 'minmax'".format(scaling))
    else:
        scaling = 'minmax'
    if 'scaling_axis' in args.keys():
        scaling_axis = args['scaling_axis']
        if scaling_axis not in [0, 1]:
            raise ValueError("scaling_axis ({!r}) must be 0 or 1".format(scaling_axis))
    else:
        scaling_axis = 0
    if 'print_results' in args.keys():
        print_results = args['print_results']
        if not isinstance(print_results, bool):
            raise ValueError("print_results ({!r}) must be a boolean variable".format(print_results))
    else:
        print_results = False


    # eliminate genes not present in both rna and sig dfs, and ensure they are in the same order:
    rna_df, sig_df = keep_common_genes(rna_df, sig_df)

    # Scale Data:
    if scaling in ['zscore', 'minmax', 'r-zscore']:
        # R implementation uses zscore scaling.
        sig_df = df_normalization(sig_df, scaling=scaling, axis=scaling_axis)
        rna_df = df_normalization(rna_df, scaling=scaling, axis=scaling_axis)

    # Convert signature to numpy array
    Sig = np.array(sig_df)

    # Check the condition number of the signature matrix:
    if check_sig:
        print("Condition number of signature matrix =", np.linalg.cond(Sig))

    # Select a patient / list of patients to solve for their cell type frequencies:
         # Patient_ID must be 'ALL' or an array of specific patient IDs.
    if patient_IDs == 'ALL':
        patient_list = rna_df.columns
    elif not isinstance(patient_IDs, type([])):
        raise ValueError("patient_IDs should be either 'ALL', or an array of IDs (not a single ID)")
    else:
        patient_list = patient_IDs

    # For each patient, run DeconRNASeq to get cell type frequencies, and save results to pandas df:
    print("Running DeconRNASeq...")
    cell_freqs_df = pd.DataFrame()
    cell_freqs_df['Patient_ID'] = sig_df.columns
    cell_freqs_df = cell_freqs_df.set_index(['Patient_ID'])
    for patient in patient_list:
        if patient in rna_df.columns:
            Mix = np.array(rna_df[patient])
            cell_freqs_df[patient] = DeconRNASeq(Sig, Mix, formulation=formulation, reg_constant=reg_constant, print_results=print_results, label=patient)
        else:
            raise ValueError("patient_ID ({!r}) not present in rna dataframe".format(patient))

    cell_freqs_df = cell_freqs_df.transpose()

    return cell_freqs_df
