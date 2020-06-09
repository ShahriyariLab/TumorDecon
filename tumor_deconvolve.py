# Import our functions:
from TumorDecon.data_utils import *
from TumorDecon.cibersort import *
from TumorDecon.DeconRNASeq import *
from TumorDecon.ssGSEA import *
from TumorDecon.SingScore import *

def tumor_deconvolve(mixture_data, method, patient_IDs='ALL', cell_signatures=None, up_genes=None, down_genes=None, args={}):
    """
        For each patient, given mixture data, cell signatures, and a deconvolution
            method, return frequencies/ranks for each cell type
        Inputs:
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
        Outputs:
            - x: dataframe of frequencies/scores for each specified patient
                and each cell type in the signature matrix, as solved for using
                the specified method
    """
    import pandas as pd
    from collections.abc import Mapping

    if not isinstance(args, Mapping):
        raise TypeError('args ({!r}) must be a dict '.format(args))

    # Check data:
    if isinstance(mixture_data, pd.DataFrame):
        # TO DO: Add check for proper indexing
        pass
    else:
        raise TypeError("mixture_data must be a pandas dataframe")

    if cell_signatures is not None:
        if isinstance(cell_signatures, pd.DataFrame):
            # TO DO: Make sure index matches that of mixture_data
            pass
        else:
            raise TypeError("cell_signatures must be a pandas dataframe")

    if up_genes is not None:
        if isinstance(up_genes, Mapping):
            # TO DO: Make sure index matches that of mixture_data
            pass
        else:
            raise TypeError("up_genes must be a dictionary")

    if down_genes is not None:
        if isinstance(down_genes, Mapping):
            # TO DO: Make sure index matches that of mixture_data
            pass
        else:
            raise TypeError("down_genes must be a dictionary")

    # don't need to check patient_IDs/up_genes/down_genes since check them within
        # each method

    # Apply the method:
    if method in ['CIBERSORT', 'Cibersort', 'cibersort']:
        if cell_signatures is not None:
            x = cibersort_main(mixture_data, cell_signatures, patient_IDs, args)
        else:
            raise ValueError("'cell_signatures' argument required for CIBERSORT method")

    elif method in ['DeconRNASeq', 'deconrnaseq', 'DeconRNAseq', 'deconRNAseq']:
        if cell_signatures is not None:
            x = DeconRNASeq_main(mixture_data, cell_signatures, patient_IDs, args)
        else:
            raise ValueError("'cell_signatures' argument required for DeconRNASeq method")

    elif method in ['ssGSEA', 'ssgsea', 'SSGSEA']:
        if up_genes is not None:
            x = ssGSEA_main(mixture_data, up_genes, patient_IDs, args)
        else:
            raise ValueError("'up_genes' argument required for ssGSEA method")

    elif method in ['SingScore', 'singscore', 'Singscore']:
        if up_genes is not None:
            x = SingScore_main(mixture_data, up_genes, down_genes, patient_IDs)
        else:
            raise ValueError("'up_genes' argument required for SingScore method")

    else:
        raise ValueError("Deconvolution Method must be one of 'CIBERSORT', 'DeconRNASeq', 'ssGSEA', or 'SingScore'")

    return x
