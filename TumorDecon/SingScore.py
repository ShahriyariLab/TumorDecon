# singscore.py

def SingScore_main(rna_df, up_genes=None, down_genes=None, patient_IDs='ALL', args={}):
    """
    * If down_genes and up_genes both defined, will run PySingScore's implementation of bidirectional singscore.
    * If just up_genes defined, will run PySingScore's implementation of unidirectional singscore.

    Inputs:
        - rna_df: pandas df of rna gene expression data.
            Rows are genes (indexed by 'Hugo_Symbol') and columns are patients
        - up_genes: dictionary. Keys are cell types, vaues are a list of up-regulated genes (hugo symbols) for
            the given cell type
        - down_genes: dictionary. Keys are cell types, vaues are a list of down-regulated genes (hugo symbols) for
            the given cell type
        - patient_IDs: list of patient IDs to run SingScore on.
            Alternatively, can use the string 'ALL' to run for all patients
        - args is empty, but included for flexibility / eventually making the methods subclasses of each other
    Outputs:
        - scores: pandas df. Contains SingScore scores for each cell type declared in 'up_genes'
            dictionary, for each patient in 'patient_IDs' list.
            Rows are indexed by cell type, columns are patient IDs
    """
    from singscore import singscore
    import pandas as pd
    from .data_utils import read_ssGSEA_up_genes
    from collections import Mapping

    # Select a patient / list of patients to solve for their cell type frequencies:
        # Patient_ID must be 'ALL' or an array of specific patient IDs.
    if patient_IDs == 'ALL':
        patient_list = rna_df.columns
    elif not isinstance(patient_IDs, type([])):
        raise ValueError("patient_IDs should be either 'ALL', or an array of IDs (not a single ID)")
    else:
        patient_list = patient_IDs

    if up_genes is None:
        up_genes = read_ssGSEA_up_genes()
    else:
        if not isinstance(up_genes, Mapping):
            raise ValueError("up_genes argument must be a dictionary")

    if down_genes is not None:
        if not isinstance(up_genes, Mapping):
            raise ValueError("up_genes argument must be a dictionary")
        elif set(up_genes.keys()) != set(down_genes.keys()):
            raise KeyError("Keys (cell type) of up_genes and down_genes must match")

    # Select only the patients in the patient_list:
    sample_df = rna_df[patient_list]

    # Use PySingScore package to generate scores for each sample in sample_df:
    print("Running SingScore...")
    scores = pd.DataFrame()
    for cell_type in up_genes:
        if down_genes is None:
            score = singscore.score(up_gene=up_genes[cell_type], sample=sample_df)
        else:
            score = singscore.score(up_gene=up_genes[cell_type], down_gene=up_genes[cell_type], sample=sample_df)
        score.rename(columns={'total_score': cell_type}, inplace=True)
        scores = pd.concat([scores, score], axis=1)
    scores.index.name='Patient_ID'
    
    return scores
