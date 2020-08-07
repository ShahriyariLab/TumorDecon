# batch_correction.py - functions originally authored by Shaya Akbarinejad

def batch_correct_datasets(data_loc, cell_file_dict, outfile="batch_corrected_datasets.txt"):
    """
    Performs batch correction on raw data files. All files must be stored in the same folder.
    Input:
        - data_loc: string. Path to folder that holds all data files
        - cell_file_dict: Dictionary. Keys = Cell Types, Values = paths to all raw files
            containing data for the given cell type
        - outfile: string. Name to save the batch corrected datasets file under.
            Default is "batch_corrected_datasets.txt"
    Output:
        - Single file 'batch_corrected_datasets.txt' containing all batch corrected
            data for all cell types in the input dictionary. Has shape (n_genes, n_samples)
    """
    import pandas as pd

    print("Performing batch correction...")

    batch_corrected_list = []
    for cell in cell_file_dict:
        file_list = cell_file_dict[cell]
        if len(file_list) > 1:
            batch_corrected_list.append(remove_batch_effect(file_list, cell, data_loc))
        else:
            file_name = file_list[0]
            data = pd.read_csv(data_loc+file_name, sep="\t", header=0)
            data.set_index(list(data)[0], inplace=True)
            cell_type_data = data.loc[:, data.columns.str.contains(cell)]
            assert(cell_type_data.shape[1] != 0), cell + " is not present in " + file_name
            cell_type_data = cell_type_data.loc[~cell_type_data.index.duplicated(keep='first')]
            batch_corrected_list.append(cell_type_data)

    batch_corrected_df = pd.concat(batch_corrected_list, axis=1, sort=True)
    batch_corrected_df.dropna(inplace=True)

    print("Saving batch corrected datasets to %s..." % outfile)
    batch_corrected_df.to_csv(outfile, sep="\t")
    # return batch_corrected_df


def differentiate_same_col_names(input_df):
    """
    There might be multiple columns of the same name. This results in Combat function to raise an error. This function
    makes sure that each column has a unique name.

    Args:
        input_df: pd DataFrame with repetitive column names
    Returns:
        pd DataFrame with unique column names, the same size as input_df
    """
    col_names = list(input_df)
    col_list = []
    col_dict = {}
    for col_name in col_names:
        if col_name in col_dict:
            col_dict[col_name] += 1
            col_list.append(col_name + "_" + str(col_dict[col_name]))
        else:
            col_dict[col_name] = 0
            col_list.append(col_name)
    input_df.columns = col_list
    return input_df


def remove_batch_effect(files_list, cell_type, address):
    """
    Args:
        files_list: List of strings: List of names of files that include cell_type
        cell_type: string:  Name of cell type we are looking for
    Returns:
        batch_corrected_data: (n_genes, n_samples) pd Dataframe:  Combined gene expression of cell types with removed
        batch effect
    """
    import pandas as pd
    from combat.pycombat import pycombat

    batch = pd.Series()
    data_dict = []
    for (i, file_name) in enumerate(files_list):
        data = pd.read_csv(address+file_name, sep="\t", header=0)
        data.set_index(list(data)[0], inplace=True)
        cell_type_data = data.loc[:, data.columns.str.contains(cell_type)]
        assert(cell_type_data.shape[1] != 0), cell_type + " is not present in " + file_name

        temp_batch = pd.Series([i for _ in range(len(list(cell_type_data)))], index=list(cell_type_data))
        batch = batch.append(temp_batch)
        cell_type_data = cell_type_data.loc[~cell_type_data.index.duplicated(keep='first')]
        data_dict.append(cell_type_data)
    data_combined = pd.concat(data_dict, axis=1, sort=True)
    data_combined = data_combined.apply(lambda row: row.fillna(row.mean()), axis=1)
    data_combined.dropna(inplace=True)
    var = data_combined.var(axis=1)
    data_combined = data_combined[var > 0]
    data_combined = differentiate_same_col_names(data_combined)
    batch.index = list(data_combined)
    batch_corrected_data = pycombat(data_combined, batch)
    return batch_corrected_data
