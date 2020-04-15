# GUI for running tumorDecon

# Add location of TumorDecon repo temporary python path:
    # (won't need once package properly installed in a location included in python's path)
import sys
sys.path.insert(0,'../../')

# Import packages
import TumorDecon as td

from tkinter import *
import tkinter as tk

from tkinter import filedialog as fd
from tkinter import ttk
from tkinter import messagebox as mb

import os

# Provided from geeksforgeeks website
from CollapsiblePane import CollapsiblePane


def browse_for_file():
    name = fd.askopenfilename()
    # global rna_df
    # rna_df = td.read_rna_file(name)
    # # Remove nans:
    # rna_df.dropna(axis=0, inplace=True)
    rna_file_selection["text"] = name

def get_custom_signature():
    global path_to_custom_sig
    path_to_custom_sig = fd.askopenfilename()
    sig_selection["text"] = path_to_custom_sig

def get_custom_up_geneset():
    global path_to_up_genes
    path_to_up_genes = fd.askopenfilename()
    upgenes_selection["text"] = path_to_up_genes

def get_custom_down_geneset():
    global path_to_up_genes
    path_to_down_genes = fd.askopenfilename()
    downgenes_selection["text"] = path_to_down_genes

def show_method_options(meth):
    # Reset the other parameters:
    ####
    if meth in ["CIBERSORT (Newman et al., 2015)", "DeconRNASeq (Gong and Szustakowski, 2013)"]:
        sig_label.grid()
        sig_dropdown.grid()
        sig_selection.grid()
        norm_label.grid()
        norm_dropdown.grid()

        upgene_label.grid_remove()
        upgene_dropdown.grid_remove()
        upgenes_selection.grid_remove()
        downgene_label.grid_remove()
        downgene_dropdown.grid_remove()
        downgenes_selection.grid_remove()
        alpha_label.grid_remove()
        alpha.grid_remove()
    else:
        sig_label.grid_remove()
        sig_dropdown.grid_remove()
        cpane2.grid_remove()
        sig_browse.grid_remove()
        sig_selection.grid_remove()
        norm_label.grid_remove()
        norm_dropdown.grid_remove()

        upgene_label.grid()
        upgene_dropdown.grid()
        upgenes_selection.grid()

        if meth == "ssGSEA (Barbie et al., 2009) (Senbabaoglu, 2016)":
            alpha_label.grid()
            alpha.grid()
            downgene_label.grid_remove()
            downgene_dropdown.grid_remove()
            downgenes_browse.grid_remove()
            downgenes_selection.grid_remove()
        else:
            alpha_label.grid_remove()
            alpha.grid_remove()
            downgene_label.grid()
            downgene_dropdown.grid()
            downgenes_selection.grid()

def show_rna_options(arg):
    if arg == 'Download from cBioPortal':
        cpane.grid_remove()
        rna_browse.grid_remove()
        rna_file_selection.grid_remove()
        rna_download_label.grid()
        xena_type_dropdown.grid_remove()
        cbio_type_dropdown.grid()
    elif arg == 'Download from UCSC Xena TCGA Hub':
        cpane.grid_remove()
        rna_browse.grid_remove()
        rna_file_selection.grid_remove()
        rna_download_label.grid()
        cbio_type_dropdown.grid_remove()
        xena_type_dropdown.grid()
    else:
        rna_download_label.grid_remove()
        cbio_type_dropdown.grid_remove()
        xena_type_dropdown.grid_remove()
        cpane.grid()
        rna_browse.grid()
        rna_file_selection.grid()

def show_sig_options(arg):
    if arg == 'LM22':
        cpane2.grid_remove()
        sig_browse.grid_remove()
        sig_selection.grid_remove()
    else:
        cpane2.grid()
        sig_browse.grid()
        sig_selection.grid()

def show_upgenes_options(arg):
    if arg == 'ssGSEA paper gene set':
        cpane3.grid_remove()
        upgenes_browse.grid_remove()
        upgenes_selection.grid_remove()
    else:
        cpane3.grid()
        upgenes_browse.grid()
        upgenes_selection.grid()

def show_downgenes_options(arg):
    if arg == 'None':
        downgenes_browse.grid_remove()
        downgenes_selection.grid_remove()
    else:
        downgenes_browse.grid()
        downgenes_selection.grid()

def get_save_location():
    global save_loc
    # Check to make sure not overriding a file...
    save_loc = fd.askdirectory()
    if os.path.exists(save_loc+"/results.csv"):
        file_exists = True
        i = 1
    else:
        save_loc = save_loc+"/results.csv"
        print_save_loc["text"] = "Results will be saved to "+save_loc
        file_exists = False
    while file_exists:
        if os.path.exists(save_loc+"/results_"+str(i)+".csv"):
            i += 1
        else:
            save_loc = save_loc+"/results_"+str(i)+".csv"
            print_save_loc["text"] = "Results will be saved to "+save_loc
            file_exists = False

# def update_status(new_text):

# def update_status():
#     # Get the current message
#     current_status = status["text"]
#     # If the message is "Working...", start over with "Working"
#     if current_status.endswith("..."):
#         current_status = "Working"
#     # If not, then just add a "." on the end
#     else:
#         current_status += "."
#     # Update the message
#     status["text"] = current_status
#     # After 1 second, update the status
#     root.after(1000, update_status)


def run_tumor_decon():
    # Read in RNA data:
    if r.get() == "Download from cBioPortal":
        rna_df = td.download_by_name('cbio', cbio_type.get())
    elif r.get() == "Download from UCSC Xena TCGA Hub":
        rna_df = td.download_by_name('xena', xena_type.get())
    else:
        rna_df = td.read_rna_file(rna_file_selection["text"])
    # Remove nans:
    rna_df.dropna(axis=0, inplace=True)

    method = 'cibersort'
    if "DeconRNASeq" in m.get():
        method = 'deconrnaseq'
    elif "ssGSEA" in m.get():
        method = 'ssgsea'
    elif "SingScore" in m.get():
        method = 'singscore'

    if method in ['cibersort','deconrnaseq']:
        if s.get() == "LM22":
            sig_df = td.read_lm22_file("../data/LM22.txt")
        else:
            sig_df = td.read_lm22_file(path_to_custom_sig)
        scaling = "None"
        if normalization.get() == "Z-Scaling":
            scaling = 'zscore'
        elif normalization.get() == 'MinMax Normalization':
            scaling = 'minmax'
        new_text = "Ran "+m.get()+" with sig matrix "+s.get()+" and "+scaling+" normalization"
        status["text"] = new_text
        solns = td.tumor_deconvolve(rna_df, method, patient_IDs='ALL', cell_signatures=sig_df, args={'scaling_axix':0, 'scaling':scaling})

    else: # ssgsea, singscore
        if u.get() == "ssGSEA paper gene set":
            up_genes = td.read_ssGSEA_up_genes(td.get_td_Home()+"data/Gene_sets.csv")
        else:
            up_genes = td.read_custom_geneset(path_to_up_genes)
        if method == 'ssgsea':
            new_text = "Ran "+m.get()+" with up genes "+u.get()
            status["text"] = new_text
            solns = td.tumor_deconvolve(rna_df, method, patient_IDs='ALL', up_genes=up_genes, args={'alpha':float(alpha.get()), 'print_progress':False})
        else: # singscore
            if d.get() != "None":
                down_genes = td.read_custom_geneset(path_to_down_genes)
                new_text = "Ran "+m.get()+" with up genes "+u.get()+" down genes "+d.get()
                status["text"] = new_text
                solns = td.tumor_deconvolve(rna_df, method, patient_IDs='ALL', up_genes=up_genes, down_genes=down_genes)
            else:
                new_text = "Ran "+m.get()+" with up genes "+u.get()
                status["text"] = new_text
                solns = td.tumor_deconvolve(rna_df, method, patient_IDs='ALL', up_genes=up_genes, down_genes=None)
    try:
        solns.to_csv(save_loc)
    except NameError:
        print("No save location specified")
    done_text = new_text+"\n Results saved to "+save_loc
    status["text"] = done_text

if __name__ == '__main__':
    root = Tk()
    # root=Frame(root,width=300,height=300)

    Label(root, text="Welcome to TumorDecon!", font="Times 24 bold").grid(row=0,column=0)
    # Upload RNA Data:
    Label(root,text="1. Choose RNA Data:", font="Times 16 bold").grid(row=1,column=0, sticky = W)
    r = StringVar()
    r.set("Download from cBioPortal")
    cbio_type, xena_type = StringVar(), StringVar()
    cbio_type.set("Acute Myeloid Leukemia")
    xena_type.set("Acute Myeloid Leukemia")
    rna_optionslist = ["Download from cBioPortal", "Download from UCSC Xena TCGA Hub", "Upload your own..."]
    rna_dropdown = OptionMenu(root, r, *rna_optionslist, command=show_rna_options)
    rna_dropdown.grid(row=2, sticky = W)
    rna_download_label = Label(root,text="Choose Data Set:")
    rna_download_label.grid(row=3, column=0, sticky=W)
    cbio_optionslist = ['Acute Myeloid Leukemia','Adrenocortical Carcinoma', 'Bladder Urothelial Carcinoma', 'Brain Lower Grade Glioma', 'Breast Invasive Carcinoma',
                 'Cervical Squamous Cell Carcinoma', 'Cholangiocarcinoma', 'Colorectal Adenocarcinoma', 'Diffuse Large B-Cell Lymphoma', 'Esophageal Adenocarcinoma',
                 'Glioblastoma Multiforme', 'Head and Neck Squamous Cell Carcinoma', 'Kidney Chromophobe', 'Kidney Renal Clear Cell Carcinoma', 'Kidney Renal Papillary Cell Carcinoma',
                 'Liver Hepatocellular Carcinoma', 'Lung Adenocarcinoma', 'Lung Squamous Cell Carcinoma', 'Mesothelioma', 'Ovarian Serous Cystadenocarcinoma', 'Pancreatic Adenocarcinoma',
                 'Pheochromocytoma and Paraganglioma', 'Prostate Adenocarcinoma', 'Sarcoma', 'Skin Cutaneous Melanoma', 'Stomach Adenocarcinoma', 'Testicular Germ Cell Tumors',
                 'Thymoma', 'Thyroid Carcinoma', 'Uterine Carcinosarcoma', 'Uterine Corpus Endometrial Carcinoma', 'Uveal Melanoma']
    xena_optionslist = ['Acute Myeloid Leukemia','Adrenocortical Cancer', 'Bile Duct Cancer', 'Bladder Cancer', 'Breast Cancer',
                 'Cervical Cancer', 'Colon and Rectal Cancer', 'Colon Cancer', 'Endometrioid Cancer', 'Esophageal Cancer',
                 'Glioblastoma', 'Head and Neck Cancer', 'Kidney Chromophobe', 'Kidney Clear Cell Carcinoma', 'Kidney Papillary Cell Carcinoma',
                 'Large B-cell Lymphoma', 'Liver Cancer', 'Lower Grade Glioma', 'Lower Grade Glioma and Glioblastoma', 'Lung Adenocarcinoma',
                 'Lung Cancer', 'Lung Squamous Cell Carcinoma', 'Melanoma', 'Mesothelioma', 'Ocular Melanomas', 'Ovarian Cancer',
                 'Pancreatic Cancer', 'Pheochromocytoma and Paraganglioma', 'Prostate Cancer', 'Rectal Cancer', 'Sarcoma', 'Stomach Cancer',
                 'Testicular Cancer', 'Thymoma', 'Thyroid Cancer', 'Uterine Carcinosarcoma']
    cbio_type_dropdown = OptionMenu(root, cbio_type, *cbio_optionslist)
    cbio_type_dropdown.grid(row=4)
    xena_type_dropdown = OptionMenu(root, xena_type, *xena_optionslist)
    xena_type_dropdown.grid(row=4)
    xena_type_dropdown.grid_remove()
    cpane = CollapsiblePane(root, 'Collapse Instructions', 'Instructions for uploading your own RNA data file:')
    cpane.grid(row=3,sticky = W)
    Label(cpane.frame,text="- Data file extension must be either .csv, .xlsx, or .txt").grid(row=4,sticky = W)
    Label(cpane.frame,text="- First column must contain a list of genes, referenced by their Hugo Symbol").grid(row=5,sticky = W)
    Label(cpane.frame,text="- Each remaining column should be the rna expression values of each gene for a patient").grid(row=6,sticky = W)
    Label(cpane.frame,text="- First row of file should be a header, containing column names (patient identifiers)").grid(row=7,sticky = W)
    rna_browse = Button(root, text='Browse...', command=browse_for_file)
    rna_browse.grid(row=8,column=0)
    rna_file_selection = tk.Label(root, text=" ")
    rna_file_selection.grid(row=9, column=0)
    cpane.grid_remove()
    rna_browse.grid_remove()
    rna_file_selection.grid_remove()

    # Pad:
    Label(root, text = "        ").grid(row=10,column=0)

    # Select Method:
    Label(root, text="2. Choose tumor deconvolution method:", font="Times 16 bold").grid(row=11, column=0, sticky = W)
    m = StringVar()
    m.set('CIBERSORT (Newman et al., 2015)')
    methods = ["CIBERSORT (Newman et al., 2015)", "DeconRNASeq (Gong and Szustakowski, 2013)", "ssGSEA (Barbie et al., 2009) (Senbabaoglu, 2016)", "SingScore (Foroutan, 2018)"]
    OptionMenu(root, m, *methods, command=show_method_options).grid(row=12, sticky = W)


    # Pad:
    Label(root, text="      ").grid(row=16,column=0)

    # Upload signature matrix (or choose lm22)
    Label(root, text="3. Additional Options:", font="Times 16 bold").grid(row=17,column=0, sticky=W)
    s, u, d = StringVar(), StringVar(), StringVar()
    s.set('LM22')
    u.set('ssGSEA paper gene set')
    d.set('None')
    sig_label = Label(root, text="Signature Matrix:")
    sig_label.grid(row=18, column=0, sticky=W)
    # sig_label.grid_remove()
    sig_mat_optionslist = ["LM22", "Upload your own..."]
    sig_dropdown = OptionMenu(root, s, *sig_mat_optionslist, command=show_sig_options)
    sig_dropdown.grid(row=19, sticky = W)
    # sig_dropdown.grid_remove()
    cpane2 = CollapsiblePane(root, 'Collapse Instructions', 'Instructions for uploading your own signature matrix data file:')
    cpane2.grid(row=20,sticky = W)
    cpane2.grid_remove()
    Label(cpane2.frame,text="- Data file extension must be either .csv, .xlsx, or .txt").grid(row=21,sticky = W)
    Label(cpane2.frame,text="- First column must contain a list of genes, referenced by their Hugo Symbol").grid(row=22,sticky = W)
    Label(cpane2.frame,text="- Each remaining column should be the signature expression value in the column's cell type, for each gene").grid(row=23,sticky = W)
    Label(cpane2.frame,text="- First row of file should be a header, containing column names").grid(row=24,sticky = W)
    sig_browse = Button(root, text="Browse...", command=get_custom_signature)
    sig_browse.grid(row=25, column=1)
    sig_browse.grid_remove()
    sig_selection = Label(root, text = "")
    sig_selection.grid(row=26,column=0)
    # sig_selection.grid_remove()

    # Normalization :
    normalization = StringVar()
    normalization.set("Z-Scaling")
    norm_label = Label(root, text="Data Normalization Method:")
    norm_label.grid(row=27, column=0, sticky=W)
    # norm_label.grid_remove()
    norm_dropdown = OptionMenu(root, normalization, "Z-Scaling", "MinMax Normalization", "No Normalization")
    norm_dropdown.grid(row=28, column=0, sticky=W)
    # norm_dropdown.grid_remove()

    # Upload gene sets (or choose the default)
    upgene_label = Label(root, text="Up-Regulated Gene Set:")
    upgene_label.grid(row=18, column=0, sticky=W)
    upgene_label.grid_remove()
    up_gene_optionslist = ["ssGSEA paper gene set", "Upload your own..."]
    upgene_dropdown = OptionMenu(root, u, *up_gene_optionslist, command=show_upgenes_options)
    upgene_dropdown.grid(row=19, sticky = W)
    upgene_dropdown.grid_remove()
    cpane3 = CollapsiblePane(root, 'Collapse Instructions', 'Instructions for uploading your own gene sets data file:')
    cpane3.grid(row=20,sticky = W)
    cpane3.grid_remove()
    Label(cpane3.frame,text="- Data file extension must be either .csv, .xlsx, or .txt").grid(row=21,sticky = W)
    Label(cpane3.frame,text="- First row should be a list of cell types (column names)").grid(row=22,sticky = W)
    # File should have columns named for each cell type,
    Label(cpane3.frame,text="- Rows should contain Hugo Symbols for genes to be considered as up-regulated for the given cell type").grid(row=23,sticky = W)
    Label(cpane3.frame,text="- If not all cell types have the same number of up-regulated genes, excess rows in that column should be coded as 'NaN'").grid(row=24,sticky = W)
    upgenes_browse = Button(root, text="Browse...", command=get_custom_up_geneset)
    upgenes_browse.grid(row=25, column=1)
    upgenes_browse.grid_remove()
    upgenes_selection = Label(root, text = "")
    upgenes_selection.grid(row=26,column=0)
    upgenes_selection.grid_remove()

    downgene_label = Label(root, text="Down-Regulated Gene Set (optional, for bidirectional SingScore):")
    downgene_label.grid(row=27, column=0, sticky=W)
    downgene_label.grid_remove()
    down_gene_optionslist = ["None", "Upload your own..."]
    downgene_dropdown = OptionMenu(root, d, *down_gene_optionslist, command=show_downgenes_options)
    downgene_dropdown.grid(row=28, sticky = W)
    downgene_dropdown.grid_remove()
    downgenes_browse = Button(root, text="Browse...", command=get_custom_down_geneset)
    downgenes_browse.grid(row=29, column=1)
    downgenes_browse.grid_remove()
    downgenes_selection = Label(root, text = "")
    downgenes_selection.grid(row=30,column=0)
    downgenes_selection.grid_remove()

    # Alpha (occupies same space as normalization & downgenes)
    alpha_label = Label(root, text="Alpha value:")
    alpha_label.grid(row=27, column=0, sticky=W)
    alpha_label.grid_remove()
    alpha = Entry(root)
    alpha.grid(row=28, column=0, sticky=W)
    alpha.insert(0, "1.0")
    alpha.grid_remove()

    # Pad:
    Label(root, text="      ").grid(row=29, column=0)

    # Save results:
    Label(root, text="4. Choose location to save results:", font="Times 16 bold").grid(row=31,column=0,sticky = W)
    Button(root, text="Choose folder...", command=get_save_location).grid(row=32,column=0)
    print_save_loc = tk.Label(root, text=" ")
    print_save_loc.grid(row=33, column=0)

    # Pad:
    Label(root, text="    ").grid(row=34, column=0)

    Label(root, text="5. Run TumorDecon:", font="Times 16 bold").grid(row=35,column=0,sticky=W)
    Button(root, text="Run!", command=run_tumor_decon).grid(row=36,column=0)

    status = tk.Label(root, text=" ")
    status.grid(row=37, column=0)

    root.mainloop()
