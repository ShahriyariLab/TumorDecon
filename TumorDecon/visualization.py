# visualization.py - authored by Sumeyye Su & Trang Le

import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from TumorDecon.data_utils import *

# combine the cells that belong to same family:
def norm_and_combine(df):
    # Check if freqs already sum to 1, and scale them first if not:
    h=np.ones(df.shape[0])
    k=df.sum(axis=1).values
    if not np.array_equal(h,k):
        df=sum_upto_1(df)

    df = combine_celltypes(df)
    return df

# scale the frequences so they sum up to 1
def sum_upto_1(df):
    df[df<0]=0
    s=df.sum(axis=1)
    df=df.divide(s,axis=0)*1
    return df


def cell_frequency_boxplot(sample_cell_freq, xsize=12, ysize=7):
    """
     Input:
        - 'sample_cell_freq': A dataframe that include cell frequency of samples
            Rows are samples id, columns are cell names
    Output:
        - cells frequency box plot in descending order
    """
    print(sample_cell_freq)
    new_cell_freq = norm_and_combine(sample_cell_freq)
    b=new_cell_freq.median(axis = 0)
    b=list(zip(b.index,b))
    b = sorted(b, key=lambda x: x[-1],reverse=True)
    sorted_cells=[x[0] for x in b]
    new_cell_freq=new_cell_freq[sorted_cells]
    sns.set(rc={'figure.figsize':(xsize,ysize)})
    sns.set(style="white")
    palette={'Macrophages':'violet','CD8 T cells':'orange','CD4 T cells':'goldenrod','Monocytes':'lightsalmon','NK cells':'olivedrab','Mast cells':'red','B cells':'darkcyan','T cells gamma delta':'dodgerblue','Dendritic cells':'gray','Plasma cells':'seagreen','Neutrophils':'navy', 'Eosinophils':'purple'}
    sns.boxplot(x="Patient_ID", y="value",order=sorted_cells, data=pd.melt(new_cell_freq),palette=palette)
    plt.xlabel('')
    plt.xticks(rotation=90)
    plt.ylabel('Frequencey')
    plt.subplots_adjust(top=0.95,bottom=0.2)
    plt.show()
    return


def cell_frequency_barchart(sample_cell_freq, title=" ", xsize=15, ysize=7):
    """
     Input:
        - 'sample_cell_freq': A dataframe that include cell frequency of samples
            Rows are samples id, columns are cell names
    Output:
        - a barchart plot that shows all cell frequency of all samples
    """
    new_cell_freq=norm_and_combine(sample_cell_freq)
    b=new_cell_freq.median(axis = 0)
    b=list(zip(b.index,b))
    b = sorted(b, key=lambda x: x[-1],reverse=True)
    sorted_cells=[x[0] for x in b]
    new_cell_freq=new_cell_freq[sorted_cells]
    palette={'Macrophages':'violet','CD8 T cells':'orange','CD4 T cells':'goldenrod','Monocytes':'lightsalmon','NK cells':'olivedrab','Mast cells':'red','B cells':'darkcyan','T cells gamma delta':'dodgerblue','Dendritic cells':'gray','Plasma cells':'seagreen','Neutrophils':'navy', 'Eosinophils':'purple'}
    sns.set(rc={'figure.figsize':(xsize,ysize)})
    sns.set_style("white")

    color=[]
    for gene in sorted_cells:
        a=palette[gene]
        color=color+[a]

    ax=new_cell_freq.plot.barh(stacked=True,color=color,width=1)
    ax.legend(loc='right', bbox_to_anchor=(1.25,0.5),
          fancybox=True, shadow=False, ncol=1)

    plt.yticks([])
    plt.margins(0, 0)
    plt.xlabel('Frequency')
    plt.ylabel('Patients')
    plt.title(title)
    plt.subplots_adjust(top=0.95,left=0.05, right=0.8) # don't cut off legend
    plt.show()
    return


def hierarchical_clustering(sample_cell_freq, xsize=20, ysize=5):
    """
     Input:
        - 'sample_cell_freq': A dataframe that include cell frequency of samples
            Rows are samples id, columns are cell names
    Output:
        - hierarchical clustered heatmap from cell frequency of samples
    """
    # Change figure size:
    plt.rcParams["figure.figsize"] = (xsize,ysize)

    new_cell_freq=norm_and_combine(sample_cell_freq)

    b=new_cell_freq.median(axis = 0)
    b=list(zip(b.index,b))
    b = sorted(b, key=lambda x: x[-1],reverse=True)
    sorted_cells=[x[0] for x in b]
    new_cell_freq=new_cell_freq[sorted_cells]
    sns.clustermap(new_cell_freq,cmap='coolwarm',row_cluster=False)
    plt.subplots_adjust(bottom=0.2)
    plt.show()
    return

def pair_plot(sample_cell_freq, xsize=5, ysize=5):
    """
     Input:
        - 'sample_cell_freq': A dataframe that include cell frequency of samples
           Rows are samples id, columns are cell names
    Output:
        - pairplot from cell frequency of samples
    """
    # Change figure size:
    plt.rcParams["figure.figsize"] = (xsize,ysize)
    new_cell_freq=norm_and_combine(sample_cell_freq)
    b=new_cell_freq.median(axis = 0)
    b=list(zip(b.index,b))
    b = sorted(b, key=lambda x: x[-1],reverse=True)
    sorted_cells=[x[0] for x in b]
    new_cell_freq=new_cell_freq[sorted_cells]

    p = sns.pairplot(new_cell_freq)
    p.fig.set_size_inches(xsize,ysize)
    plt.subplots_adjust(bottom=0.1, top=1.0, left=0.05, right=0.95)
    plt.show()
    return


def stack_barchart(methods, results, true_freqs, cell_types, colors, fig_size, fig_name):
    sample_labels = np.arange(1,len(true_freqs)+1)

    fig, axs = plt.subplots(1, len(methods)+1, sharey=True, figsize=fig_size)
    # Remove vertical space between axes
    fig.subplots_adjust(wspace=0.1)

    for i in range(len(methods)):
        results[methods[i]][cell_types].plot.barh(ax = axs[i], stacked=True, color=colors, width=0.9)
        axs[i].get_legend().remove()
        axs[i].set_yticklabels(sample_labels)
        axs[i].set_title(methods[i])
        axs[i].set_ylabel('Sample')
        if methods[i] in ['ssGSEA', 'singscore']:
            axs[i].set_xlabel('Score')
        else:
            axs[i].set_xlabel('Frequency')

    true_freqs[cell_types].plot.barh(ax = axs[-1], stacked=True, color=colors, width=0.9)
    axs[-1].get_legend().remove()
    axs[-1].set_yticklabels(sample_labels)
    axs[-1].set_xlabel('Frequency')
    axs[-1].set_title('Ground_truth')

    fig_leg = plt.figure(figsize=(3, 3))
    ax_leg = fig_leg.add_subplot(111)
    # add the legend from the previous axes
    ax_leg.legend(*axs[-1].get_legend_handles_labels(), loc='center')
    # hide the axes frame and the x/y labels
    ax_leg.axis('off')
    # fig_leg.savefig('Figures/'+fig_name+'_stack_barchart_legend.eps')

    # fig.savefig('Figures/'+fig_name+'_stack_barchart.eps')
