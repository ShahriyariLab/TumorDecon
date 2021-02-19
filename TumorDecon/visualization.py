# visualization.py - originally authored by Sumeyye Su & Trang Le

import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from TumorDecon.data_utils import *
from collections.abc import Mapping

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


def cell_frequency_boxplot(sample_cell_freq, title="", title_fontsize=10, save_as=None, axes_style=None, font_scale=1, xylabel_fontsize=10, rcParams={'figure.figsize':(15,7)}):
    """
     Input:
        - 'sample_cell_freq': A dataframe that include cell frequency of samples
            Rows are samples id, columns are cell names
        - 'title': string.
        - 'title_fontsize': int. Font-size for title.
        - 'save_as': string. Filename to save plot as (full or relative path). File format is
            inferred from filename's extension ('.png', '.pdf', '.eps', etc)
                If None, plot is not saved.
        - 'axes_style': a dictionary of parameters to pass to sns.set_style, to customize plot
            Valid dictionary keys can be found by running "sns.axes_style()" with no arguments
        - 'font_scale': float. Scaling for plot elements such as tick labels
        - 'xylabel_fontsize': int. Font-size for axis labels.
        - 'rcParams': a dictionary of parameters to pass as the rc argument to
            sns.set function
                Valid dictionary keys can be found at https://matplotlib.org/stable/tutorials/introductory/customizing.html
    Output:
        - cells frequency box plot in descending order
    """
    new_cell_freq = norm_and_combine(sample_cell_freq)
    b=new_cell_freq.median(axis = 0)
    b=list(zip(b.index,b))
    b = sorted(b, key=lambda x: x[-1],reverse=True)
    sorted_cells=[x[0] for x in b]
    new_cell_freq=new_cell_freq[sorted_cells]
    # Customize the plots with user-defined parameters
    sns.set(rc=rcParams, font_scale=font_scale)
    if axes_style is not None:
        if not isinstance(axes_style, Mapping):
            raise TypeError('axes_style must be a dict')
        sns.set_style("white", axes_style)
    else:
        sns.set(style="white", font_scale=font_scale)

    palette={'Macrophages':'violet','CD8 T cells':'orange','CD4 T cells':'goldenrod','Monocytes':'lightsalmon','NK cells':'olivedrab','Mast cells':'red','B cells':'darkcyan','T cells gamma delta':'dodgerblue','Dendritic cells':'gray','Plasma cells':'seagreen','Neutrophils':'navy', 'Eosinophils':'purple'}
    try:
        sns.boxplot(x="Patient_ID", y="value",order=sorted_cells, data=pd.melt(new_cell_freq),palette=palette)
    except ValueError:
        sns.boxplot(x="variable", y="value",order=sorted_cells, data=pd.melt(new_cell_freq),palette=palette)
    plt.xlabel('')
    plt.xticks(rotation=90)
    plt.ylabel('Frequency',fontsize=xylabel_fontsize)
    plt.title(title, fontsize=title_fontsize)
    plt.subplots_adjust(top=0.95,bottom=0.2)
    if save_as is not None:
        plt.savefig(save_as)
    plt.show()
    return


def cell_frequency_barchart(sample_cell_freq, title=" ", title_fontsize=10, save_as=None, axes_style=None, font_scale=1, xylabel_fontsize=10, rcParams={'figure.figsize':(15,7)}):
    """
     Input:
        - 'sample_cell_freq': A dataframe that include cell frequency of samples
            Rows are samples id, columns are cell names
        - 'title': string.
        - 'title_fontsize': int. Font-size for title.
        - 'save_as': string. Filename to save plot as (full or relative path). File format is
            inferred from filename's extension ('.png', '.pdf', '.eps', etc)
                If None, plot is not saved.
        - 'axes_style': a dictionary of parameters to pass to sns.set_style, to customize plot
            Valid dictionary keys can be found by running "sns.axes_style()" with no arguments
        - 'font_scale': float. Scaling for plot elements such as tick labels
        - 'xylabel_fontsize': int. Font-size for axis labels.
        - 'rcParams': a dictionary of parameters to pass as the rc argument to
            sns.set function
                Valid dictionary keys can be found at https://matplotlib.org/stable/tutorials/introductory/customizing.html

    Output:
        - a barchart plot that shows all cell frequency of all samples
    """
    new_cell_freq=norm_and_combine(sample_cell_freq)
    b=new_cell_freq.median(axis = 0)
    b=list(zip(b.index,b))
    b = sorted(b, key=lambda x: x[-1],reverse=True)
    sorted_cells=[x[0] for x in b]
    new_cell_freq=new_cell_freq[sorted_cells]
    # Customize the plots with user-defined parameters
    sns.set(rc=rcParams, font_scale=font_scale)
    if axes_style is not None:
        if not isinstance(axes_style, Mapping):
            raise TypeError('axes_style must be a dict')
        else:
            sns.set_style("white", axes_style)
    else:
        sns.set(style="white", font_scale=font_scale)
    palette={'Macrophages':'violet','CD8 T cells':'orange','CD4 T cells':'goldenrod','Monocytes':'lightsalmon','NK cells':'olivedrab','Mast cells':'red','B cells':'darkcyan','T cells gamma delta':'dodgerblue','Dendritic cells':'gray','Plasma cells':'seagreen','Neutrophils':'navy', 'Eosinophils':'purple'}

    color=[]
    for gene in sorted_cells:
        a=palette[gene]
        color=color+[a]

    ax=new_cell_freq.plot.barh(stacked=True,color=color,width=1)
    ax.legend(loc='right', bbox_to_anchor=(1.25,0.5),
          fancybox=True, shadow=False, ncol=1)

    plt.yticks([])
    plt.margins(0, 0)
    plt.xlabel('Frequency', fontsize=xylabel_fontsize)
    plt.ylabel('Patients', fontsize=xylabel_fontsize)
    plt.title(title, fontsize=title_fontsize)
    plt.subplots_adjust(top=0.95,left=0.05, right=0.8) # don't cut off legend
    if save_as is not None:
        plt.savefig(save_as)
    plt.show()
    return


def hierarchical_clustering(sample_cell_freq, title="", title_fontsize=10, font_scale=1, save_as=None, figsize=(20,5)):
    """
    Input:
        - 'sample_cell_freq': A dataframe that include cell frequency of samples
            Rows are samples id, columns are cell names
        - 'title': string.
        - 'title_fontsize': int. Font-size for title.
        - 'font_scale': float. Scaling for plot elements such as tick labels
        - 'save_as': string. Filename to save plot as (full or relative path). File format is
            inferred from filename's extension ('.png', '.pdf', '.eps', etc)
                If None, plot is not saved.
        - 'figsize': tuple. (x,y) size of the figure
    Output:
        - hierarchical clustered heatmap from cell frequency of samples
    """

    new_cell_freq=norm_and_combine(sample_cell_freq)

    b=new_cell_freq.median(axis = 0)
    b=list(zip(b.index,b))
    b = sorted(b, key=lambda x: x[-1],reverse=True)
    sorted_cells=[x[0] for x in b]
    new_cell_freq=new_cell_freq[sorted_cells]
    sns.set(font_scale=font_scale)
    sns.clustermap(new_cell_freq,cmap='coolwarm',row_cluster=False, figsize=figsize)
    plt.subplots_adjust(bottom=0.2)
    plt.title(title, fontsize=title_fontsize)
    if save_as is not None:
        plt.savefig(save_as)
    plt.show()
    return

def pair_plot(sample_cell_freq, title="", title_fontsize=10, font_scale=1, save_as=None, figsize=(20,20)):
    """
     Input:
        - 'sample_cell_freq': A dataframe that include cell frequency of samples
           Rows are samples id, columns are cell names
       - 'title': string.
       - 'title_fontsize': int. Font-size for title.
       - 'font_scale': float. Scaling for plot elements such as tick labels
       - 'save_as': string. Filename to save plot as (full or relative path). File format is
            inferred from filename's extension ('.png', '.pdf', '.eps', etc)
                If None, plot is not saved.
       - 'figsize': tuple. (x,y) size of the figure
     Output:
        - pairplot from cell frequency of samples
    """
    new_cell_freq=norm_and_combine(sample_cell_freq)
    b=new_cell_freq.median(axis = 0)
    b=list(zip(b.index,b))
    b = sorted(b, key=lambda x: x[-1],reverse=True)
    sorted_cells=[x[0] for x in b]
    new_cell_freq=new_cell_freq[sorted_cells]

    sns.set(font_scale=font_scale)
    p = sns.pairplot(new_cell_freq)
    p.fig.set_size_inches(*figsize)
    plt.subplots_adjust(bottom=0.1, top=1.0, left=0.05, right=0.95)
    plt.title(title, fontsize=title_fontsize)
    if save_as is not None:
        plt.savefig(save_as)
    plt.show()
    return


def stack_barchart(methods, results, true_freqs, cell_types, colors, figsize=(10,10), save_as=None):
    """
    Plot for comparing predicted frequencies of each method to true data
     Input:
        - 'methods': 1D string array.
                All methods with results to plot.
        - 'results': 3D float array (num_methods x num_cell_types x num_samples)
                All TumorDecon results from each method.
        - 'true_freqs': 2D float array (num_cell_types x num_samples)
                Experimental cell frequencies to compare results to.
        - 'cell_types': 1D string array.
                All cell_types listed in the results
        - 'colors': 1D string array of length num_cell_types.
                Enties must be valid matplotlib colors string (to associate with each given cell_type)
        - 'save_as': string. Filename to save plot as (full or relative path). File format is
            inferred from filename's extension ('.png', '.pdf', '.eps', etc)
                If None, plot is not saved.
        - 'figsize': tuple. (x,y) size of the figure
     Output:
        - pairplot from cell frequency of samples
    """
    sample_labels = np.arange(1,len(true_freqs)+1)

    fig, axs = plt.subplots(1, len(methods)+1, sharey=True, figsize=figsize)
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

    if save_as is not None:
        plt.savefig(save_as)
