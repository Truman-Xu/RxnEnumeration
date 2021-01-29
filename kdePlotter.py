# KDE Plotting Script for Reaction Enumeration Project
# Plots as seen in https://doi.org/10.1038/s41586-020-2142-y

# BSD 3-Clause License
# Copyright (c) 2020-2021, Cernak Lab at the University of Michigan

import pickle
import numpy as np
import pandas as pd

import matplotlib
import matplotlib.pyplot as plt
matplotlib.use('Agg') # use it for writing file only, it disables plt.show()
matplotlib.matplotlib_fname()

from rdkit import Chem
from rdkit.Chem import AllChem, Draw
from rdkit.Chem import rdMolDescriptors, Descriptors

import seaborn as sns
from collections import namedtuple

AllMolProps = namedtuple('AllMolProps', ['MolWt','LogP','HBA','HBD','PSA',
                                         'RotB','FSP3','FC','QED','Arom'])
def getAllProps(mol):
    
    return AllMolProps(
        round(Chem.Descriptors.MolWt(mol),2),
        round(Chem.rdMolDescriptors.CalcCrippenDescriptors(mol)[0],2),
        Chem.rdMolDescriptors.CalcNumHBA(mol),
        Chem.rdMolDescriptors.CalcNumHBD(mol),
        Chem.rdMolDescriptors.CalcTPSA(mol),
        Chem.rdMolDescriptors.CalcNumRotatableBonds(mol),
        Chem.rdMolDescriptors.CalcFractionCSP3(mol),
        Chem.rdmolops.GetFormalCharge(mol),
        round(Chem.QED.qed(mol),2),
        Chem.rdMolDescriptors.CalcNumAromaticRings(mol)
    )

def smiListToDF(smi_list):
    '''
    convert SMILES to Mol
    Calculates the selected Mol Props and make a dataframe
    Requires getAllProps() function
    '''
    mols = [Chem.MolFromSmiles(smi) for smi in smi_list]
    df = pd.DataFrame([getAllProps(m) for m in mols])
    df['SMILES'] = smi_list
    
    return df

def smiFileToDF(filename):
    '''
    Opens a text file of SMILES and convert them to DataFrame
    Uses smiListToDF(smi_list)
    '''
    with open(filename,'r') as f:
        smi_list = [l.rstrip() for l in f.readlines()]
    
    return smiListToDF(smi_list)


def ridgeplot_preview(df, save = None):
    '''
    Fast preview all the KDE plots with default sns theme and original data structure
    '''
    # Selecte the plottable datatypes
    df = df.select_dtypes(include=['float64','int64'])
    # Determine how many rows needed
    n_rows = len(df.columns)//3 + (len(df.columns) % 3 > 0)
    
    fig, axes = plt.subplots(n_rows, 3, figsize=(10,10))

    i = 0
    for row in axes:
        for ax in row:
            if i < len(df.columns):
                sns.kdeplot(ax = ax, data = df, x = df.columns[i], fill = True)
                ax.set_title(df.columns[i], fontweight='bold', fontsize=13, fontname='Times New Roman')
                d_min = df[df.columns[i]].min()
                d_max = df[df.columns[i]].max()
                d_d = d_max - d_min
                if d_min != d_max:
                    ax.set(xlim = [d_min-d_d,d_max+d_d])

                ax.set(xlabel=None)
                ax.set(ylabel=None)
                ax.set(yticks=[])
                ticks = df[df.columns[i]].unique()
                if len(ticks) > 3:
                    ticks = [d_min, d_max]
                ax.set(xticks = ticks)
                ax.tick_params(axis='both', which='major', labelsize=10)
                i+=1
            else: 
                ax.axis('off')
                
    plt.tight_layout()
    print('Plotting Finished.')
    
    if save:
        plt.savefig(save)
        print('Plot save at "{}"'.format(save))
    
    return fig


def ridgeplot_bo(df, save = None):
    '''
    Ridgeplot theme and data structure designed by Bo
    '''
    # Transform to a stacked df
    df2 = df.drop(columns = ['Arom','SMILES'])
    for c in df2.columns:
        df2[c] = df2[c].transform(lambda x: np.interp(x, [df2[c].min(),df2[c].max()], (-0, +1)))
    df2 = df2.stack().droplevel(0)
    df2 = df2.reset_index().rename(columns ={'index':'attribute', 0:'value'})
    
    #### Code below copied from the original script ridgeplot.py ####
    
    sns.set(rc={'figure.figsize':(10,10)})
    sns.set(style="white", rc={"axes.facecolor": (0, 0, 0, 0)},font_scale=1.5)
    pal = sns.cubehelix_palette(10, rot=-.25, light=.7)

    # set up the facet grid
    g = sns.FacetGrid(df2, col="attribute", hue="attribute", aspect=1, height=2.25, palette=pal, col_wrap=3, sharex=False)

    # plot the kde plots
    g.map(sns.kdeplot,"value",clip_on=True,shade=True,alpha=1,lw=1.5,bw_method=0.2)
    g.map(sns.kdeplot, "value", clip_on=False, color="w", lw=2, bw_method=0.2)
    g.map(plt.axhline, y=0, lw=2, clip_on=False)

    g.set(yticks=[])
    g.set(xticks=[0.0,1.0])
    g.despine(top=False, right=False)

    for i, ax in enumerate(g.axes.flat): # set every-other axis for testing purposes

        ax.set_title(df.columns[i], fontweight='bold', fontsize=11, fontname='Times New Roman')
        ax.tick_params(axis='both', which='major', labelsize=11)
        ax.set_xticklabels([round(df[df.columns[i]].min(),1),round(df[df.columns[i]].max(),1)])
        ax.set_xlim(-.5,1.5)

    g.fig.tight_layout()
    
    print('Plotting Finished.')
    if save:
        g.savefig(save)
        print('Plot save at "{}"'.format(save))
        
    
    return g.fig

if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(prog='ridgeplot')
    parser.add_argument("in_file", metavar = 'I', type=str, 
                        help="SMILES file for making the plot. The file should be one SMILES \
                        per line with no header or other information")
    parser.add_argument("-o","--out_file", type=str,
                        help="file name for the plot, default to 'ridgeplot.png'", 
                        default = "./ridgeplot.png")
    parser.add_argument("-s","--plot_style", type=str,
                        help="Choose the ridgeplot style between 'bo' and 'preview', \
                        default to 'bo'", 
                        default = "bo")
    a = parser.parse_args()
    
    df = smiFileToDF(a.in_file)
    
    if a.plot_style == 'bo':
        ridgeplot_bo(df, save = a.out_file)
    elif a.plot_style == 'preview':
        ridgeplot_preview(df, save = a.out_file)
    else: 
        raise 'Invalid style option "{}", please choose from "bo" or "simple"'.format(a.plot_style)