#!/opt/Python/2.7.3/bin/python
import sys
from collections import defaultdict
import re
import os
import argparse
from Bio import SeqIO
import numpy as np
from scipy.cluster.hierarchy import dendrogram, linkage
import matplotlib as mp
mp.use('pdf')
import matplotlib.pyplot as plt
from matplotlib import gridspec


def usage():
    test="name"
    message='''
python CircosConf.py --input circos.config --output pipe.conf

    '''
    print message


def set_ticks_XY_Right(ax):
    fig = plt.gcf()
    #fig.set_size_inches(8, 11)

    # turn off the frame
    ax.set_frame_on(True)
    ax.yaxis.set_ticks_position('right')

    # put the major ticks at the middle of each cell
    #ax.set_yticks(np.arange(xlen) + 0.5, minor=False)
    #ax.set_xticks(np.arange(ylen) + 0.5, minor=False)

    # want a more natural, table-like display
    #ax.invert_yaxis()
    ax.xaxis.tick_top()

    # Set the labels
    #ax.set_xticklabels(xlab, minor=False, fontsize=8)
    #ax.set_yticklabels(ylab, minor=False)

    # rotate the
    plt.xticks(rotation=0)

    #ax.grid(False)

    # Turn off all the ticks
    ax = plt.gca()

    for t in ax.xaxis.get_major_ticks():
        t.tick1On = False
        t.tick2On = False
    for t in ax.yaxis.get_major_ticks():
        t.tick1On = False
        t.tick2On = False

    return ax



def set_ticks_XY_empty(ax):
    fig = plt.gcf()
    #fig.set_size_inches(8, 11)

    # turn off the frame
    ax.set_frame_on(False)

    # put the major ticks at the middle of each cell
    #ax.set_yticks(np.arange(xlen) + 0.5, minor=False)
    #ax.set_xticks(np.arange(ylen) + 0.5, minor=False)

    # want a more natural, table-like display
    ax.invert_yaxis()
    ax.xaxis.tick_top()

    # Set the labels
    ax.set_xticklabels([], minor=False)
    ax.set_yticklabels([], minor=False)

    # rotate the
    plt.xticks(rotation=90)

    ax.grid(False)

    # Turn off all the ticks
    ax = plt.gca()

    for t in ax.xaxis.get_major_ticks():
        t.tick1On = False
        t.tick2On = False
    for t in ax.yaxis.get_major_ticks():
        t.tick1On = False
        t.tick2On = False

    return ax



def LOD_curve(mat):
    # clustering
    dist_mat = mat
    linkage_matrix = linkage(dist_mat, "single")

    # Create a figure.
    figsize=(8,11)
    fig = plt.figure(figsize=figsize)
    plt.subplots_adjust(bottom=0.2)
  
    # split plot into 4 rows X 2 cols. gs[1, 0] is 2th row and 1th column, which is left bottom region.
    # gs[0:, 0]  is all rows of 1th column
    gs = gridspec.GridSpec(4, 2, wspace=0.1, hspace=0.1, width_ratios=[0.25,1], height_ratios=[0.25,0.25,0.25,0.25])

    #dendrogram
    ax0 = fig.add_subplot(gs[0:, 0])
    ddata = dendrogram(linkage_matrix,
                   color_threshold=1,
                   orientation='right',
                   labels=["a", "b", "c", "d"])
    
    # Assignment of colors to labels: 'a' is red, 'b' is green, etc.
    label_colors = {'a': 'r', 'b': 'g', 'c': 'b', 'd': 'm'}
    
    # set axis
    ax0 = set_ticks_XY_empty(ax0) 


    # Draw curve for each panel
    for i in range(0, 4):
        ax1 = fig.add_subplot(gs[i, 1])
        ax1.yaxis.set_ticks_position('right')
        ddata = dendrogram(linkage_matrix,
                   color_threshold=1,
                   labels=["a", "b", "c", "d"])
        #ax1 = set_ticks_XY_Right(ax1)

    #ax = plt.gca()
    #xlbls = ax.get_xmajorticklabels()
    #for lbl in xlbls:
    #    lbl.set_color(label_colors[lbl.get_text()])

    fig.savefig('%s.pdf' %('QTL_LOD_curve'), bbox_inches='tight')

def readtable(infile):
    data = defaultdict(str)
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2: 
                unit = re.split(r'\t',line)
                if not data.has_key(unit[0]):
                    data[unit[0]] = unit[1]
    return data


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input')
    parser.add_argument('-o', '--output')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()
    try:
        len(args.input) > 0
    except:
        usage()
        sys.exit(2)
    
    mat = np.array([[1.0,  0.5,  0.0],
                [0.5,  1.0, -0.5],
                [1.0, -0.5,  0.5],
                [0.0,  0.5, -0.5]])

    LOD_curve(mat)

if __name__ == '__main__':
    main()

