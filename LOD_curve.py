#!/opt/Python/2.7.3/bin/python
import sys
from collections import defaultdict
import re
import os
import argparse
from Bio import SeqIO
import numpy as np
from scipy.cluster.hierarchy import dendrogram, linkage
import scipy.spatial.distance as distance
import matplotlib as mp
mp.use('pdf')
import matplotlib.pyplot as plt
from matplotlib import gridspec
import pandas as pd


def usage():
    test="name"
    message='''
python CircosConf.py --input circos.config --output pipe.conf

    '''
    print message

def set_ticks_XY(ax, ypos, ylim, chrs):
    fig = plt.gcf()
    #fig.set_size_inches(8, 11)

    # turn off the frame
    ax.set_frame_on(True)
    ax.yaxis.set_ticks_position('right')
    ax.xaxis.set_ticks(chrs[2], minor=True)

    # put the major ticks at the middle of each cell
    ax.set_yticks(ypos, minor=False)
    ax.set_xticks(chrs[1], minor=False)

    # want a more natural, table-like display
    #ax.invert_yaxis()
    ax.xaxis.tick_top()

    # Set the labels
    ax.set_xticklabels(chrs[0], minor=False)
    ax.set_yticklabels(ylim, minor=False)

    # rotate the
    #plt.xticks(rotation=0)

    ax.grid(False)

    # Turn off all the ticks
    ax = plt.gca()


    for t in ax.xaxis.get_major_ticks():
    #    t.tick1On = False
        t.tick2On = False
    for t in ax.yaxis.get_major_ticks():
        t.tick1On = False
        t.tick2On = False
 

    return ax


def set_ticks_XY_Right(ax, ypos, ylim):
    fig = plt.gcf()
    #fig.set_size_inches(8, 11)

    # turn off the frame
    ax.set_frame_on(True)
    ax.yaxis.set_ticks_position('right')

    # put the major ticks at the middle of each cell
    ax.set_yticks(ypos, minor=False)
    #ax.set_xticks(np.arange(ylen) + 0.5, minor=False)

    # want a more natural, table-like display
    #ax.invert_yaxis()
    ax.xaxis.tick_top()

    # Set the labels
    ax.set_xticklabels([], minor=False)
    ax.set_yticklabels(ylim, minor=False)

    # rotate the
    #plt.xticks(rotation=0)

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



def set_ticks_XY_empty(ax):
    fig = plt.gcf()
    #fig.set_size_inches(8, 11)

    # turn off the frame
    ax.set_frame_on(False)

    # put the major ticks at the middle of each cell
    #ax.set_yticks(np.arange(xlen) + 0.5, minor=False)
    #ax.set_xticks(np.arange(ylen) + 0.5, minor=False)

    # want a more natural, table-like display
    #ax.invert_yaxis()
    ax.xaxis.tick_top()

    # Set the labels
    ax.set_xticklabels([], minor=False)
    #ax.set_yticklabels([], minor=False)

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



def LOD_curve(lod ,mat, labs, chrs, lodcut):
    # clustering
    dist_mat = mat
    linkage_matrix = linkage(dist_mat, "single")

    # Create a figure.
    figsize=(20,10)
    fig = plt.figure(figsize=figsize)
    plt.subplots_adjust(bottom=0.2)
  
    # split plot into 4 rows X 2 cols. gs[1, 0] is 2th row and 1th column, which is left bottom region.
    # gs[0:, 0]  is all rows of 1th column
    ncurve = 10.00
    gs = gridspec.GridSpec(int(ncurve), 2, wspace=0.2, hspace=0.15, width_ratios=[0.15,1], height_ratios=[1/ncurve]*int(ncurve))

    #dendrogram
    ax0 = fig.add_subplot(gs[0:, 0])
    ddata = dendrogram(linkage_matrix,
                   color_threshold=1,
                   orientation='right',
                   labels=labs)
    
    # Assignment of colors to labels: 'a' is red, 'b' is green, etc.
    #label_colors = {'a': 'r', 'b': 'g', 'c': 'b', 'd': 'm'}
    
    # set axis
    ax0 = set_ticks_XY_empty(ax0) 


    # Draw curve for each panel
    #for i in range(0, int(ncurve)):
    #    ax1 = fig.add_subplot(gs[i, 1])
    #    ax1.yaxis.set_ticks_position('right')
    #    ddata = dendrogram(linkage_matrix,
    #               color_threshold=1,
    #               labels=labs)
        #ax1 = set_ticks_XY_Right(ax1)
    count = 0
    for i in ddata['leaves']:
        ax1 = fig.add_subplot(gs[int(ncurve)-count-1, 1])
        ax1.yaxis.set_ticks_position('right')
        trait = labs[i]
        #print count, i, trait
        #plt.xlim=(10,20)
        x = lod['Position']
        y = lod[trait]
        xmin = 0 
        xmax = int (max(x))
        ymax = int (max(y) * 1.1) if int (max(y) * 1.1) > 5 else 5
        ymin = 0
        ax1.set_xlim([xmin, xmax])
        ax1.set_ylim([ymin, ymax])
        ax1.plot(x,y)
        print count, i, trait, float(lodcut[trait])
        ax1.axhline(y=float(lodcut[trait]), color='Orange')
        count += 1 
        if count <= int(ncurve)-1:
            ax1 = set_ticks_XY_Right(ax1, [ymin, ymax*0.9], [ymin, ymax])
        else:
            ax1 = set_ticks_XY(ax1, [ymin, ymax*0.9], [ymin, ymax], chrs)
        #count += 1
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

    lodfile = '../input/MPR.cross.uniq.QTL.mr.table.new'
    lodcutfile = '../input/MPR.cross.uniq.QTL.mr.table.LOD_threshold.new'
    chrfile = '../input/MSU7.Chr.midpoint'
    pdist = pd.read_table(args.input, index_col=0, header=None)
    midpoint = pd.read_table(chrfile, header=None)
    lod   = pd.read_table(lodfile)
    lodcut= pd.read_table(lodcutfile)   
 
    mat = np.array([[100,  80,  50],
                [1,  0.8, 0.5],
                [30, 80,  80],
                [0.5,  0.8, 1]])

    pairwise_dists = distance.squareform(distance.pdist(mat))
    LOD_curve(lod, pdist, pdist.index, midpoint, lodcut)

if __name__ == '__main__':
    main()

