#!/opt/Python/2.7.3/bin/python
import sys
from collections import defaultdict
from numpy import *
import re
import os
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC

def usage():
    test="name"
    message='''
Python LOD_Genome.py --input ../input/MPR.cross.uniq.QTL.mr.table

Convert position from chromosome into single genome position.

    '''
    print message

def fasta(fastafile):
    fasta_seq = defaultdict(str)
    p = re.compile(r'(\d+)')
    for record in SeqIO.parse(fastafile,"fasta"):
        m = p.search(record.id)
        rank = m.groups(0)[0] if m else 0
        fasta_seq[rank] = str(record.seq)
    return fasta_seq

'''
Convert LOD
""      "chr"   "pos"   "Heading.Days"  "Plant.Height..cm..in.Field"    "Biomass"       "Number.of.Tillers"     "Single.Plant..Grain.Yield...g."        
"0100222046"    "1"     0       0.00726849104738392     1.75434974017193        0.706883828535311       0.745160416472634       1.50828253471261        
"0100500860"    "1"     0.38023 0.00249508279722122     2.48456985459131        0.606550483403985       0.662605759329817       1.29921364448124        
"0100609676"    "1"     2.584996        7.66947829689002e-08    2.94093242870076        0.514032120297038       0.648206340956266       1.18183467128868
'''
def convert_LOD(infile, outfile, chr_start):
    data = defaultdict(lambda : str)
    last_pos = defaultdict(lambda : int)
    ofile = open(outfile, 'w')
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            line = line.replace('"', '')
            if line[0:1].isdigit(): 
                unit = re.split(r'\t',line)
                pos = int(unit[0][2:]) + chr_start[unit[1]]
                if last_pos.has_key(unit[1]):
                    pos = (pos - last_pos[unit[1]])/2 + last_pos[unit[1]] 
                    last_pos[unit[1]] = pos
                else:
                    last_pos[unit[1]] = chr_start[unit[1]]
                    pos = (pos - last_pos[unit[1]])/2 + last_pos[unit[1]]
                    last_pos[unit[1]] = pos
                #print int(unit[0][2:]), pos
                unit[0] = str(pos)
                line = '\t'.join(unit)
                print >> ofile, line
            else:
                unit = re.split(r'\t',line)
                unit[0] = 'Position'
                line = '\t'.join(unit)
                print >> ofile, line
    ofile.close()
    return data

def get_start(sequences):
    chr_start = defaultdict(lambda : int)
    last_start= 0
    for seq_id in sorted(sequences.keys(), key=int):
        chr_start[seq_id] = last_start
        #print seq_id, last_start
        length = len(sequences[seq_id])
        last_start += length 
        midpoint = int((last_start-chr_start[seq_id])/2 + chr_start[seq_id])
        chrn = 'Chr%s' %(seq_id)
        print '%s\t%s\t%s' %(chrn, str(midpoint), str(last_start))
    return chr_start 

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

    try:
        len(args.output) > 0
    except:
        
        args.output = args.input + '.new'

    ref = '/rhome/cjinfeng/BigData/00.RD/seqlib/MSU_r7.fa'
    fasta_seq = fasta(ref)
    chr_start = get_start(fasta_seq)
    if not os.path.isfile(args.output):
        convert_LOD(args.input, args.output, chr_start)
 
if __name__ == '__main__':
    main()

