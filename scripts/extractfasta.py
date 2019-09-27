#! /usr/bin/env python

"""
Extract full or partial sequences from a fasta file
Target specifications can be a bed file or entered on the command line (single entry)
"""

import sys
from string import maketrans

def readbed(filename):
    bed = {}
    with open(filename, 'r') as fin:
        for line in fin:
            line_l = line.strip().split()
            if len(line_l) > 2:
                chrom, start, stop = line_l[0:3]
                name = '{}_{}_{}'.format(chrom, start, stop)
                score = 0
                strand = '+'
            if len(line_l) > 3:
                name = line_l[3]
            if len(line_l) > 4:
                score = int(line_l[4])
            if len(line_l) > 5:
                strand = line_l[5]
            if chrom not in bed:
                bed[chrom] = []
            bed[chrom].append([int(start), int(stop), name, strand])
    return bed

try:
    bed = readbed(sys.argv[1])
except IOError:
    feature, start, stop = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
    name = '{}_{}_{}'.format(feature, start, stop)
    strand = '+'
    bed = {sys.argv[1]: [[start, stop, name, strand]]}
sequence = ''
header = ''
head = ''
trans = maketrans('ACGTacgt', 'TGCAtgca')
for line in sys.stdin:
    if line.startswith('>'):
        if header != '':
            if header in bed or head in bed:
                if header in bed:
                    useheader = header
                else:
                    useheader = head
                for start, stop, name, strand in bed[useheader]:
                    if strand != '-':
                        sys.stdout.write('>{}\n{}\n'.format(name, sequence[start:stop]))
                    else:
                        out = sequence[-(start+1):-(stop+1):-1].translate(trans)
                        sys.stdout.write('>{}\n{}\n'.format(name, out))
        header = line[1:].strip()
        head = header.split()[0]
        sequence = ''
    else:
        sequence += line.strip()
else:
    if header != '':
        if header in bed or head in bed:
            if header in bed:
                useheader = header
            else:
                useheader = head
            for start, stop, name, strand in bed[useheader]:
                if strand != '-':
                    sys.stdout.write('>{}\n{}\n'.format(name, sequence[start:stop]))
                else:
                    out = sequence[-(start+1):-(stop+1):-1].translate(trans)
                    sys.stdout.write('>{}\n{}\n'.format(name, out))

