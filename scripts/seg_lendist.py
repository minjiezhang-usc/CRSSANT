#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# RNA crosslinking, proximity ligation and high throughput sequencing produces non-continuous reads that indicate base pairing and higher order interactions,\
# either in RNA secondary structures or intermolecular complexes. CRSSANT (pronounced 'croissant') is a computational pipeline for analyzing non-continuous/gapped \
# reads from a variety of methods that employ the crosslink-ligation principle, including PARIS, LIGR, SPLASH, COMRADES, hiCLIP, etc. CRSSANT optimizes short-read \
# mapping, automates alignment processing, and clusters gap1 and trans alignments into duplex groups (DG) and non-overlapping groups (NG). More complex arrangments \
# are assembled into higher level structures. In particular gapm alignments with 2 gaps or 3 segments are assembled into tri-segment groups (TGs). Overlapping alignments \
# are used to discover homotypic interactions (RNA homodimers).

# Briefly, the CRSSANT pipeline operates as follows. First, sequencing reads that have been processed to remove adapters are mapped references with STAR and a new \
# set of optimized options. Second, alignments are filtered, rearranged and classified into different types (gaptypes.py and gapfilter.py). Third, we use network \
# analysis methods to cluster non-continuous alignments into DGs and calculate the confidence for each DG. The DGs are used as the foundation for the assembly of TGs.

# CRSSANT is written in Python and available as source code that you can download and run directly on your own machine (no compiling needed). An earlier version of \
# the DG assembly method is available here: (https://github.com/ihwang/CRSSANT). For more about the CRSSANT pipeline, please see the bioRxiv preprint by Zhang et al. 2021.

"""
Created on Tue Nov 01 2021
This script calculates the size distribution of each segment in STAR alignments. 
"""

import matplotlib.pyplot as plt
import sys, re, matplotlib, matplotlib.pyplot as plt

if len(sys.argv)< 3:
    print("Usage: ")
    print("python  seg_lendist.py  inputfile  outprefix")
    print("inputfile:   gap1 sam file")
    print("outprefix:   output prefix")
    sys.exit()

inputfile = open(sys.argv[1], 'r')
outprefix = sys.argv[2]


##part1: save the size list as space delimited numbers in one line in a file. 
sizelist = []
for line in inputfile:
    if line[0] == "@": continue
    CIGAR = line.split()[5]
    segs=[i.rstrip('0123456789') for i in CIGAR.split('N')]
    Mlens=[sum([int(i[:-1])for i in re.findall('\d+[M=X]',s)])for s in segs]
    sizelist += Mlens
#save the size list as space delimited numbers in one line in a file. 
inputfile.close()


##part2: plot the segment size distribution
fig, ax = plt.subplots(figsize=(3.0,2))
plt.subplots_adjust(bottom=0.32)
plt.subplots_adjust(top=0.75)
plt.subplots_adjust(left=0.25)
plt.subplots_adjust(right=0.95)

n, bins, patches = plt.hist(sizelist, [x+0.5 for x in range(0,200)], \
                                histtype='step', cumulative=True, density=1) 
plt.setp(patches, 'facecolor', 'b', 'alpha', 0.75)
plt.xlim(0, 100)
plt.ylim(0, 1.1)
max_yticks = 5
max_xticks = 10
yloc = plt.MaxNLocator(max_yticks)
xloc = plt.MaxNLocator(max_xticks)
ax.yaxis.set_major_locator(yloc)
ax.set_xlabel("segment length (nt)") #, fontsize=15
ax.set_ylabel("cumulative frequency") #, fontsize=15
plt.savefig(outprefix+'_SegLen_distribution.pdf')
#plt.show()

