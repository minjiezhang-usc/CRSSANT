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
Created on Tue Nov 1 2021

@author: Minjie Zhang (minjiez@usc.edu)

Step 1: merger gap1sam and transsam
Step 2: prepare bedgraphs for DG assembly
"""

import os
import argparse
import time

def merger_gap_trans():
    gap1align = open(gap1sam,"r")
    transalign = open(transsam,"r")
    outalign = open(outprefix+'_crssant.sam',"w")
    for line in gap1align: 
        if line[0]=="@": outalign.write(line)
        else:
            if "SA:Z:" not in line.split()[-1]: outalign.write(line)
            if "SA:Z:" in line.split()[-1]: outalign.write('\t'.join(line.split()[:-2])+'\n')
    for line in transalign:
        if line[0]!="@": outalign.write(line)
    gap1align.close()
    transalign.close()
    outalign.close()
    
def get_bedgraph():
    sam_2_bam_cmd = 'samtools view -bS -o {} {}'
    sam_2_bam_cmd = sam_2_bam_cmd.format(outprefix+'_crssant.bam', outprefix+'_crssant.sam')
    
    bam_sort_cmd = 'samtools sort -o {} {}'
    bam_sort_cmd = bam_sort_cmd.format(outprefix+'_crssant_sorted.bam', outprefix+'_crssant.bam')
    
    get_bedgraph_plus_cmd = 'bedtools genomecov -bg -split -strand + -ibam {} > {}'
    get_bedgraph_plus_cmd = get_bedgraph_plus_cmd.format(outprefix+'_crssant_sorted.bam', outprefix+'_crssant_plus.bedgraph')
    
    get_bedgraph_minus_cmd = 'bedtools genomecov -bg -split -strand - -ibam {} > {}'
    get_bedgraph_minus_cmd = get_bedgraph_minus_cmd.format(outprefix+'_crssant_sorted.bam', outprefix+'_crssant_minus.bedgraph')
    
    os.system(sam_2_bam_cmd)
    os.system(bam_sort_cmd)
    os.system(get_bedgraph_plus_cmd)
    os.system(get_bedgraph_minus_cmd)
    os.system('rm %s %s' %(outprefix+'_crssant.bam', outprefix+'_crssant_sorted.bam'))
    

if __name__ == '__main__':
    # Inputs from user
    parser = argparse.ArgumentParser(description='Prepare files for CRSSANT')
    parser.add_argument('gap1sam', help='Name of gap1 filtered sam')
    parser.add_argument('transsam', help='Name of trans sam')
    parser.add_argument('outprefix', help='Name of mergered sam')

    args = parser.parse_args()
    gap1sam = args.gap1sam
    transsam = args.transsam
    outprefix = args.outprefix
    
    merger_gap_trans()
    get_bedgraph()