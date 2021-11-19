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
gapmcluster.py
In CRSSANT, the clustering is performed on two arms, we can also cluster them
based on multiple segments, here using 3M as an example.
1. Get all DGs, build an interval tree
2. Build all TGs based on overlapping alternative DGs
3. Find the biggest overlaps of each alignment to all TGs


1. Analysis of RN7SK data
cd /Users/lu/Documents/lulab/projects/psoralen/HEK1/HEK1_hg38
Here we used manually curated DG info for the TG assembly because the automatic
assembly of DGs does not work well for the alternatives in the 5' end. The DGs
are first made on the RN7SK coordinates and then lifted to the hg38 coordinates
for easy processing. Together, ...
"""


import sys, itertools, subprocess, datetime, re
# for use on changrila
import intervaltree
from datetime import datetime

if len(sys.argv) < 3:
    print("Usage: python gapmcluster.py bedpe gapm.sam outprefix")
    print("Assume the DGs are all on the same chromsome and strand")
    sys.exit()

bedpe = open(sys.argv[1], 'r')
gapmsam = open(sys.argv[2], 'r')
gapmout = open(sys.argv[3]+'_tg.sam', 'w')
over_threshold = 0.5 #overlap of overlapped arm in alternative structures
alt_threshold = 0.2 #overlap of the other arm in alternative structures


########1. function definitions
################################################################################
def overlap(al, ar, bl, br):
    #two intervals a [al, ar] and b [bl, br]
    #returns the ratio of overlap
    if ar < bl or al > br: return 0
    else:
        alength, blength = ar - al, br - bl
        coords = sorted([al, ar, bl, br])
        overlaplength = float(coords[2] - coords[1])
        overlap = min(overlaplength/alength, overlaplength/blength)
        return overlap
def alternativecheck(intvl1, intvl2, over_threshold, alt_threshold):
    #returns whether the two intervals overlap
    #use the overlap() function defined above
    if (overlap(intvl1[0],intvl1[1],intvl2[0],intvl2[1]) > over_threshold) and \
       (overlap(intvl1[2],intvl1[3],intvl2[2],intvl2[3]) < alt_threshold) or \
       (overlap(intvl1[0],intvl1[1],intvl2[0],intvl2[1]) < alt_threshold) and \
       (overlap(intvl1[2],intvl1[3],intvl2[2],intvl2[3]) > over_threshold) or \
       (overlap(intvl1[0],intvl1[1],intvl2[2],intvl2[3]) > over_threshold) or \
       (overlap(intvl1[2],intvl1[3],intvl2[0],intvl2[1]) > over_threshold):
        return 1
    else: return 0
    
def CIGARdiv(CIGAR):
    #considers all CIGAR operations [MINDSPH=X]
    #divide a cigar string to N-separated segments, used in functions below
    #return these results:
    #gaps, a list of N strings, e.g. ['100N', '200N']
    #segs, a list of segments, each with all possible operations except N
    #gaplens, a list of N lengths, e.g. [100, 200]
    #Glen, genomic length of the entire alignment (consumed ref) [MDN=X]
    #Mlens, segment lengths of matches only [M=X]
    #Rlens, segment lengths of consumed Reference [MD=X]
    #Qlens, segment lengths of consumed Query [MIS=X]
    #example: 
    #Glen,gaps,gaplens,segs,Mlens,Qlens,Rlens = CIGARdiv(CIGAR)
    #or replace unwanted output with '_' 
    gaps=re.findall('\d+N', CIGAR)
    gaplens=[int(gap[:-1]) for gap in gaps] #gap lengths 
    segs=[i.rstrip('0123456789') for i in CIGAR.split('N')]
    Mlens=[sum([int(i[:-1]) for i in re.findall('\d+[M=X]',s)]) for s in segs]
    Rlens=[sum([int(i[:-1]) for i in re.findall('\d+[MD=X]',s)]) for s in segs]
    Qlens=[sum([int(i[:-1]) for i in re.findall('\d+[MIS=X]',s)]) for s in segs]
    Glen=sum(gaplens+Rlens)
    return Glen, gaps, gaplens, segs, Mlens, Qlens, Rlens


########2. construct a DG tree for each chrom, the bedpe file may not be sorted
################################################################################
print(str(datetime.now())[:-7]+'\t'+"Started TGs assembly"+'\n')
#This analysis does not deal with structures on different chroms/strands
dgtrees = {} #each dgtree is named after the chrom name
intervals = {} #positions for each chrom is stored in a list
for line in bedpe:
    bed = line.strip('\n').split()
    chrom = bed[0]
    begin1, end1 = int(bed[1]), int(bed[2])
    begin2, end2 = int(bed[4]), int(bed[5]) 
    if chrom not in dgtrees:
        dgtrees[chrom] = intervaltree.IntervalTree()
        intervals[chrom] = []
    dgtrees[chrom].addi(begin1, end1, line)
    dgtrees[chrom].addi(begin2, end2, line)
    for i in range(begin1, end1):
        intervals[chrom].append(i)
    for i in range(begin2, end2):
        intervals[chrom].append(i)    
bedpe.close()


#########3. convert dgtrees to armpairs
################################################################################
treesize, intervallength, armpairscount = 0, 0,0
chroms = dgtrees.keys()
armpairs = {}
#example armpair: (Interval(454, 475, 'bedrecord'), Interval(454, 476, 'bedrecord'))
for chrom in chroms:
    treesize += len(dgtrees[chrom])
    intervals[chrom] = sorted(list(set(intervals[chrom])))
    intervallength += len(intervals[chrom])
    #make all dgclusters from this chromosome
    dgclusters = set()
    for i in intervals[chrom]:
        dgcluster = tuple(dgtrees[chrom][i])
        if dgcluster and len(dgcluster) > 1: dgclusters.add(dgcluster)
    #make all armpairs from this chromosome
    armpairs[chrom] = set()
    for dgcluster in dgclusters:
        for dg1, dg2 in list(itertools.combinations(dgcluster, 2)):
            armpairs[chrom].add(tuple(sorted((dg1, dg2))))
    armpairscount += len(armpairs[chrom])
#print("\n", str(datetime.datetime.today()))
print("                       Number of all DGs:", treesize/2)
print("                 Length of all intervals:", intervallength)
print("                      Number of armpairs:", armpairscount)

#########4. convert alternatives to tgs (triple-segment groups)
################################################################################
tgtrees = {} # [start1,end1, [start2,end2,start3,end3,otherinfo]]
altpairs = []
altpairsbed = [] #bed12_1 + structure_1 + bed12_2 + structure_2 + maxoverlap.
for chrom in chroms:
    for armpair in armpairs[chrom]:
        dg1, dg2 = armpair # a tuple of two intervals
        #print(dg1, dg2); sys.exit()
        dg1info = dg1.data.strip('\n').split("\t")
        dg2info = dg2.data.strip('\n').split("\t")
        intvl1  = [int(i) for i in dg1info[1:3]+dg1info[4:6]]
        intvl2  = [int(i) for i in dg2info[1:3]+dg2info[4:6]]
        ###test alternative structures
        #maxoverlap,structure1,structure2=bpcheck(chrom,intvl1,intvl2,fastadict)
        if alternativecheck(intvl1,intvl2,over_threshold,alt_threshold):
            ll = overlap(intvl1[0],intvl1[1],intvl2[0],intvl2[1])
            rl = overlap(intvl1[2],intvl1[3],intvl2[0],intvl2[1])
            lr = overlap(intvl1[0],intvl1[1],intvl2[2],intvl2[3])
            rr = overlap(intvl1[2],intvl1[3],intvl2[2],intvl2[3])
            intvls = [] #[start1, end1, start2, end2, start3, end3]
            if ll == max(ll,rl,lr,rr):
                intvls = [min(intvl1[0],intvl2[0]),max(intvl1[1],intvl2[1]),
                          min(intvl1[2],intvl2[2]),min(intvl1[3],intvl2[3]),
                          max(intvl1[2],intvl2[2]),max(intvl1[3],intvl2[3])]
            elif rr == max(ll,rl,lr,rr):
                intvls = [min(intvl1[0],intvl2[0]),min(intvl1[1],intvl2[1]),
                          max(intvl1[0],intvl2[0]),max(intvl1[1],intvl2[1]),
                          min(intvl1[2],intvl2[2]),max(intvl1[3],intvl2[3])]
            elif rl == max(ll,rl,lr,rr):
                intvls = [intvl1[0],intvl1[1],min(intvl1[2],intvl2[0]),
                          max(intvl1[3],intvl2[1]),intvl2[2],intvl2[3]]
            elif lr == max(ll,rl,lr,rr):
                intvls = [intvl2[0],intvl2[1],min(intvl1[0],intvl2[2]),
                          max(intvl1[1],intvl2[3]),intvl1[2],intvl1[3]]
            if chrom not in tgtrees: tgtrees[chrom]=intervaltree.IntervalTree()
            tginfo = intvls[2:6]+dg1info[6:8]+dg2info[6:8]
            tgtrees[chrom].addi(intvls[0], intvls[1], tginfo)
           

#########5. assign alignments to TGs. 
################################################################################
#right now only check 3-segment alignments. For 
for line in gapmsam:
    #QNAME 0 US47 1 255 12M6N28M4N39M4S * 0 0 SEQ QUAL
    if line[0] == "@": gapmout.write(line); continue
    align = line.split()
    RNAME, POS, CIGAR = align[2], int(align[3]), align[5]
    #if RNAME !="chr6": continue
    STRAND = '-' if '{0:012b}'.format(int(align[1]))[-5] == '1' else '+'
    _,gaps,gaplens,segs,_,_,Rlens = CIGARdiv(CIGAR)
    intvls = [] 
    for i in range(len(segs)): #n gaps and n+1 matches
        intvls.append([POS+sum(gaplens[:i])+sum(Rlens[:i]), \
                       POS+sum(gaplens[:i])+sum(Rlens[:i+1])])
    tgs = '' #get the TGs (triple-segment groups)
    if len(intvls) == 3 and RNAME in tgtrees:
        tgs = list(tgtrees[RNAME][intvls[0][0]:intvls[0][1]])
    #which TG to associate with?
    overlaplist = []
    for tg in tgs:
        overlap1 = overlap(intvls[0][0],intvls[0][1],tg.begin,tg.end)
        overlap2 = overlap(intvls[1][0],intvls[1][1],tg.data[0],tg.data[1])
        overlap3 = overlap(intvls[2][0],intvls[2][1],tg.data[2],tg.data[3])
        nreads = int(tg.data[5])*int(tg.data[7])
        overlaplist.append([overlap1*overlap2*overlap3, nreads, tg])
    overlaplist.sort() #sort by the overlap and numbers of reads in the 2 DGs
    #if overlaplist: print('\n', overlaplist[-1])
    #output a TG list in bed12 format. Cannot use the bedpe format due to 3 segs
    if overlaplist and overlaplist[-1][0]>0:
        tgid = '_'.join(overlaplist[-1][2].data[4:8])
        gapmout.write(line.strip('\n')+ '\tTG:Z:' + tgid + '\n')

gapmsam.close()
gapmout.close()

print("############################################################################"+'\n')







