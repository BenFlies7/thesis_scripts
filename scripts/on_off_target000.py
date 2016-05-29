#!/usr/bin/python
# -*- coding: utf-8 -*-

'''
This script calculates the ratio of reads aligning to genes
defined in a BED file comparing to the total aligned reads.
It then plots the result as boxplots
'''

#import modules
from collections import namedtuple, defaultdict
import pysam
import re
import sys
import matplotlib.pyplot as plt
import numpy as np
import glob
import seaborn as sns
import pandas as pd

#Load files
REFERENCE = "/media/partition/hg19_broadinstitute/ucsc.hg19.fasta"
INTERVALS_BED = "/media/partition/00100-1407755742_Regions.bed"
#INTERVALS_BED = "/media/partition/haloplex_csc_genes.bed"
#directory = "/media/partition/hpx_csc_velona/"
directory = "/media/partition/test"



#Prepare BED file
IntervalColumns = namedtuple('bed', ['chr', 'start', 'end'])
intervals_list = []

if INTERVALS_BED:
    with open(INTERVALS_BED, "r") as fin:
        for line in fin.readlines():
            line = line.rstrip('\n')
            line = re.split(r'\t+', line.rstrip('\t'))
            prc = ()

            #For Agilent Haloplex BED file
            if len(line) == 4:
                line[1] = int(line[1])
                line[2] = int(line[2])
                bed_line = IntervalColumns(*(line[0], int(line[1]), int(line[2])))
                intervals_list.append(bed_line)

            #For Illumina Trusight Tumor 15 BED files
            if len(line) == 12:
                line[1] = int(line[1])
                line[2] = int(line[2])
                bed_line = IntervalColumns(*(line[0], int(line[1]), int(line[2])))
                intervals_list.append(bed_line)
else:
    print("ERROR: Provide an interval list (bed format)")

fig = plt.figure()
ax = fig.add_subplot(111)

bam_file_list = [f for f in glob.iglob(directory+"/*.bam")]

boxes = []
#boxes2 = http://stackoverflow.com/questions/26540035/rotate-label-text-in-seaborn-factorplothttp://stackoverflow.com/questions/26540035/rotate-label-text-in-seaborn-factorplot[]
collected = defaultdict(list)
ratio_list = []

file_list = []

collected = {
'Sample' : [],
'Ratio' : []
}

for file in bam_file_list:

    print('currently in file %s') %file

    basename = re.sub('.bam$','',file)
    basename = re.sub(directory,'',basename)

    file_list.append(basename)

    total_mapped_counter = 0
    on_target_chr_counter = 0

    samfile = pysam.AlignmentFile(file, "rb")
    chr_lengths = samfile.lengths

    chr_list = ['chrM','chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY']

    for i, chr in enumerate(chr_list):

        on_target_chr_list = []

        total_mapped_chr = samfile.count(chr_list[i], 1, chr_lengths[i])
        total_mapped_counter += total_mapped_chr
        #total_mapped_chr_list.append(total_mapped_chr)

        for interval in intervals_list:
            if chr == interval.chr:
                on_target_chr = samfile.count(interval.chr, interval.start, interval.end)
                on_target_chr_counter += on_target_chr
                #on_target_chr_list.append(float(on_target_chr / total_mapped_chr))

    if total_mapped_chr > 0:
        ratio = on_target_chr_counter / float(total_mapped_counter)
        collected['Sample'].append(basename)
        collected['Ratio'].append(ratio)

bla = pd.DataFrame(collected)

print bla

'''
g = sns.barplot('Sample','Ratio',data=bla,color='b')
g, labels = plt.xticks()
plt.setp(labels, rotation=45)
plt.show()
'''

g = sns.boxplot(y='Ratio',data=bla)
g, labels = plt.xticks()
plt.setp(labels)
plt.show()
