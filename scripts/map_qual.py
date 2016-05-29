#!/usr/bin/python
# -*- coding: utf-8 -*-

#import modules
import os
import pysam
import matplotlib.pyplot as plt
import numpy
import csv
import sys
from collections import namedtuple, defaultdict
import re

#input_bamfile =  "/media/partition/TST15_Test_1_Early_February/Base_Space/BRAF_15018040/Libraries/15018040_S1.bam"
#INTERVALS = "/media/partition/TST15_Test_1_Early_February/TST_15-A-manifest.bed"

input_bamfile = "/media/partition/Haloplex_Test_1_Late_January/Velona/15061857_S2.bam"
INTERVALS = "/media/partition/Haloplex_Test_1_Late_January/00100-1407755742_Regions.bed"

bamfile = pysam.AlignmentFile(input_bamfile, "rb")

IntervalColumns = namedtuple('bed', ['chr', 'start', 'end', 'name'])
intervals_list = []


if INTERVALS:
    # NOTE BED format is 0 based (open interval)
    with open(INTERVALS, "r") as fin:
        for line in fin.readlines():
            line = line.rstrip('\n')
            line = re.split(r'\t+', line.rstrip('\t'))
            if len(line) == 4:
                line[1] = int(line[1])
                line[2] = int(line[2])
                bed_line = IntervalColumns(*(line[0], int(line[1]), int(line[2]), line[3]))
                intervals_list.append(bed_line)
else:
    print("ERROR: Provide either an interval list (bed format) or a vcf file")
print(intervals_list)

quals = []
for interval in intervals_list:
    start = int(interval.start)
    end = int(interval.end)
    length = end - start
    pos_key = "%s:%s:%s" % (interval.chr, interval.start, interval.end)
    for aligned_pos in range(start, start + length):
        for read in bamfile.fetch(reference=interval.chr, start=start, end=end):
            qual = read.mapping_quality
            quals.append(qual)

#print numpy.mean(quals)
plt.hist(quals)
plt.title("Mapping Quality Histograms")
plt.xlabel("Mapping Quality")
plt.ylabel("Count")
plt.show()
