#!usr/bin/env python
# -*- coding: utf-8 -*-

'''
This script plots the normalized coverage distribution per amplicon in a
given directory. The average coverage per target region is normalized by
dividing this value by the overall average coverage of each target sequences.
'''

#Load modules
import pybedtools
import re
import sys
import glob
from collections import namedtuple, defaultdict, OrderedDict
import numpy as np
import matplotlib.pyplot as plt
from operator import itemgetter
import seaborn as sns
import pandas as pd

#INTERVALS_BED = '/media/partition/TST_15-A-manifest.bed'
#INTERVALS_BED = '/media/partition/TST_15-B-manifest.bed'
INTERVALS_BED = '/media/partition/00100-1407755742_Regions.bed'

COVERAGE_THRESHOLD = 1000

#directory = '/media/partition/tst15/MixA'
#directory = '/media/partition/tst15/MixB'
directory = '/media/partition/collected/test'

fig,ax1 = plt.subplots()
plt.hold = True
boxes = []

bam_file_list = [f for f in glob.iglob(directory+"/*.bam")]

collected = defaultdict(list)

for file in bam_file_list:

    almnt = pybedtools.BedTool(file)

    IntervalColumns = namedtuple('bed', ['chr', 'start', 'end', 'gene'])
    intervals_list = []

    if INTERVALS_BED:
        with open(INTERVALS_BED, "r") as fin:
            for line in fin.readlines():
                line = line.rstrip('\n')
                line = re.split(r'\t+', line.rstrip('\t'))

                # For Agilent Haloplex BED files
                if len(line) == 4:
                    line[1] = int(line[1])
                    line[2] = int(line[2])
                    bed_line = IntervalColumns(*(line[0], int(line[1]), int(line[2]), line[3]))
                    intervals_list.append(bed_line)

                # For Illumina Trusight Tumor 15 BED files
                if len(line) == 12:
                    line[1] = int(line[1])
                    line[2] = int(line[2])
                    bed_line = IntervalColumns(*(line[0], int(line[1]), int(line[2]), line[3]))
                    intervals_list.append(bed_line)
    else:
        print("ERROR: Provide an interval list (bed format)")

    coverage_result = almnt.coverage(intervals_list).sort()

    print('\nFile: %s \n' %file)

    cov_counter = []

    for interval in coverage_result:
        cov_counter.append(float(interval[4]))

    cov_mean = np.mean(cov_counter)

    for interval in coverage_result:
        collected[interval[3].encode('ascii','ignore')].append(float(interval[4]) / cov_mean)


for key, value in collected.items():
    median = np.median(collected[key])
    collected[key]['Median'] += median
    collected_sorted = OrderedDict(sorted(collected.items(), key=lambda t: t[1]['Median']))

map(collected_sorted.pop, ['Median'])

results_data_frame = pd.DataFrame(collected)

print results_data_frame

g = sns.boxplot(data=results_data_frame,color='b')

plt.show()
