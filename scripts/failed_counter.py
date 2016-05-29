#!usr/bin/env python
# -*- coding: utf-8 -*-

#Load modules
from collections import namedtuple, defaultdict, OrderedDict
import csv
import glob
import re
import sys

import matplotlib.pyplot as plt
import numpy as np
from operator import itemgetter
import pybedtools

#Load BED file & define directory to parse over
INTERVALS_BED = "/media/partition/Haloplex/00100-1407755742_Regions.bed"
DIRECTORY = '/media/partition/Haloplex/Haloplex_Test_2_Mid_February/Velona'

#Define coverage thresholds
THRESHOLD_1 = 1
THRESHOLD_2 = 50
THRESHOLD_3 = 100
THRESHOLD_4 = 250
THRESHOLD_5 = 500
THRESHOLD_6 = 1000

#Prepare BED file
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
print('Preparing BED file:... DONE')

#Search BAM files in the input directory
bam_file_list = [f for f in glob.iglob(DIRECTORY+"/*.bam")]
print('Searching for BAM files:... DONE')

#Count how many amplicons do not reach a certain threshold
collected = defaultdict(list)
print('Check for failed amplicons')
for file in bam_file_list:

    almnt = pybedtools.BedTool(file)

    coverage_result = almnt.coverage(intervals_list).sort()

    for interval in coverage_result:
        if interval[3] not in collected.keys():
            collected[interval[3].encode('ascii','ignore')] = {
            'Failed 1' : 0,
            'Failed 1-50' : 0,
            'Failed 50-100' : 0,
            'Failed 100-250' : 0,
            'Failed 250-500' : 0,
            'Failed 500-1000' : 0
            }

        if float(interval[4]) <= THRESHOLD_1 :
            collected[interval[3]]['Failed 1'] += 1
        elif THRESHOLD_1 < float(interval[4]) <= THRESHOLD_2:
            collected[interval[3]]['Failed 1-50'] += 1
        elif THRESHOLD_2 < float(interval[4]) <= THRESHOLD_3:
            collected[interval[3]]['Failed 50-100'] += 1
        elif THRESHOLD_3 < float(interval[4]) <= THRESHOLD_4:
            collected[interval[3]]['Failed 100-250'] += 1
        elif THRESHOLD_4 < float(interval[4]) <= THRESHOLD_5:
            collected[interval[3]]['Failed 250-500'] += 1
        elif THRESHOLD_5 < float(interval[4]) <= THRESHOLD_6:
            collected[interval[3]]['Failed 500-1000'] += 1
print('Check for failed amplicons:... DONE')

print('Write data into a CSV')
#Write result in a CSV file
with open('counter.csv',"w") as f_out:
    writer = csv.writer(f_out)
    writer.writerow(['Amplicon', 'Failed 1', 'Failed 1-50', 'Failed 50-100', 'Failed 100-250', 'Failed 250-500', 'Failed 500-1000'])
    for key, values in collected.items():
        data = [key, values['Failed 1'], values['Failed 1-50'], values['Failed 50-100'], values['Failed 100-250'], values['Failed 250-500'], values['Failed 500-1000']]
        writer.writerow(data)
