#!/usr/bin/python
# -*- coding: utf-8 -*-

#import modules
import pybedtools
from collections import namedtuple, defaultdict, OrderedDict
import re
import matplotlib.pyplot as plt
import glob
import numpy as np

IntervalColumns = namedtuple('bed', ['chr', 'start', 'end', 'gene'])
intervals_list = []

INTERVALS_BED = "/media/partition/00100-1407755742_Regions.bed"
#INTERVALS_BED = '/media/partition/tst15/TST_15-A-manifest.bed'

if INTERVALS_BED:
    with open(INTERVALS_BED, "r") as fin:
        for line in fin.readlines():
            line = line.rstrip('\n')
            line = re.split(r'\t+', line.rstrip('\t'))
            if len(line) == 4:
                line[1] = int(line[1])
                line[2] = int(line[2])
                bed_line = IntervalColumns(*(line[0], int(line[1]), int(line[2]), line[3]))
                intervals_list.append(bed_line)
            if len(line) == 12:
                line[1] = int(line[1])
                line[2] = int(line[2])
                bed_line = IntervalColumns(*(line[0], int(line[1]), int(line[2]), line[3]))
                intervals_list.append(bed_line)
else:
    print("ERROR: Provide an interval list (bed format)")

directory = '/media/partition/collected/hpx_csc_surecall'
#directory = '/media/partition/tst15/MixA'

fig,ax1 = plt.subplots()
plt.hold = True

bam_file_list = [f for f in glob.iglob(directory+"/*.bam")]

boxes = []
collected = defaultdict(list)

for file in bam_file_list:

    print file

    coverage_list = []
    almnt = pybedtools.BedTool(file)

    coverage_result = almnt.coverage(intervals_list).sort()

    basename = re.sub('.bam$','',file)
    basename = re.sub(directory,'',basename)

    collected[basename] = {
    'Coverages' : [],
    'Median' : 0
    }

    cov_counter = []

    for interval in coverage_result:
        cov_counter.append(float(interval[4]))

    cov_mean = np.mean(cov_counter)
    cov_median = np.median(cov_counter)
    '''
    for interval in coverage_result:
        collected[basename]['Coverages'].append(float(interval[4]))
    '''
    for interval in coverage_result:
        collected[basename]['Coverages'].append(float(interval[4]) / cov_mean)

for key, value in collected.items():
    median = np.median(collected[key]['Coverages'])
    collected[key]['Median'] += median
    collected_sorted = OrderedDict(sorted(collected.items(), key=lambda t: t[1]['Median']))


bla = pd.DataFrame(collected)

print bla

g = sns.barplot('Sample','Ratio',data=bla,color='b')
g, labels = plt.xticks()
plt.setp(labels, rotation=45)
plt.show()

'''

for key, value in collected_sorted.items():
    boxes.append(collected_sorted[key]['Coverages'])

plt.boxplot(boxes,sym='')
xtickNames = plt.setp(ax1, xticklabels = collected.keys())
plt.setp(xtickNames, rotation=90, fontsize=7)
plt.ylabel('Coverage (x)')
plt.title('Normalized Coverage Distribution per Sample')

plt.show()

'''
