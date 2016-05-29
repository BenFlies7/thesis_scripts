#!usr/bin/env python
# -*- coding: utf-8 -*-

#TODO: fix

#Load modules
import pybedtools
import re
import sys
from collections import namedtuple, defaultdict, OrderedDict
import numpy as np
import matplotlib.pyplot as plt
from operator import itemgetter

COVERAGE_THRESHOLD = 1000

#BAM_FILE = "/media/partition/Haloplex/Haloplex_Test_2_Mid_February/Velona/15016513_S5.bam"
#BAM_FILE = "/media/partition/Haloplex/Haloplex_Test_2_Mid_February/Velona/15028422_S6.bam"
BAM_FILE = "/media/partition/Haloplex/Haloplex_Test_2_Mid_February/Velona/ECD_S9.bam"

INTERVALS_BED = "/media/partition/Haloplex/00100-1407755742_Regions.bed"

almnt = pybedtools.BedTool(BAM_FILE)

IntervalColumns = namedtuple('bed', ['chr', 'start', 'end', 'gene'])
intervals_list = []

if sys.argv[1] == 'all':
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
elif sys.argv[1] == 'interest':
    if INTERVALS_BED:
        with open(INTERVALS_BED, "r") as fin:
            for line in fin.readlines():
                line = line.rstrip('\n')
                line = re.split(r'\t+', line.rstrip('\t'))
                if len(line) == 4:
                    if (re.sub('_\d$','',line[3]) == 'EGFR') \
                    or (re.sub('_\d$','',line[3]) == 'KRAS') \
                    or (re.sub('_\d$','',line[3]) == 'NRAS') \
                    or (re.sub('_\d$','',line[3]) == 'BRAF'):
                        line[1] = int(line[1])
                        line[2] = int(line[2])
                        #line[3] = re.sub('_\d$','',line[3])
                        bed_line = IntervalColumns(*(line[0], int(line[1]), int(line[2]), line[3]))
                        intervals_list.append(bed_line)
                if len(line) == 12:
                    if (re.sub('_\d$','',line[3]) == 'EGFR') \
                    or (re.sub('_\d$','',line[3]) == 'KRAS') \
                    or (re.sub('_\d$','',line[3]) == 'NRAS') \
                    or (re.sub('_\d$','',line[3]) == 'BRAF'):
                        line[1] = int(line[1])
                        line[2] = int(line[2])
                        bed_line = IntervalColumns(*(line[0], int(line[1]), int(line[2]), line[3]))
                        intervals_list.append(bed_line)
    else:
        print("ERROR: Provide an interval list (bed format)")

coverage_result = almnt.coverage(intervals_list).sort()

coverage_list = []

for interval in coverage_result:
    coverage_list.append(float(interval[4]))

mean = np.mean(coverage_list)

collected_normalize = defaultdict(list)
collected_normalize = {
'chr' : [],
'start' : [],
'end' : [],
'name' : [],
'normalized coverage' : []
}

for interval in coverage_result:
    if interval[3] not in collected_normalize:
        collected_normalize[interval[3].encode('ascii','ignore')] = {
        'chr' : 0,
        'start' : 0,
        'end' : 0,
        'normalized coverage' : 0,
        'median' : 0
        }
    collected_normalize[interval[3].encode('ascii','ignore')]['chr'] == interval[0].encode('ascii','ignore')
    collected_normalize[interval[3].encode('ascii','ignore')]['start'] += float(interval[1])
    collected_normalize[interval[3].encode('ascii','ignore')]['end'] += float(interval[2])
    collected_normalize[interval[3].encode('ascii','ignore')]['normalized coverage'] += float(float(interval[4])/float(mean))

norm_cov = []
for key, value in collected_normalize.items():
    norm_cov.append(collected_normalize[key]['normalized coverage'])
print norm_cov


cov_list = np.asarray(norm_cov).astype(np.float)
cov_list.sort()
N = len(cov_list)
x = range(N)
plt.bar(x, cov_list)
plt.show()
