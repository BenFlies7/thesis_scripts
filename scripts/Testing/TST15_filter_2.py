#!usr/bin/env python
# -*- coding: utf-8 -*-

# TODO: do samtools pileup
# TODO: get reads covering each individual position
# TODO: check if forward or reverse read, if forward, consider first element of up_down_tuple
# if reverse, consider second element of up_down_tuple
# TODO: check if there are any mismatches, if not, write the read in a new file






from __future__ import print_function

from collections import namedtuple, defaultdict, Counter, OrderedDict
from operator import itemgetter
import json
import re
import time

import matplotlib.pyplot as plt
import numpy
import pysam
import scipy.stats.stats as st


REFERENCE = "/media/partition/hg19/ucsc.hg19.fasta"
BAM_FILE = "/media/partition/TST15_Test_1_Early_February/Base_Space/BRAF_15018040/Libraries/15018040_S1.bam"
INTERVALS_BED = "/media/partition/TST15_Test_1_Early_February/TST_15-A-manifest.bed"

IntervalColumns = namedtuple('bed', ['chr', 'oligo_name', 'start', 'end_first_interval', 'start_second_interval', 'end'])
intervals_list = []

if INTERVALS_BED:
    with open(INTERVALS_BED, "r") as fin:
        for line in fin.readlines():
            line = line.rstrip('\n')
            line = re.split(r'\t+', line.rstrip('\t'))
            prc = ()
            if len(line) == 12:
                up_down = line[10].split(',')
                (b1,b2) = (int(up_down[0]),int(up_down[1]))
                up_down_tuple = (b1,b2)
                line[1] = int(line[1])
                line[2] = int(line[2])
                end_1st_interval = line[1] + b1
                start_2nd_interval = line[2] - b2
                bed_line = IntervalColumns(*(line[0], line[3], line[1], end_1st_interval, start_2nd_interval, line[2]))
                intervals_list.append(bed_line)
else:
    print("ERROR: Provide an interval list (bed format)")
print(intervals_list)

def only_good_reads(reference, start, end):

    samfile = pysam.AlignmentFile(BAM_FILE, "rb")

    for pileupcolumn in samfile.pileup(
        reference=reference,
        start=start,
        end=end + 1,
        truncate=True):

        reference_pos = pileupcolumn.reference_pos # 0 based
        for pileupread in pileupcolumn.pileups:
            strand_key = None
            is_reverse = pileupread.alignment.is_reverse
            if is_reverse:
                strand_key = "r"
            else:
                strand_key = "f"

                ref_query_alignment = pileupread.alignment.get_aligned_pairs(with_seq=True)
                alignment_tuple_for_position = None
                index_of_tuple_for_position = None
                for index_of_tuple, alignment_tuple in enumerate(ref_query_alignment):
                    if alignment_tuple[1] == reference_pos:
                        alignment_tuple_for_position = alignment_tuple
                        index_of_tuple_for_position = index_of_tuple
                        break
