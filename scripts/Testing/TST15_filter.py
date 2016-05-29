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
"""
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
                bed_line = IntervalColumns(*(line[0], line[3], int(line[1]), int(end_1st_interval), int(start_2nd_interval), int(line[2])))
                intervals_list.append(bed_line)
else:
    print("ERROR: Provide an interval list (bed format)")
print(intervals_list)
"""

#Open Bam file
samfile = pysam.AlignmentFile(BAM_FILE, "rb")

count = 0
intervals_list = [IntervalColumns('chr1', '1NRASxxE02TF004SR004', 115252121, 115252148, 115252271, 115252291)]
#Do a samtools pileup for each position
for interval in intervals_list:
    # for read in samfile.fetch(interval.chr, interval.start, interval.end):
    # TODO: do a pileup for each position
    for pileupcolumn in samfile.pileup(
        reference = interval.chr,
        start = interval.start,
        end = interval.end,
        truncate=True):
            for pileupread in pileupcolumn.pileups:
                # TODO: check if read is forward or reverse
                is_reverse = pileupread.alignment.is_reverse
                if is_reverse:
                    print("read is reverse")
                    # TODO: search reads that start with the end position
                    #if pileupread.alignment.query_alignment_sequence[0] == reference
                    # TODO: go over each base in the interval and check if there is a mismatches
                    # TODO: if there is no mismatch, write the read to a new file
                else:
                    print("read is forward")
                    # TODO: search reads thns
                    at start at the start position
                    # TODO: go over each base in the interval and check if there is a mismatches
                    # TODO: if there is no mismatch, write the read to a new file
                    """
                        #Get the reference position
                        reference_pos = pileupcolumn.reference_pos # 0 based

                        #Now iterate over each read found in the region
                        for pileupread in pileupcolumn.pileups:
                            strand_key = None

                            #Determine whether the read is in forward or reverse direction
                            is_reverse = pileupread.alignment.is_reverse
                            if is_reverse:
                                strand_key = "r"
                            else:
                                strand_key = "f"

                    intervals_list = [IntervalColumns('chr1', '1NRASxxE02TF004SR004', '115252121', '115252148', '115252271', '115252291')]

                    for interval in intervals_list:
                    only_good_reads(interval)


                                ref_query_alignment = pileupread.alignment.get_aligned_pairs(with_seq=True)
                                alignment_tuple_for_position = None
                                index_of_tuple_for_position = None
                                for index_of_tuple, alignment_tuple in enumerate(ref_query_alignment):
                                    if alignment_tuple[1] == reference_pos:
                                        alignment_tuple_for_position = alignment_tuple
                                        index_of_tuple_for_position = index_of_tuple
                                        break
                    """
