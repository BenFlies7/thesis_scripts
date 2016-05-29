#!usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function

from collections import namedtuple, defaultdict
import re
import time

import pysam

BAM_FILE = "/mnt/files-bioseq/bioseq-temporary/ben/testfiles_bamstats/assembly.bam"
INTERVALS = "/mnt/files-bioseq/bioseq-temporary/ben/testfiles_bamstats/Haloplex Cancer.bed"
REFERENCE = "/mnt/files-bioseq/bioseq-temporary/ben/testfiles_bamstats/ucsc.hg19.fasta"
AT_EDGE_THRESHOLD = 5

result_dict = defaultdict(dict)

# Parse intervals
# NOTE 1-based and half open, meaning something like 10, 11 will only get ONE position
# This is the position as I see it in IGV
IntervalColumns = namedtuple('bed', ['chr', 'start', 'end'])
intervals_list = []
with open(INTERVALS, "r") as fin:
    for line in fin.readlines():
        line = line.rstrip('\n')
        line = re.split(r'\t+', line.rstrip('\t'))
        if len(line) == 3:
            bed_line = IntervalColumns(*line)

            intervals_list.append(bed_line)

# Manual intervals for testing
# intervals_list = [IntervalColumns("chr1", 115256524, 115256525)]

# Load bam
samfile = pysam.AlignmentFile(BAM_FILE, "rb")

# start
start = time.time()

for interval in intervals_list:
    start = int(interval.start)
    end = int(interval.end)
    length = end - start
    pos_key = "%s:%s:%s" % (interval.chr, interval.start, interval.end)
    for aligned_pos in range(start, start + length):
        # NOTE: remember that python ranges are half-open meaning the last one
        # is not provided e.g. 100, 1 will only give 100
        result_dict[pos_key][aligned_pos] = {
            "fA": 0,
            "fC": 0,
            "fT": 0,
            "fG": 0,
            "fN": 0,
            "rA": 0,
            "rC": 0,
            "rT": 0,
            "rG": 0,
            "rN": 0,
            "fA_edge_start": 0,
            "fC_edge_start": 0,
            "fT_edge_start": 0,
            "fG_edge_start": 0,
            "fN_edge_start": 0,
            "rA_edge_start": 0,
            "rC_edge_start": 0,
            "rT_edge_start": 0,
            "rG_edge_start": 0,
            "rN_edge_start": 0,
            "fA_edge_end": 0,
            "fC_edge_end": 0,
            "fT_edge_end": 0,
            "fG_edge_end": 0,
            "fN_edge_end": 0,
            "rA_edge_end": 0,
            "rC_edge_end": 0,
            "rT_edge_end": 0,
            "rG_edge_end": 0,
            "rN_edge_end": 0,
            "f_read_count": 0,
            "r_read_count": 0,
        }

        for read in samfile.fetch(reference=interval.chr, start=aligned_pos - 1, end=aligned_pos - 1):
            # FIXME: I am still not sure about the above start and end points !!!
            # NOTE: aligned pos here is the same as in the interval above, so 1-based
            # NOTE: but pysam wants 0 based and halfopen, so I substract 1 everywhere I used aligned_pos
            # TODO: I only get one position here, because I have a dict storing things for each position
            # try getting over the whole interval using pysam and then only iterating over each positon
            # this should make less pysam calls and maybe be faster??

            # NOTE (BEN): You converted the start position to zero-based coordinates, but you forgot to do
            # the same for the end position. This probably causes the string index error.
            # The script works if you put "end=aligned_pos - 1" instead of "end=aligned_pos"

            """Get nuclotide in the read relative to the aligned position"""
            """Also check if nucleotide is at edge of read"""

            try:
                one_based_reference_start = read.reference_start + 1 # TESTED IN IGV -> OK !
                pos_in_string = aligned_pos - one_based_reference_start - 1
                nucleotide = read.query_sequence[pos_in_string]
            except IndexError:
                print("---------------------------------------------------")
                print("ERROR getting relative position:")
                print("Aligned pos", aligned_pos)
                print("Sequence %s, length %s" % (read.query_sequence, len(read.query_sequence)))
                print("one_based_reference_start", one_based_reference_start)
                print("my pos_in_string", pos_in_string)
                raise

            if read.is_reverse is False:
                result_dict[pos_key][aligned_pos]["f_read_count"] += 1
                nucleotide_key = "f%s" % nucleotide
            else:
                result_dict[pos_key][aligned_pos]["r_read_count"] += 1
                nucleotide_key = "r%s" % nucleotide

            result_dict[pos_key][aligned_pos][nucleotide_key] += 1


            """Is nucleotide at edge of read?"""
            edge_key = "%s_edge" % nucleotide_key

            distance_to_start = aligned_pos - read.reference_start - 1 # OK tested: == pos_in_string, how many jumps from start to position
            if distance_to_start <= AT_EDGE_THRESHOLD:
                result_dict[pos_key][aligned_pos][edge_key + "_start"] += 1

            if read.reference_end is None:
                print("Warning: No reference_end. Read is unmapped or no cigar alignment present: %s" % read)
            else:
                distance_to_end = read.reference_end - aligned_pos + 1 # OK tested !!!
                if distance_to_end <= AT_EDGE_THRESHOLD:
                    result_dict[pos_key][aligned_pos][edge_key + "_end"] += 1



        """Calculate stats for that postion"""
        # f_edge_perc = float(result_dict[pos_key][aligned_pos]["f_edge"]) / float(result_dict[pos_key][aligned_pos]["f_read_count"]) * 100.0
        # r_edge_perc = float(result_dict[pos_key][aligned_pos]["r_edge"]) / float(result_dict[pos_key][aligned_pos]["r_read_count"]) * 100.0
        # result_dict[pos_key][aligned_pos]["f_edge_perc"] = float("{0:.2f}".format(f_edge_perc))
        # result_dict[pos_key][aligned_pos]["r_edge_perc"] =  float("{0:.2f}".format(r_edge_perc))

        for aligned_pos, stats in result_dict[pos_key].items():
            print(interval.chr, aligned_pos, stats)
    print("%s done" % pos_key)

samfile.close()

print("Took {:.2f} seconds".format(time.time() - start))
