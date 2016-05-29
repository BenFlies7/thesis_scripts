#!usr/bin/env python
# -*- coding: utf-8 -*-

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


"""
This script calculates quality metrics of alignment using pysam

BAM Index required
Give vcf or bed file for positions to analyse
If both are given, vcf file is taken
"""

# TODO: get base on reference sequence

REFERENCE = "/media/partition/hg19/ucsc.hg19.fasta"
BAM_FILE = "/media/partition/TST15_Test_1_Early_February/Base_Space/EGFR_15001181/Libraries/15001181_S1.bam"
INTERVALS_BED = "/media/partition/TST15_Test_1_Early_February/TST_15-A-manifest.bed"
VCF_FILE = None#"/home/christophe/workspace/bioseq-scripts/testfiles_bamstats/variants_short.vcf"

AT_EDGE_THRESHOLD = 5
MAP_QUAL_THRESHOLD = 0 # reads with quality below this value will be excluded, reads with qualities equal to this value will be included
MAX_DEPTH = 8000 # pysam has default of 8000

THREADS = 3


def calculate_for_column(reference, zero_based_pos):
    result_dict = {
        "fA": 0,
        "fC": 0,
        "fT": 0,
        "fG": 0,
        "fN": 0,
        "fD": 0,
        "rA": 0,
        "rC": 0,
        "rT": 0,
        "rG": 0,
        "rN": 0,
        "rD": 0,
        "fI": 0,
        "rI": 0,
        "fA_edge_start": 0,
        "fC_edge_start": 0,
        "fT_edge_start": 0,
        "fG_edge_start": 0,
        "fN_edge_start": 0,
        "fD_edge_start": 0,
        "rA_edge_start": 0,
        "rC_edge_start": 0,
        "rT_edge_start": 0,
        "rG_edge_start": 0,
        "rN_edge_start": 0,
        "rD_edge_start": 0,
        "fA_edge_end": 0,
        "fC_edge_end": 0,
        "fT_edge_end": 0,
        "fG_edge_end": 0,
        "fN_edge_end": 0,
        "fD_edge_end": 0,
        "rA_edge_end": 0,
        "rC_edge_end": 0,
        "rT_edge_end": 0,
        "rG_edge_end": 0,
        "rN_edge_end": 0,
        "rD_edge_end": 0,
        "f_read_count": 0,
        "r_read_count": 0,
        "total_read_count": None,
        "base_qualities": [],
        "base_qualities_average": None,
        "mapping_qualities": [],
        "mapping_qualities_average": None
    }


    samfile = pysam.AlignmentFile(BAM_FILE, "rb")
    for pileupcolumn in samfile.pileup(
        reference=reference,
        start=zero_based_pos,
        end=zero_based_pos + 1,
        max_depth=MAX_DEPTH,
        truncate=True):


        """
        Get reads covering position
        Truncate=True gives us only the column asked not all columns
        covered by the reads for the position
        """
        reference_pos = pileupcolumn.reference_pos # 0 based
        for pileupread in pileupcolumn.pileups:
            strand_key = None
            nucleotide_key = None
            ref_query_alignment = None


            """
            Quality filter
            """
            if pileupread.alignment.mapping_quality < MAP_QUAL_THRESHOLD:
                continue

            result_dict['mapping_qualities'].append(pileupread.alignment.mapping_quality)

            """
            Forward / Reverse
            """
            is_reverse = pileupread.alignment.is_reverse
            if is_reverse:
                strand_key = "r"
            else:
                strand_key = "f"
            read_count_key = strand_key + "_read_count"
            result_dict[read_count_key] += 1


            """
            Get alignment detail for my position
            """
            # from aligned_pairs, get the tuple of the current position
            # pileupread.query_position (first memeber of tuple) is None if the read
            # has a deletion compared to reference at that position.
            # e.g. (95, None, None), (96, 7577626, 'G'), (None, 7577627, 'C'), (None, 7577628, 'A'), (97, 7577629, 'a'), (98, 7577630, 'T'), (99, 7577631, 'C')]
            #       (None, 7577627, 'C') read has deletion compared to ref
            #       (95, None, None) read has insertion compared to ref
            # For deletions i take the index to later look for insertions around
            # the position and calculate edge distances.
            # So if no deletion, index_of_tuple_for_position == query_position
            ref_query_alignment = pileupread.alignment.get_aligned_pairs(with_seq=True)
            alignment_tuple_for_position = None
            index_of_tuple_for_position = None
            for index_of_tuple, alignment_tuple in enumerate(ref_query_alignment):
                if alignment_tuple[1] == reference_pos:
                    alignment_tuple_for_position = alignment_tuple
                    index_of_tuple_for_position = index_of_tuple
                    break

            """
            Get base in reference
            """
            # TODO


            """
            Get nucleotide (or deletion) for that read at the given reference position
            """
            if not pileupread.is_del:
                # query position is None if is_del is True.
                # here, the read has a deletion compared to reference
                base_in_read = pileupread.alignment.query_sequence[pileupread.query_position]
                nucleotide_key = strand_key + base_in_read
                result_dict["base_qualities"].append(pileupread.alignment.query_qualities[pileupread.query_position])
            else:
                # Deletions at a postion do not count towards coverage
                nucleotide_key = strand_key + "D"
                result_dict[read_count_key] -= 1
                # print("* %s has deletion compared to reference" % pileupread.alignment.query_name)

            result_dict[nucleotide_key] += 1


            """
            Does read have insertion around the position?

            It could be that the position is in the middle of the read but
            all subsequent bases are not aligned to ref. This is not an insertion
            but a partial alignment
            """
            before = None
            after = None

            if index_of_tuple_for_position - 1 >= 0:
                before = ref_query_alignment[index_of_tuple_for_position - 1]
            if index_of_tuple_for_position + 1 < len(ref_query_alignment):
                after = ref_query_alignment[index_of_tuple_for_position + 1]

            if before and before[1] is None and any(x[1] is not None for x in ref_query_alignment[0:index_of_tuple_for_position - 1]):
                result_dict[strand_key + "I"] += 1
                # print("* %s has insertion before reference pos" % pileupread.alignment.query_name, before, index_of_tuple_for_position)
            if after and after[1] is None and any(x[1] is not None for x in ref_query_alignment[index_of_tuple_for_position + 1: len(ref_query_alignment)]):
                result_dict[strand_key + "I"] += 1
                # print("* %s has insertion after reference pos" % pileupread.alignment.query_name, after, index_of_tuple_for_position)


            """
            Is nucleotide at edge of read?
            """
            edge_key = "%s_edge" % nucleotide_key
            # For num_to_end, compare to actual read length, not the alignment. They can differ significantly
            num_to_start = index_of_tuple_for_position
            num_to_end = len(pileupread.alignment.query_sequence) - index_of_tuple_for_position - 1
            if num_to_start <= AT_EDGE_THRESHOLD:
                result_dict[edge_key + "_start"] += 1
            if num_to_end <= AT_EDGE_THRESHOLD:
                result_dict[edge_key + "_end"] += 1


    """
    Clean up for that position
    """
    samfile.close()

    """
    Calculate stats for that position
    """
    result_dict["total_read_count"] = result_dict["f_read_count"] + result_dict["r_read_count"]
    result_dict['base_qualities_average'] = numpy.mean(result_dict['base_qualities'])
    result_dict['mapping_qualities_average'] = numpy.mean(result_dict['mapping_qualities'])
    result_dict["base_qualities_average"] = float("{0:.1f}".format(result_dict["base_qualities_average"]))
    result_dict["mapping_qualities_average"] = float("{0:.1f}".format(result_dict["mapping_qualities_average"]))

    """
    Histogram of mapping qualities
    """
    map_qual_count = {n: float(result_dict['mapping_qualities'].count(n) / float(len(result_dict['mapping_qualities']))) for n in result_dict['mapping_qualities']}
    map_qual_count = OrderedDict(sorted(map_qual_count.items(), key = itemgetter(0)))
    print(map_qual_count)
    plt.bar(map_qual_count.keys(), map_qual_count.values())
    plt.show()

    """
    Histogram of base qualities
    """
    base_qual_count = {n: float(result_dict['base_qualities'].count(n) / float(len(result_dict['base_qualities']))) for n in result_dict['base_qualities']}
    base_qual_count = OrderedDict(sorted(base_qual_count.items(), key = itemgetter(0)))
    print(base_qual_count)
    plt.bar(base_qual_count.keys(), base_qual_count.values())
    plt.show()

    del result_dict['base_qualities']
    del result_dict['mapping_qualities']

    result = "Pos(1b) %s:%s A %s %s C %s %s T %s %s G %s %s N %s %s INS %s %s DEL %s %s COUNT %s %s" % (
        reference,
        zero_based_pos + 1,
        result_dict['fA'],
        result_dict['rA'],
        result_dict['fC'],
        result_dict['rC'],
        result_dict['fT'],
        result_dict['rT'],
        result_dict['fG'],
        result_dict['rG'],
        result_dict['fN'],
        result_dict['rN'],
        result_dict['fI'],
        result_dict['rI'],
        result_dict['fD'],
        result_dict['rD'],
        result_dict['f_read_count'],
        result_dict['r_read_count']
        )
    print(result)

    return ("%s:%s" % (reference, zero_based_pos + 1), result_dict)

"""
Main program    base_qual_count = sorted(base_qual_count.iteritems())

http://sebastianraschka.com/Articles/2014_multiprocessing_intro.html
"""

"""
Get intervals (or positions) to calculate stats
"""
IntervalColumns = namedtuple('bed', ['chr', 'start', 'end'])
intervals_list = []

if VCF_FILE:
    # NOTE: vcf coordinates are 1-based, but pysam converts to 0-based
    bcf_in = pysam.VariantFile(VCF_FILE)  # auto-detect input format
    for x in bcf_in.fetch():
        interval  = IntervalColumns(x.contig, int(x.start), int(x.start) + 1) # It's correct with this
        intervals_list.append(interval)

elif INTERVALS_BED:
    # NOTE BED format is 0 based (open interval)
    with open(INTERVALS_BED, "r") as fin:
        for line in fin.readlines():
            line = line.rstrip('\n')
            line = re.split(r'\t+', line.rstrip('\t'))
            if len(line) == 3:
                line[1] = int(line[1])
                line[2] = int(line[2])
                bed_line = IntervalColumns(*line)
                intervals_list.append(bed_line)
else:
    print("ERROR: Provide either an interval list (bed format) or a vcf file")


intervals_list = [IntervalColumns("chr17", 7577625, 7577626)] # end of reads

results = []
overall_result = defaultdict(dict)
start = time.time()

for interval in intervals_list:
    """Loop over all intervals"""
    # NOTE: we use 1-based closed interval to indicate region, intervals_list contains 0-based coordinates
    region_key = "%s:%s:%s" % (interval.chr, interval.start + 1, interval.end)
    print("Interval: %s (1 based, inclusive) added" % region_key)
    for zero_based_pos in range(interval.start, interval.end):
        results.append(calculate_for_column(interval.chr, zero_based_pos))

for result in results:
    overall_result[result[0]] = result[1]
# print(json.dumps(overall_result))
print(time.time() - start)
