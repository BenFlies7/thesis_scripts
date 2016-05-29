#!usr/bin/env python
# -*- coding: utf-8 -*-

import pysam

BAM_FILE = "/home/ben/work_directory/assembly.bam"
samfile = pysam.AlignmentFile(BAM_FILE, "rb")

REFERENCE = "/media/partition/hg19/ucsc.hg19.fasta"

AT_EDGE_THRESHOLD = 5
MAP_QUAL_THRESHOLD = 0 # reads with quality below this value will be excluded, reads with qualities equal to this value will be included
MAX_DEPTH = 8000 # pysam has default of 8000

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

            print result_dict['mapping_qualities']
