#!usr/bin/env python
# -*- coding: utf-8 -*-

import pybedtools
import time
import urllib
import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec

#BAM_file = "/media/partition/Haloplex_Test_1_Late_January/Velona/15038519_S3.bam"
BED_file = "/media/partition/Haloplex/00100-1407755742_Regions.bed"
#BAM_FILE = "/media/partition/Haloplex_Test_1_Late_January/Velona/15038519_S3.bam"
#INTERVALS_BED = "/media/partition/Haloplex_Test_1_Late_January/00100-1407755742_Regions.bed"

#BAM_FILE = "/media/partition/Haloplex_Test_1_Late_January/Velona/15038519_S3.bam"
#INTERVALS_BED = "/media/partition/Haloplex_Test_1_Late_January/00100-1407755742_Regions.bed"

#BAM_FILE = "/media/partition/TST15_Test_1_Early_February/Base_Space/BRAF_15018040/Libraries/15018040_S1.bam"
#INTERVALS_BED = "/media/partition/TST15_Test_1_Early_February/TST_15-A-manifest.bed"
#BAM_FILE = "/media/partition/15001181_S1.bam"
#BAM_file = "/media/partition/Haloplex_Test_2_Mid_February/Velona/15028422_S6.bam"
BAM_file = "/media/partition/Haloplex/Haloplex_Test_2_Mid_February/Velona/ECD_S9.bam"
#BAM_FILE = "/media/partition/Haloplex_Test_2_Mid_February/Velona/15051669_S1.bam"
#BAM_FILE = "/media/partition/Haloplex_Test_2_Mid_February/Velona/15020056_S2.bam"
#BAM_FILE = "/media/partition/Haloplex_Test_2_Mid_February/Velona/15010800_S3.bam"
#BAM_FILE = "/media/partition/Haloplex_Test_2_Mid_February/Velona/15039121_S7.bam"

def get_coverage(bam_file,bed_file):
    """
    Compute the coverage from an input BAM file,
    by only taking into consideration the regions
    provided in the BED file. Plot the resulting
    histogram & boxplot.
    """
    almnt = pybedtools.BedTool(bam_file)
    regions = pybedtools.BedTool(bed_file)
    print "Calculating coverage over regions ..."
    sys.stdout.flush()
    t0 = time.time()
    coverage_result = almnt.coverage(regions).sort()
    print(coverage_result)
    print coverage_result
    coverage_array = np.array([i[-4] for i in coverage_result], dtype=float)
    print coverage_array
    t1 = time.time()
    print "completed in %.2fs" % (t1-t0)
    sys.stdout.flush()
    fig = plt.figure(figsize=(12, 6))
    gs = gridspec.GridSpec(1, 2, width_ratios=[2.5, 1])
    ax1 = plt.subplot(gs[0])
    ax1.set_title("Coverage Histograms")
    bins = len(coverage_array)
    ax1.hist(coverage_array,bins,orientation="horizontal")
    ax1.set_xlabel("Count")
    ax1.set_ylabel("Coverage (X)")
    ax2 = plt.subplot(gs[1])
    ax2.set_title("Coverage Distribution")
    ax2.boxplot(coverage_array)
    ax2.set_ylabel("Coverage (X)")
    plt.show()

get_coverage(BAM_file, BED_file)
