#!usr/bin/env python
# -*- coding: utf-8 -*-

#Load modules
import pybedtools
import matplotlib.pyplot as plt
from pylab import *
import numpy as np
import scipy

bam_file_tst = '/media/partition/collected/tst15_velona/mixB/15001181_S5_velona.bam'
bam_file_hpx = '/media/partition/collected/hpx_csc_velona/15001181_S1_velona.bam'

bed_file_tst = '/media/partition/TST_15-B-manifest.bed'
bed_file_hpx = '/media/partition/00100-1407755742_Regions.bed'


almnt_tst = pybedtools.BedTool(bam_file_tst)
almnt_hpx = pybedtools.BedTool(bam_file_hpx)

regions_tst = pybedtools.BedTool(bed_file_tst)
regions_hpx = pybedtools.BedTool(bed_file_hpx)

coverage_result_tst = almnt_tst.coverage(regions_tst, d = True)
coverage_result_hpx = almnt_hpx.coverage(regions_hpx, d = True)

collected_coverage_tst = []
collected_coverage_hpx = []

for line in coverage_result_tst:
    if line[3] == '2BRAFxxE15TF001SR001':
        collected_coverage_tst.append(line[13])

for line in coverage_result_hpx:
    if line[3] == 'BRAF_1':
        collected_coverage_hpx.append(line[5])

#print collected_coverage
#print collected_position

plt.plot(collected_coverage_tst)
plt.plot(collected_coverage_hpx)


#fill_between(x, y1, y2, color='cyan')

#fig, ax1 = plt.subplots(1,1)
#ax1.fill_between(collected_coverage, 0, x)

plt.show()
