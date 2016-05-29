#!usr/bin/env python
# -*- coding: utf-8 -*-

import pysam
import glob
import matplotlib.pyplot as plt
import numpy as np

DIRECTORY = '/media/partition/Haloplex/Haloplex_Test_2_Mid_February/Velona'
#DIRECTORY = '/media/partition/TST15/TST15_Test_1_Early_February/Base_Space'

bam_file_list = [f for f in glob.iglob(DIRECTORY+"/*.bam")]

ratio_list = []
for file in bam_file_list:
    almnt = pysam.AlignmentFile(file,'rb')
    mapped = almnt.mapped
    unmapped = almnt.unmapped
    total = mapped + unmapped
    ratio = mapped / float(total)
    ratio_list.append(ratio)

mean = np.mean(ratio_list)

print mean

N = len(ratio_list)
x = range(N)
plt.bar(x, ratio_list)
plt.show()
