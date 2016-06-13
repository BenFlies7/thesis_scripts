#! usr/bin/env python
# -*- coding: utf-8 -*-

"""
This script calculates the percentage of reads above a certain
threshold for each position and plots the result
"""

#import modules
from Bio import SeqIO
from collections import defaultdict
import numpy as np
import matplotlib.pyplot as plt
import pylab as p
import glob
import gzip

### Input
directory = '/media/partition/'

fastq_file_list = [f for f in glob.iglob(directory+"/*.fastq.gz")]

### Open file & parse over it
for file in fastq_file_list:
    handle = gzip.open(file)
    recs = SeqIO.parse(handle,"fastq")
    d1 = defaultdict(list)
    for rec in recs:
        pos = 0
        for i,qual in enumerate(rec.letter_annotations['phred_quality']):
            d1[pos].append(qual)
            pos = pos + 1
    means = []
    mean = 0
    for key, value in d1.items():
        mean = np.mean(value)
        means.append(mean)
    print means

    x = np.arange(len(means))

    plt.plot(x, means)

plt.show()
