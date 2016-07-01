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
import seaborn as sns

directory = '/media/partition/test/'

fastq_file_list = [f for f in glob.iglob(directory+"*.fastq.gz")]

print fastq_file_list

for file in fastq_file_list:
    handle = gzip.open(file)
    recs = SeqIO.parse(handle,"fastq")

    means = []
    counter = 0

    for rec in recs:
        counter += 1
        values = []
        #if counter % 5000 == 0:
        #    print counter
        for i,qual in enumerate(rec.letter_annotations['phred_quality']):
            values.append(qual)
        mean = np.mean(values)
        means.append(mean)

    bad_counter = 0
    for mean in means:
        if mean < float(35):
            bad_counter += 1

    means.append(40)

    print('Total reads: %s') %(counter)
    print('Reads >Q35: %s') %(bad_counter)

    binsize = 120
    y,binEdges = np.histogram(means, binsize, normed=True, density=True)
    unity_y = y / y.sum()

    #Consider only the bin centers
    bincenters = 0.5*(binEdges[1:]+binEdges[:-1])

    #Align the bin centers
    p.plot(bincenters,unity_y,'-')

plt.ylabel('Frequency')
plt.xlabel('Quality (Q score)')
plt.show()
