#!/usr/bin/python
# -*- coding: utf-8 -*-

#import modules
import os
import pysam
import matplotlib.pyplot as plt
import numpy

#import bam file and index if necessary
bamfile = pysam.Samfile("assembly.bam","rb")
if not os.path.exists("assembly.bam.bai"):
    pysam.index("assembly.bam")
"""
#count reads
sizeCount = {}
for read in bamfile.fetch():
    sizeCount[len(read.seq)] = sizeCount.get(len(read.seq), 0) + 1

N = 0
for key, count in sizeCount.iteritems():
    N += count
print "Total: %d reads" % N

print bamfile.references
print bamfile.lengths

N = 0
for read in bamfile.fetch("chr7", 55241674,55241712):
    N += 1
print N, "reads in region"
"""

x = []
y = []

for column in bamfile.pileup("chr7", 55242415, 55242515, max_depth = 1000000):
    x.append(column.pos)
    n = 0
    for read in column.pileups:
        if (not read.is_del):
            n += 1
    y.append(n)
plt.figure(figsize=(15, 5))
plt.plot(x, y, 'b')
plt.plot([x[0], x[-1]], [numpy.mean(y[50:-50]), numpy.mean(y[50:-50])], ':r')
plt.title("Coverage Plot")
plt.xlabel("Position (bp)")
plt.ylabel("Coverage (X)")
plt.show()
