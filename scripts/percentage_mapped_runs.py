#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
This script sorts all .bam files of a given input directory
"""

import sys
import os

def aligned_counts(ifile1):
    '''
    count how many alignments are aligned back to genome,
    ifile1 is a sorted bam file
    '''
    import HTSeq
    sortedbamfile= HTSeq.BAM_Reader(ifile1)
    aligned_counts=0
    unaligned_counts=0
    for almnt in sortedbamfile:
        if almnt.aligned:
            aligned_counts+= 1
        else:
            unaligned_counts+=1
    sum = float(aligned_counts) + float(unaligned_counts)
    ratio = (aligned_counts / sum) * float(100)
    #print "number of aligned tags of %s is %d " % (ifile1, aligned_counts)
    #print "number of unaligned tags of %s is %d "% (ifile1, unaligned_counts)
    return ratio

directories = sys.argv[1:]

for arg in directories:

    os.chdir("/home")
    path_fastq = arg
    os.chdir(os.path.dirname(os.path.abspath(arg)))
    files = os.listdir(path_fastq)
    almnt = 0
    #search fastq files:
    f = [os.path.join(root,name)
                 for root, dirs, files in os.walk(path_fastq)
                 for name in files
                 if name.endswith(".bam")]

    f.sort()

    for file in f:
        print file
        aligned = aligned_counts(file)
        print aligned
        almnt += aligned


    percentage = almnt / len(f)
    print "%s : average percentage mapped: %s" %(arg,percentage)
