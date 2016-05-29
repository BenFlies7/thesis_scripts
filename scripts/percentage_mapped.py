#!/usr/bin/python
#!/usr/bin/python
# -*- coding: utf-8 -*-


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

file = '/media/usb/Thesis/Haloplex_2_Mid_February/SURECALL/15010800_S3_L001_R1_001.fastq_04Mar2016_11_53_24_15_Sorted.bam'

aligned = aligned_counts(file)
print "File %s : %s" %(file,aligned)
