#!usr/bin/env python
# -*- coding: utf-8 -*-

#Load modules
import re
import sys
from collections import namedtuple


#BAM_FILE = "/media/partition/Haloplex/Haloplex_Test_2_Mid_February/Velona/15016513_S5.bam"
#BAM_FILE = "/media/partition/Haloplex/Haloplex_Test_2_Mid_February/Velona/15028422_S6.bam"
VCF_FILE = "/media/usb/Thesis/Haloplex_2_Mid_February/VCF/ECD_S9_04Mar2016_13_04_39_281.vcf"

INTERVALS_BED = "/media/partition/Haloplex/00100-1407755742_Regions.bed"

IntervalColumns = namedtuple('vcf', ['chr', 'pos', 'id', 'ref', 'alt', 'qual', 'filter', 'info', 'format', 'sample'])
intervals_list = []

if VCF_FILE:
    with open(VCF_FILE, "r") as fin:
        for line in fin.readlines():
            line = line.rstrip('\n')
            line = re.split(r'\t+', line.rstrip('\t'))
            info = []
            if len(line) == 10:
                line[7] = line[7].split(';')
                vcf_line = IntervalColumns(*(line[0], line[1], line[2], line[3], line[4], line[5], line[6], line[7], line[8], line[9]))
                intervals_list.append(vcf_line)
else:
    print("ERROR: Provide an interval list (bed format)")

count = 0
for line in intervals_list[1:]:
    count += 1
print intervals_list[3]
print('Number of variants in the file: %s' %count)
