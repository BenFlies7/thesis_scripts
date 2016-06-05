#!usr/bin/env python
# -*- coding: utf-8 -*-

#Load modules
import re
import sys
import glob
from collections import namedtuple
from matplotlib_venn import venn3
import matplotlib.pyplot as plt

#directory = '/media/partition/variants'
directory = '/Volumes/KING_BEN/variants'


vcf_file_list = [f for f in glob.iglob(directory+"/*.vcf")]

set_varscan = []
set_freebayes = []
set_GATKHC =  []
set_mutect = []
set_lofreq = []

for VCF_FILE in vcf_file_list:

    basename = re.sub('.vcf$','',VCF_FILE)
    basename = re.sub(directory,'',basename)
    basename = re.sub('_freebayes','',basename)
    basename = re.sub('_varscan','',basename)
    basename = re.sub('_GATKHC','',basename)
    basename = re.sub('_mutect','',basename)
    basename = re.sub('_lofreq','',basename)

    IntervalColumns = namedtuple('vcf', ['chr', 'pos', 'ref', 'alt', 'sample'])
    intervals_list = []

    if VCF_FILE:
        with open(VCF_FILE, "r") as fin:
            for line in fin.readlines():
                line = line.rstrip('\n')
                line = re.split(r'\t+', line.rstrip('\t'))

                if len(line) == 10 or len(line)==8 or len(line)==11:
                    if line[0] != '#CHROM':
                        if line[3] == 'A' or line[3] == 'T' or line[3] == 'C' or line[3] == 'G':
                            if line[4] == 'A' or line[4] == 'T' or line[4] == 'C' or line[4] == 'G':

                                line[7] = line[7].split(';')
                                #vcf_line = IntervalColumns(*(line[0], line[1], line[3], line[4], basename))
                                #intervals_list.append(vcf_line)
                                variant_ID = basename+'_'+line[0]+'_'+line[1]+'_'+line[3]+'_'+line[4]
                                #print variant_ID

                                if 'freebayes' in VCF_FILE:
                                    set_freebayes.append(variant_ID)
                                if 'varscan' in VCF_FILE:
                                    set_varscan.append(variant_ID)
                                if 'GATKHC' in VCF_FILE:
                                    set_GATKHC.append(variant_ID)
                                if 'mutect' in VCF_FILE:
                                    set_mutect.append(variant_ID)
                                if 'lofreq' in VCF_FILE:
                                    set_lofreq.append(variant_ID)


print set_lofreq
print('#####')
print set_mutect



set_freebayes = set(set_freebayes)
set_varscan = set(set_varscan)
set_GATKHC = set(set_GATKHC)
set_mutect = set(set_mutect)
set_lofreq = set(set_lofreq)

venn3([set_varscan, set_mutect, set_GATKHC], ('VarScan2', 'MuTect1.1.7', 'GATK_HaplotypeCaller'))
plt.show()
