#!usr/bin/env python
# -*- coding: utf-8 -*-

#Load modules
import re
import sys
import glob
from collections import namedtuple
from matplotlib_venn import venn3
import matplotlib.pyplot as plt

directory = '/media/partition/variants'

vcf_file_list = [f for f in glob.iglob(directory+"/*.vcf")]
print vcf_file_list
caller_list = ['varscan', 'freebayes', 'GATKHC', 'mutect', 'lofreq']

collected = {
'freebayes' : [],
'varscan' : [],
'GATKHC' : [],
'mutect' : [],
'lofreq' :  []
}

variant_list = []

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

                if len(line) == 10 or len(line) == 8 or len(line) == 11:
                    if line[0] != '#CHROM':
                        if line[3] == 'A' or line[3] == 'T' or line[3] == 'C' or line[3] == 'G':
                            if line[4] == 'A' or line[4] == 'T' or line[4] == 'C' or line[4] == 'G':
                                line[7] = line[7].split(';')
                                #vcf_line = IntervalColumns(*(line[0], line[1], line[3], line[4], basename))
                                #intervals_list.append(vcf_line)
                                variant_ID = basename+'_'+line[0]+'_'+line[1]+'_'+line[3]+'_'+line[4]
                                #print variant_ID

                                variant_list.append(variant_ID)

                                if 'freebayes' in VCF_FILE:
                                    collected['freebayes'].append(variant_ID)
                                if 'varscan' in VCF_FILE:
                                    collected['varscan'].append(variant_ID)
                                if 'GATKHC' in VCF_FILE:
                                    collected['GATKHC'].append(variant_ID)
                                if 'mutect' in VCF_FILE:
                                    collected['mutect'].append(variant_ID)
                                if 'lofreq' in VCF_FILE:
                                    collected['lofreq'].append(variant_ID)

for key1, value1 in collected.items():
    res = 0
    coll = len(value1)
    for key2, value2 in collected.items():
        counter = 0
        for value in value1:
            if value in value2:
                counter +=1
        if coll != 0:
            res = counter / float(coll)
        else:
            res = 0
        print('%s vs %s : %s') %(key1, key2, res)
