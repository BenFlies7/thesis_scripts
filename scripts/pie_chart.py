#!usr/bin/env python
# -*- coding: utf-8 -*-

#Load modules
import re
import sys
import glob
from collections import namedtuple, Counter
from matplotlib_venn import venn3
import matplotlib.pyplot as plt

directory = '/media/partition/variants'

vcf_file_list = [f for f in glob.iglob(directory+"/*.vcf")]

caller_list = ['varscan', 'freebayes', 'GATKHC', 'mutect', 'lofreq']

collected = {
'freebayes' : [],
'varscan' : [],
'GATKHC' : [],
'mutect' : [],
'lofreq' :  []
}

variants_list = []

variants_collected = []

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

                                if variant_ID not in variants_collected:
                                    variants_collected.append(variant_ID)

                                variants_list.append(variant_ID)

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

five = 0
four = 0
three = 0
two = 0
one = 0

sizes = []

for item in variants_collected:
    a = variants_list.count(item)

    if a ==5:
        five +=1
        print item
    if a == 4:
        four += 1
    if a == 3:
        three +=1
    if a == 2:
        two +=1
    if a == 1:
        one +=1

def make_autopct(sizes):
    def my_autopct(pct):
        total = sum(sizes)
        val = int(round(pct*total/100))
        return '{p:.2f}%  ({v:d})'.format(p=pct,v=val)
    return my_autopct

sizes.append(five)
sizes.append(four)
sizes.append(three)
sizes.append(two)
sizes.append(one)

labels = '5 callers', '4 callers', '3 callers', '2 callers', '1 caller'
#colors = ['#191970','#001CF0','#0038E2','#0055D4','#0071C6','#008DB8','#00AAAA','#00C69C','#00E28E','#00FF80',]
colors = ['yellow', 'green','firebrick','orange','slategrey']
#explode = (0.05, 0.05, 0.05, 0.05, 0)

#plt.pie(sizes, explode=explode, labels=labels, colors=colors,autopct=make_autopct(sizes))

plt.pie(sizes, labels=labels, colors=colors,autopct=make_autopct(sizes))

plt.axis('equal')
plt.show()
