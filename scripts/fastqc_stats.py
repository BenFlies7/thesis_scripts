#!/usr/bin/python
# -*- coding: utf-8 -*-

#import modules
from collections import namedtuple
import os
import re
import shutil
from shutil import copyfile

WarningColumns = namedtuple('warnings', ['warning', 'analysis'])
warnings_list = []

ResultColumns = namedtuple('results', ['parameter', 'result'])
results_list = []

OverrepColumns = namedtuple('overrepresented_sequences', ['sequence', 'count', 'percentage', 'hit'])
overrep_list = []

'''
directory = '/Desktop/15001181_S1_L001_R1_001_fastqc'

summary = os.path.join(directory,'summary.txt')

data = os.path.join(directory,'fastqc_data.txt')

images = os.path.join(directory,'Images')
'''

summary = 'summary.txt'

results = 'fastqc_data.txt'

with open(summary, "r") as f:
    for line in f.readlines():
        line = line.rstrip('\n')
        line = re.split(r'\t+', line.rstrip('\t'))
        warning_line = WarningColumns(*(line[0], line[1]))
        warnings_list.append(warning_line)

with open(results,"r") as res:
    for line in res.readlines():
        line = line.rstrip('\n')
        line = re.split(r'\t+', line.rstrip('\t'))

        if len(line) == 2:
            if line[0]=='Total Sequences' or line[0]=='Filtered Sequences' or line[0]=='Sequence length' or line[0]=='%GC':
                result_line = ResultColumns(*(line[0], line[1]))
                results_list.append(result_line)

        if len(line) == 4:
            if line[0].startswith('A') or line[0].startswith('T') or line[0].startswith('C') or line[0].startswith('G'):
                if float(line[2]) > 1:
                    overrep_line = OverrepColumns(*(line[0], line[1], line[2], line[3]))
                    overrep_list.append(overrep_line)


shutil.copy('Images/per_base_quality.png', 'per_base_quality.png')
os.remove('fastqc_report.html')
shutil.rmtree('Icons')
shutil.rmtree('Images')

print warnings_list
print results_list
print overrep_list
