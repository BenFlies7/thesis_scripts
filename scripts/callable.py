#!/usr/bin/python
# -*- coding: utf-8 -*-

from collections import namedtuple
import re

INTERVALS_BED = '/home/ben/callable_status.bed'

IntervalColumns = namedtuple('bed', ['chr', 'start', 'end', 'status'])
intervals_list = []

if INTERVALS_BED:
    with open(INTERVALS_BED, "r") as fin:
        for line in fin.readlines():
            line = line.rstrip('\n')
            line = re.split(r'\t+', line.rstrip('\t'))

            # For Agilent Haloplex BED files
            if len(line) == 4:
                line[1] = int(line[1])
                line[2] = int(line[2])
                bed_line = IntervalColumns(*(line[0], int(line[1]), int(line[2]), line[3]))
                intervals_list.append(bed_line)

counter = 0
for interval in intervals_list:
    if interval[3] == 'CALLABLE':
        counter += 1

print counter
