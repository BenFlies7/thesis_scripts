#!usr/bin/env python
# -*- coding: utf-8 -*-

#Load modules

directory = '/media/partition/variants'

vcf_file_list = [f for f in glob.iglob(directory+"/*.vcf")]

VARIANTS_LIST = '/media/partition/variants/variants.txt'
