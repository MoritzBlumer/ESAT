
# coding: utf-8

# In[32]:


#!/usr/bin/python

import csv
import sys
import glob
import re

def sort_list( l ):
    """ Sorts the given iterable in alphanumerically.
 
    Required arguments:
    l -- The iterable to be sorted.
 
    """
    convert = lambda text: int(text) if text.isdigit() else text
    alphanum_key = lambda key: [convert(c) for c in re.split('([0-9]+)', key)]
    return sorted(l, key = alphanum_key)


files = sort_list(glob.glob((str(str(sys.argv[1]))+'*')))

outfilename = "merged_stats.csv"


with open (outfilename, 'w', newline='') as outfile:

    for file in files:
        title = ((str(file)).strip('stats.')).strip('.txt')

        with open(file) as f:
            content = f.readlines()
        content = [x.strip() for x in content] 


        a = csv.writer(outfile, delimiter =',')
        data = [[title,'detections','true_positives','false_negatives','false_positives'],
                ['Mobster:ALL',content[10],content[12],content[14],content[16]],
                ['MELT:ALL',content[26],content[28],content[30],content[32]],
                ['MELT:PASS_ONLY',content[38],content[40],content[42],content[44]],
                ['MELT:PASS_ONLY_HOMOZYGOTES',content[50],content[52],content[54],content[56]],
                ['MELT:PASS_ONLY_HETEROZYGOTES',content[62],content[64],content[66],content[68]]]

        a.writerows(data)

