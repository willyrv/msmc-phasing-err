#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 21 20:31:17 2015

@author: willy
"""

import random

HAP = "ACAC"
LEN = 2000000
TOTAL_SNPS = 100

if __name__ == "__main__":
    pos = [random.randint(1, LEN) for i in range(TOTAL_SNPS)]
    pos.sort()
    for i in range(1, 11):
        curr_pos = 0
        f = open("./phased_data/data{}.txt".format(i), 'w')
        for j in range(len(pos)):
            diff = pos[j]-curr_pos
            f.write("{}\t{}\t{}\t{}\n".format(i, pos[j], diff, HAP))
            curr_pos = pos[j]
        f.close()

    