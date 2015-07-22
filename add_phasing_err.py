#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Created on Tue Jun 30 00:26:48 2015

@author: willy
"""

import os
import argparse
import numpy as np
import glob

def read_haplotypes(msmc_input):
    """
    Returns a tuple (value1, value2, value3)
    value1 : The chromossome name 
    value 2 : a list of tuples, each containing
    the positions of the SNPS (the column 2 and 3) of the msmc input
    
    value3: a list of tuples:
    Each tuple is a diploid individual. Ex: if
    the input of msmc has 4 haplotypes, returns [(hap1, hap2), (hap3, hap4)]
    """
    a = open(msmc_input, 'r')
    lines = a.readlines()
    l1 = lines[-1].split('\t')
    value_1 = l1[0]
    hap = l1[-1].strip('\n')
    value_2 = []
    haplotypes_list = [[] for i in list(hap)]

    for l in lines:
        values = l.split('\t')
        value_2.append((values[1], values[2]))
        hap = values[3].strip('\n')
        ignore_output = [haplotypes_list[i].append(hap[i]) 
                            for i in range(len(hap))]
    value_3 = [(''.join(haplotypes_list[2*i]), ''.join(haplotypes_list[2*i+1]))
                for i in range(len(haplotypes_list)/2)]
    return (value_1, value_2, value_3)
    
def write_msmc_input(filename, value_1, value_2, value_3):
    """
    Given the type of values returned by the 'read_haplotypes' method, write
    to filename the values in the msmc input format
    """
    f = open(filename, 'w')
    for i in range(len(value_2)):
        temp_line = "{}\t{}\t".format(value_1, '\t'.join(value_2[i]))
        hap_line = ''.join([value_3[j][0][i]+value_3[j][1][i] 
                    for j in range(len(value_3))])
        f.write("{}{}\n".format(temp_line, hap_line))
    f.close()

def add_phasing_err(chrA, chrB, snp_pos, errors_positions):
    """
    Given two haplotypes with the snp positions (as in the msmc input format)
    as well as the positions of some phasing errors, this method will output
    two new chromosomes switched at positions where phasing errors happen.
    """
    snp_pos = [int(p) for p in snp_pos]
    if errors_positions == []:
        return (chrA, chrB)
    newA = []
    newB = []
    A = chrA
    B = chrB
    start_pos = 0
    curr_pos = 0
    err_index = 0
    
    while curr_pos < len(snp_pos):
        if (snp_pos[curr_pos] >= errors_positions[err_index]):
        # if phasing error
            newA.append(A[start_pos:curr_pos])
            newB.append(B[start_pos:curr_pos])
            # Invert the chromosomes
            (A, B) = (B, A)
            start_pos = curr_pos
            err_index += 1
            if err_index == len(errors_positions):
                break
        curr_pos +=1
    # Copy the rest of the haplotypes to the new chromosomes
    newA.append(A[start_pos:])
    newB.append(B[start_pos:])
    return(''.join(newA), ''.join(newB))

def generate_phasing_err_positions(seq_len, mean):
    """
    Generates the positions where phasing errors occur. It assumes 
    phasing errors are a geometrical random variable with some rate
    """
    positions = []
    curr_position = 0
    while curr_position <= seq_len:
        err_position = np.random.exponential(mean)
        err_position = np.trunc(err_position) + 1
        curr_position += err_position
        positions.append(curr_position)
    
    return [int(p) for p in positions[:-1]]
        
        
def add_phasing_err_folder(input_folder, output_folder, seq_len, err_mean):
    """
    Assuming input_folder contains some .msmcin files (the input of msmc)
    and that these files are haplotypes of size seq_len,
    for each msmcin file in the input folder we will create a new msmcin file
    in the output_folder containing some phasing errors.
    We assume that phasing errors are distributed following a Poison process
    over the annalyzed genome.
    For introducing phasing errors, we simulate values with an exponential
    distribution.
    """

    input_filenames = []
    # We take all the msmcin files in the input_folder
    for filename in glob.glob(os.path.join(input_folder, '*.txt')):
        input_filenames.append(filename)
    
    for f in input_filenames:
        (v1, v2, v3) = read_haplotypes(f)
        snp_pos = [p[0] for p in v2]
        new_v3 = []
        for v in v3:
            err_positions = generate_phasing_err_positions(seq_len, err_mean)
            new_v3.append(add_phasing_err(v[0], v[1], snp_pos, err_positions))
        output_name = f.split('/')[-1]
        f_output = os.path.join(output_folder, output_name)
        write_msmc_input(f_output, v1, v2, new_v3)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Add phasing errors to haplotypes")
    
    parser.add_argument("input_folder", help = "The input folder ")
    parser.add_argument("output_folder", help="The output folder")
    parser.add_argument("-l", "--sequence_length", help="Length of the \
                        sequence where pahsing errors will be added", type=int,
                        default = 2000000)
    parser.add_argument("-e", "--error_rate", help = "The rate of the phasing\
                        errors", type=float, default=1e6)
    args = parser.parse_args()
    
    l = args.sequence_length
    err_mean = args.error_rate
    
    add_phasing_err_folder(args.input_folder, args.output_folder, l, err_mean)
    
    print("Done")