# -*- coding: utf-8 -*-
"""
Created on Wed Mar 24 18:25:37 2021

@author: stuart
"""

import pandas as pd
import numpy as np
import os
from Bio import SeqIO

# Set working directory
os.chdir("C:\\Users\\stuart\\Documents\\SubalpineMyco")

# Read phd file with each line as an element in a list
phd = open("Data\\raw_sequence_files\\cleaned_77-ITS1F.phd.1").readlines()

# Create empty list to store trace indices of base calls
stringdex=[]

# Create empty list to store called bases
base_call = []

# Start a counter to locate DNA sequence lines
counter=0

# Loop through lines adding trace indices to list
for line in phd:
    if line=='END_DNA\n': break
    if line=='BEGIN_DNA\n':
        counter+=1
        continue
    if counter==1:
        base_call.append(line.split(' ')[0])
        stringdex.append(int(line.split(' ')[2].rstrip('\n')))

# Read ab1 file
abif = SeqIO.read("Data\\raw_sequence_files\\cleaned_77-ITS1F.ab1", "abi")

# Create pandas dataframe of relevant channels
channel_set = {'Chnl9': abif.annotations["abif_raw"]["DATA9"],
               'Chnl10': abif.annotations["abif_raw"]["DATA10"],
               'Chnl11': abif.annotations["abif_raw"]["DATA11"],
               'Chnl12': abif.annotations["abif_raw"]["DATA12"]}
trace_df = pd.DataFrame(channel_set,
                        columns = ['Chnl9', 'Chnl10', 'Chnl11', 'Chnl12'])

# Subset trace dataframe to basecall indices
trace_df = trace_df.iloc[stringdex]

# Create column showing channel with highest signal per base
trace_df['max_chnl'] = trace_df.idxmax(axis="columns")

# Add basecalls to trace_df
trace_df['Basecall'] = base_call

# Define function for making new base call using trace
def new_call (row):
    traces = [row['Chnl9'], row['Chnl10'], row['Chnl11'], row['Chnl12']]
    traces.sort()
    if traces[-2] != 0 and traces[-1] / traces[-2] < 2:
        return 'N'
    else:
        if row['max_chnl'] == 'Chnl9':
            return 'G'
        if row['max_chnl'] == 'Chnl10':
            return 'A'
        if row['max_chnl'] == 'Chnl11':
            return 'T'
        if row['max_chnl'] == 'Chnl12':
            return 'C'

# Make new base calls
trace_df['new_call'] = trace_df.apply (lambda row: new_call(row), axis=1)

# Sum number of discrepancies between new and original calls
sum(trace_df['Basecall'] != trace_df['new_call'])

# View discrepancies
diffs = np.where(trace_df['Basecall'] != trace_df['new_call'])
diffs_df = trace_df.iloc[diffs]
