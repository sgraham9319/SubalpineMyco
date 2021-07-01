# -*- coding: utf-8 -*-
"""
Created on Wed Mar 24 18:25:37 2021

@author: stuart
"""

import pandas as pd
# import numpy as np
import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


# Define function for making new base call using trace columns
def new_call(row):
    traces = [row['Chnl9'], row['Chnl10'], row['Chnl11'], row['Chnl12']]
    traces.sort()
    if traces[-2] != 0 and traces[-1]/traces[-2] < 2:
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


# Specify sequence id
# sample = "79"

# Set working directory
os.chdir("C:\\Users\\stuart\\Documents\\SubalpineMyco\\Data\\bioassay")

# Loop through samples
for i in range(1, 90):
    sample = str(i)
    if os.path.exists("raw_data\\sequence_files\\cleaned_" + sample +
                      "-ITS1F.phd.1"):

        # Read phd file with each line as an element in a list
        phd = open("raw_data\\sequence_files\\cleaned_" + sample +
                   "-ITS1F.phd.1").readlines()

        # Create empty list to store trace indices of base calls
        stringdex = []

        # Create empty list to store called bases
        base_call = []

        # Start a counter to locate DNA sequence lines
        counter = 0

        # Loop through lines adding trace indices to list
        for line in phd:
            if line == 'END_DNA\n':
                break
            if line == 'BEGIN_DNA\n':
                counter += 1
                continue
            if counter == 1:
                base_call.append(line.split(' ')[0])
                stringdex.append(int(line.split(' ')[2].rstrip('\n')))

        # Read ab1 file
        abif = SeqIO.read("raw_data\\sequence_files\\cleaned_" + sample +
                          "-ITS1F.ab1", "abi")

        # Create pandas dataframe of relevant channels
        channel_set = {'Chnl9': abif.annotations["abif_raw"]["DATA9"],
                       'Chnl10': abif.annotations["abif_raw"]["DATA10"],
                       'Chnl11': abif.annotations["abif_raw"]["DATA11"],
                       'Chnl12': abif.annotations["abif_raw"]["DATA12"]}
        trace_df = pd.DataFrame(channel_set,
                                columns=['Chnl9', 'Chnl10',
                                         'Chnl11', 'Chnl12'])

        # Subset trace dataframe to basecall indices and reset index
        trace_df = trace_df.iloc[stringdex]
        trace_df = trace_df.reset_index(drop=True)

        # Find channel with highest signal per base and value of signal
        max_chnl = trace_df.idxmax(axis=1)
        max_signal = trace_df.max(axis=1)

        # Add to data frame
        trace_df['max_chnl'] = max_chnl
        trace_df['max_signal'] = max_signal

        # Add basecalls to trace_df
        trace_df['Basecall'] = base_call

        # Make new base calls
        trace_df['new_call'] = trace_df.apply(lambda row: new_call(row),
                                              axis=1)

        # Trim first 20 bases and reset index
        trace_df = trace_df.drop(range(20)).reset_index(drop=True)

        # Create moving sum of max signal values to find where quality
        # deteriorates
        window = 5
        move_sum = [sum(trace_df.iloc[i:i+window, 5]) for
                    i in range(trace_df.shape[0] - (window - 1))]
        trace_df['move_sum'] = move_sum + ([0] * (window - 1))

        # Trim end of sequence where moving sum reaches threshold
        seq_end = trace_df.query('move_sum < 1000').index.tolist()[0]
        trace_df = trace_df.iloc[range(seq_end), ]

        # Sum number of discrepancies between new and original calls
        # sum(trace_df['Basecall'] != trace_df['new_call'])

        # View discrepancies
        # diffs = np.where(trace_df['Basecall'] != trace_df['new_call'])
        # diffs_df = trace_df.iloc[diffs]

        # Write sequence to file in fasta format
        new_seq = ''.join(list(trace_df['new_call']))
        test_seq = SeqRecord(Seq(new_seq), id="SampleID" + sample,
                             description="No description")
        SeqIO.write(test_seq, "annotated_sequences\\sample" + sample +
                    "_antd.faa", "fasta")
