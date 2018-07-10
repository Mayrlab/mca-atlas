#!/usr/bin/env python

import pandas as pd
from multiprocessing import Pool
from functools import partial
import argparse


# FUNCTIONS
def applyParallel(grouped, func, cores):
    # Credit: https://stackoverflow.com/a/29281494/570918
    with Pool(cores) as p:
        results = p.map(func, [group for name, group in grouped])
    return pd.concat(results)


def merge_to_mode(df, epsilon):
    idx = 0
    while idx < len(df):
        cur = df.iloc[idx]
        neighborhood = df.loc[(df.position >= (cur.position - epsilon)) & (df.position <= (cur.position + epsilon)), 'reads']
        max_name = neighborhood.idxmax()
        if max_name == cur.name:
            idx += 1
        else:
            df.loc[max_name, 'reads'] = cur.reads + df.loc[max_name, 'reads']
            df.drop(cur.name, inplace=True)

    return df


# PARSE ARGUMENTS

parser = argparse.ArgumentParser(description="Merge neighboring read counts within an interval.")
parser.add_argument('-e', '--epsilon', type=int, default=12,
                    help="Number of nucleotides to search upstream and downstream to merge with.")
parser.add_argument('-p', '--processes', type=int, default=1, help="Number of cores to use.")
parser.add_argument('inFile')
parser.add_argument('outFile')

args = parser.parse_args()


# IMPORT FILE
counts = pd.read_table(args.inFile, header=None, names=['seqname', 'position', 'reads'])
print("Imported %d genomic locations" % len(counts))

merged_counts = applyParallel(counts.groupby('seqname'), partial(merge_to_mode, epsilon=args.epsilon), args.processes)
print("Exporting %d genomic locations" % len(merged_counts))

merged_counts.to_csv(args.outFile, sep='\t', header=False, index=False, compression='gzip')
