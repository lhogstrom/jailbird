#! /usr/bin/env python
'''
command line tool to run drug-gene anlaysis. Checks the connection 
between a drug and the knockdown of its gene target

input: drug names (BRDs) - gene targets (gene IDs)
cell lines of interest - (default is all)

'''


#three options for running query:
# 1) pre-slice database
# 2) specify subset of sig IDs to query
# 3) run query in full space and ignore extra results