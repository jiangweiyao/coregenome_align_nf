#!/usr/bin/env python

import sys
import os
import argparse
import pandas as pd
from tabulate import tabulate

def main():

    cli = argparse.ArgumentParser()
    cli.add_argument('-i', '--InputFile', help="Input file is Mash Screen Output File", required=True)
    cli.add_argument('-o', '--OutputFile', help=f"Output file for parsed result", required=True)
    cli.add_argument('-n', '--number', help="Display top N matches", type=int, required=False, default=5)
    cli.add_argument('-m', '--message', help="Display Message", required=False, default="Mash Screen Results:")

    args = cli.parse_args()
    
    InputFile = os.path.expanduser(args.InputFile)
    print(InputFile)
    OutputFile = os.path.expanduser(args.OutputFile)
    print(OutputFile)
    table1 = pd.read_csv(InputFile,sep="\t", header = None)
    table1.columns = ['identity', 'shared_hashes', 'multiplicity', 'p_value', 'query_ID', 'ID']

    table1['Reference'] = table1['ID'].str.slice(0,80)
    table2 = table1[~table1.ID.str.contains("phage")]
    table1 = table2[["Reference", "shared_hashes", "identity", "multiplicity"]]
    #print(tabulate(table1.head(10), floatfmt=".4f", headers="keys", showindex=False, tablefmt="psql"))
    if os.path.exists(OutputFile):
        append_write = 'a' # append if already exists
    else:
        append_write = 'w' # make a new file if not
    dff = open(OutputFile, append_write) 
    dff.write(f"{args.message}\n")
    dff.write(tabulate(table1.head(args.number), floatfmt=".4f", headers="keys", showindex=False, tablefmt="psql")) 
    dff.write("\n\n")
    dff.close()

if __name__ == "__main__":
    sys.exit(main())
