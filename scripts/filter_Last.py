#!/usr/bin/env python3
"""
Filter a parsed Last report file, attempting to remove invalid lines.
Input should prepared using the 'parseLast.py' script.
"""

__author__ = "Harald Grove"
__version__ = "0.1.0"
__license__ = "MIT"

import argparse
import time
import sys
import pandas as pd

def process_queries(args):
    """
    Collects all alignments from a single query sequence and keep the highest scoring
    alignment for each set of overlapping alignments.
    Each new segment is compared to the last overlapping segment that was approved.
    """
    idlim, covlim, lenlim = [int(a) for a in args.limits.split(",")]
    df = pd.read_table(args.infile, header=0)
    with open(args.transcriptfile, 'r') as fin:
        for line in fin:
            transcript = line.strip()
            df1 = df[df['name2']==transcript].sort_values(by='start2+')
            if len(df1) == 0:
                continue
            df1.reset_index(drop=True, inplace=True)
            if len(df1) > 1:
                df1['Keep'] = pd.Series(len(df1))
                index2, index3 = 0,1
                while index3 < len(df1):
                    score2, score3 = df1.loc[index2,'score'], df1.loc[index3,'score']
                    start3 = df1.loc[index3,'start2+']
                    end2 = df1.loc[index2,'end2+']
                    # Next segment is not overlapping (based on maximum allowable overlap).
                    if end2-start3 < args.overlap: 
                        df1.loc[index2, 'Keep'] = True
                        df1.loc[index3, 'Keep'] = True
                        index2 = index3
                    # Pick the highest scoring segment.
                    else:
                        if score2 > score3:
                            df1.loc[index2, 'Keep'] = True
                            df1.loc[index3, 'Keep'] = False
                        else:
                            df1.loc[index2, 'Keep'] = False
                            df1.loc[index3, 'Keep'] = True
                            index2 = index3
                    index3 += 1
                df2 = df1[df1['Keep']]
                df2 = df2.drop(columns=['Keep'])
            else:
                df2 = df1
            outfile = f'{args.out}.{transcript}.tsv'
            df2.to_csv(outfile, sep='\t', index=None)


def main(args):
    last = Last(args.infile)
    process_queries(args)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    # Required positional argument
    parser.add_argument("infile", help="Last alignment file")
    parser.add_argument("transcriptfile", help="List of query sequences to process.")

    # Optional argument which requires a parameter (eg. -d test)
    parser.add_argument(
        "-o", "--out", help="Prefix for outputfile.", default = ""
    )
    parser.add_argument(
        "-l", "--limits", default="0,0,0", help="Minimum idpct, covpct and length"
    )
    parser.add_argument(
        "--overlap", default=7, help="Max overlap between neighboring query segments."
    )

    # Optional verbosity counter
    parser.add_argument(
        "-v", "--verbose", action="count", default=0, help="Verbosity (-v, -vv, etc)"
    )

    # Specify output of '--version'
    parser.add_argument(
        "--version",
        action="version",
        version="%(prog)s (version {version})".format(version=__version__),
    )

    args = parser.parse_args()
    main(args)
