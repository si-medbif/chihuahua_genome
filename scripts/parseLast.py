#!/usr/bin/env python3
"""
Calculates id% and cov% for all last hits.
Adds columns for end coordinate and for positive strand coordinates for Query.
"""

__author__ = "Harald Grove"
__version__ = "0.1.0"
__license__ = "MIT"

import argparse
import time
import sys
import pandas as pd


class Last(object):

    def __init__(self, lastfile):
        self.lastfile = lastfile
        self.lines = []
        self.header = []

    def _calc_seqid(self, score, blocks):
        """
        Calculates the sequence identity from the alignment
        """
        alnSize = 0
        gaps = 0
        gapSize = 0
        for b in blocks.split(","):
            if ":" in b:
                gap = max([int(a) for a in b.split(":")])
                gaps += gap
                gapSize += 21 + 9 * gap
            else:
                alnSize += int(b)
        mismatch = ((alnSize * 6 - gapSize) - score) / 24
        return 1 - (mismatch + gaps) / (gaps + alnSize)

    def read_last(self, idlim=0, covlim=0, lenlim=0):
        """
        Reads each alignment from the last-alignment file.
        Calculates and adds sequence identity and coverage.
        :param idlim: Minimum identity, in percent
        :param covlim: Minimum coverage, in percent
        :param lenlim: Minimum alignment length
        """
        head = True
        with open(self.lastfile, "r") as fin:
            for line in fin:
                # Add all comment lines at start of file to 'header'
                if line.startswith("#"):
                    if head:
                        self.header.append(line)
                    continue
                head = False
                l = line.strip().split()
                if len(l) != 14:
                    continue
                score, start1, alnSize1, seqSize1, start2, alnSize2, seqSize2 = [
                    int(i) for i in [l[0], l[2], l[3], l[5], l[7], l[8], l[10]]
                ]
                name1, strand1, name2, strand2, blocks, *e = [
                    l[1], l[4], l[6], l[9], l[11], l[12:]
                ]
                end1 = start1 + alnSize1
                end2 = start2 + alnSize2
                if strand2 == "-":
                    start2p = seqSize2 - end2
                    end2p = seqSize2 - start2
                else:
                    start2p = start2
                    end2p = end2
                seqid = 100 * self._calc_seqid(score, blocks)
                if seqSize2 <= seqSize1:
                    seqcov = 100 * (alnSize2 / seqSize2)
                else:
                    seqcov = 100 * (alnSize1 / seqSize1)
                # Deciding whether to keep the alignment or not
                if seqid < idlim:
                    continue
                if alnSize2 < lenlim:
                    continue
                if seqcov < covlim:
                    continue
                self.lines.append(
                    [ score, seqid, seqcov, name1, start1, alnSize1, end1, strand1, seqSize1,
                        name2, start2, alnSize2, end2, strand2, seqSize2, blocks,
                        start2p, end2p ]
                )

    def write_last(self, fout=sys.stdout):
        for line in self.header:
            l = line.strip().split()
            if len(l) > 1 and l[1] == "score":
                l.insert(2, "idpct")
                l.insert(3, "covpct")
                l.insert(7, "end1")
                l.insert(13, "end2")
                l.append("start2+")
                l.append("end2+")
            else:
                continue
            line = "{}\n".format("\t".join(l[1:]))
            fout.write(line)
        for line in self.lines:
            output = "\t".join([str(b) for b in line])
            fout.write("{}\n".format(output))


def main(args):
    last = Last(args.infile)
    idlim, covlim, lenlim = [int(a) for a in args.limits.split(",")]
    last.read_last(idlim=idlim, covlim=covlim, lenlim=lenlim)
    if args.outfile is not None:
        last.write_last(open(args.outfile, "w"))
    else:
        last.write_last()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    # Required positional argument
    parser.add_argument("infile", help="Last alignment file")

    # Optional argument which requires a parameter (eg. -d test)
    parser.add_argument("-o", "--outfile", help="Output, parsed last alignment.")
    parser.add_argument(
        "-c", "--limits", default="0,0,0", help="Minimum idpct, covpct and length"
    )
    # parser.add_argument("-n", "--name", action="store", dest="name")

    # Optional verbosity counter (eg. -v, -vv, -vvv, etc.)
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
