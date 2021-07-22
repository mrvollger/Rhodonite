#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Author: Mitchell R. Vollger
import argparse
import sys
import pandas as pd

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="", formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        "infiles", nargs="+", help="list of input .dat file(s) from trf"
    )
    parser.add_argument("-o", "--outfile", help="output bed file", default=sys.stdout)
    args = parser.parse_args()

    trf = []
    header = "#chr start end PeriodSize CopyNumber ConsensusSize PercentMatches PercentIndels Score A C G T Entropy Motif Sequence".split()

    for dat_f in args.infiles:
        chrom = None
        sys.stderr.write("\r" + dat_f)
        with open(dat_f, "r") as dat:
            for line in dat:
                splitline = line.split()
                if line.startswith("Sequence:"):
                    chrom = int(line.split()[1].strip())
                    # sys.stderr.write(chrom + "\n")
                elif line.startswith("@"):
                    chrom = splitline[0][
                        1:
                    ].strip()  # grab everything after the @ in the first word
                else:
                    # Catch index errors when line is blank
                    try:
                        # Check if in header sequence (all non-header lines start with an int: start pos)
                        try:
                            int(splitline[0])
                        except ValueError:
                            continue
                        trf.append([chrom] + splitline[0 : (len(header) - 1)])
                    except IndexError:
                        pass

    sys.stderr.write("\n")
    trf = pd.DataFrame(trf, columns=header)
    trf.to_csv(args.outfile, sep="\t", index=False)
