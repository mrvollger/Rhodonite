#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Author: Mitchell R. Vollger
import argparse
import pysam
import textwrap

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="", formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("infile", help="input fasta file")
    parser.add_argument(
        "--outputs", nargs="+", help="list of output files", required=True
    )
    parser.add_argument(
        "-w", "--width", help="width of fasta output", type=int, default=60
    )
    args = parser.parse_args()
    N_IDS = len(args.outputs)

    fasta = pysam.FastaFile(args.infile)

    outs = [open(f, "w+") for f in args.outputs]
    out_idx = 0
    for name in fasta.references:
        seq = fasta.fetch(name)
        seq_fold = "\n".join(textwrap.wrap(seq, args.width)).strip()
        outs[out_idx].write(">{}\n{}\n".format(name, seq_fold))
        out_idx += 1
        if out_idx == N_IDS:
            out_idx = 0

    for out in outs:
        out.close()
