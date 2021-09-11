#!/usr/bin/env python
import argparse
import os
import sys
import pysam
from numba import njit
import multiprocessing
from functools import partial


@njit
def get_gaps(name, seq, soft=False):
    counter = 0
    idx = 0
    out = []
    for char in seq:
        if char.upper() == "N":
            counter += 1
        elif soft and char.islower():
            counter += 1
        else:
            if counter > 0:
                out.append((idx - counter, idx))
            counter = 0
        idx += 1
    return out


def fetch_seq(name, fasta, soft):
    fasta = pysam.FastaFile(fasta)
    seq = fasta.fetch(name)
    intervals = get_gaps(name, seq, soft=soft)
    return (name, intervals)


# global var for inputs
args = None

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="", formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        "infile", help="positional input", default=snakemake.input.fasta
    )
    parser.add_argument(
        "outfile", help="positional output bed", default=snakemake.output.bed
    )
    parser.add_argument(
        "-s",
        "--soft",
        help="in addition to Ns count soft masked sequence as well",
        action="store_true",
        default=False,
    )
    parser.add_argument("-t", "--threads", help="n threads to use", type=int, default=8)
    parser.add_argument(
        "-d", help="store args.d as true if -d", action="store_true", default=False
    )
    args = parser.parse_args()

    fasta = pysam.FastaFile(args.infile)
    refs = fasta.references
    fasta.close()
    # sys.stderr.write(f"{refs}\n")
    out_bed = open(args.outfile, "w+")
    with multiprocessing.Pool(args.threads) as pool:
        fetch_seq_extra = partial(fetch_seq, fasta=args.infile, soft=args.soft)
        for i, rtn in enumerate(pool.imap(fetch_seq_extra, refs)):
            contig, intervals = rtn
            for start, end in intervals:
                out_bed.write(f"{contig}\t{start}\t{end}\n")
            sys.stderr.write(f"{contig} done\n")
    out_bed.close()
    # for NotADirectoryError in fasta:
    # sys.stderr.write(f"{rec.name}\n")
    old = """
        counter=0
        idx = 0
        for char in rec.sequence:
            if(char.upper() == "N"):
                counter +=1
            else:
                if( counter > 0 ):
                    print("{}\t{}\t{}".format(rec.name, idx-counter, idx))
                counter = 0
            idx += 1
        """
