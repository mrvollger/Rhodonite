#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Author: Mitchell R. Vollger
import argparse
import pysam
import textwrap

print(snakemake.input)
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="", formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        "--infile", help="input fasta file", default=snakemake.input.fasta
    )
    parser.add_argument(
        "--outputs",
        nargs="+",
        help="list of output files",
        default=snakemake.output.fasta,
    )
    parser.add_argument(
        "-w", "--width", help="width of fasta output", type=int, default=60
    )
    parser.add_argument(
        "-r",
        "--reheader",
        help="Make new fasta names and save them in this file",
        default=None,
    )
    parser.add_argument(
        "-m",
        "--maxheader",
        help="Set a maximum header length for the input fasta",
        default=40,
        type=int,
    )

    args = parser.parse_args()
    N_IDS = len(args.outputs)

    if args.reheader is not None:
        conversion = open(args.reheader, "w+")
    fasta = pysam.FastaFile(args.infile)

    outs = [open(f, "w+") for f in args.outputs]
    out_idx = 0
    for idx, name in enumerate(fasta.references):
        seq = fasta.fetch(name)

        if args.reheader is not None:
            new_name = f"RM_NAME_{idx}"
            conversion.write(f"{new_name}\t{name}\n")
            name = new_name

        if args.maxheader is not None and len(name) > args.maxheader:
            raise IOError(
                f"fasta header {name} is longer than {args.maxheader}. Please change your headers to be shorter than {args.maxheader} characters. (https://github.com/rmhubley/RepeatMasker/issues/59)"
            )

        seq_fold = "\n".join(textwrap.wrap(seq, args.width)).strip()
        outs[out_idx].write(">{}\n{}\n".format(name, seq_fold))
        out_idx += 1
        if out_idx == N_IDS:
            out_idx = 0

    for out in outs:
        out.close()
