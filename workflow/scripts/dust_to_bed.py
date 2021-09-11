#!/usr/bin/env python3
import sys
import argparse

USAGE = """
convert-dustout-to-bed.py  -i <sddust file> -o <output file>

convert output of SDDust to bed format
"""
###############################################################################
parser = argparse.ArgumentParser(
    description=USAGE, formatter_class=argparse.ArgumentDefaultsHelpFormatter
)
parser.add_argument("-i", "--sddustFile", default=snakemake.input.intervals)
parser.add_argument("-o", "--outbed", default=snakemake.output.bed)
args = parser.parse_args()

outfile = open(args.outbed, "w+")
for line in args.sddustFile.open():
    line = line.rstrip()
    if line[0] == ">":
        name = line[1:]
    else:
        line = line.split()
        b = int(line[0])
        e = int(line[2])
        outfile.write("%s\t%i\t%i\n" % (name, b, e + 1))

outfile.close()
