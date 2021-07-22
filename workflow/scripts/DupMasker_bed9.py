#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Author: Mitchell R. Vollger
import argparse
import pandas as pd


def hex_to_rgb(h):
    h = h.lstrip("#")
    return ",".join(tuple(str(int(h[i : i + 2], 16)) for i in (0, 2, 4)))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="", formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("infile", help="input fasta file")
    parser.add_argument("outfile", help="input fasta file")
    args = parser.parse_args()

    color = pd.read_csv(args.infile, sep="\s+")
    # these rows are no good, they come from contigs that have messed up results. fix TODO
    bad = (color["chr"] == "0") & (color["chrEnd"] == "#BEBEBE")
    color.drop(color[bad].index, inplace=True)

    color["strand"] = "+"
    color.loc[color["orient"] == "R", "strand"] = "-"

    color["rgb"] = color["color"].map(hex_to_rgb)
    color["score"] = 0

    out = color[
        [
            "chr",
            "chrStart",
            "chrEnd",
            "Repeat",
            "score",
            "strand",
            "chrStart",
            "chrEnd",
            "rgb",
        ]
    ]

    out = out.rename(columns={"chr": "#chr"})

    out.to_csv(args.outfile, sep="\t", index=False, compression='gzip')
