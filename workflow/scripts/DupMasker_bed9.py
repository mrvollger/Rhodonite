#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Author: Mitchell R. Vollger
import argparse
import pandas as pd
import sys


def hex_to_rgb(h):
    h = h.lstrip("#")
    return ",".join(tuple(str(int(h[i : i + 2], 16)) for i in (0, 2, 4)))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="", formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("infile", help="input fasta file")
    parser.add_argument("--outfile", help="input fasta file", default=sys.stdout)
    args = parser.parse_args()

    # chr     chrStart        chrEnd  orient  Repeat  color   width   offset
    d_types = {
        "chr": str,
        "chrStart": "Int64",
        "chrEnd": "Int64",
        "offset": "Int64",
        "width": "Int64",
        "color": str,
    }
    color = pd.read_csv(args.infile, sep="\t", dtype=d_types)
    # these rows are no good, they come from contigs that have messed up results. fix TODO
    bad = (color["chr"] == "0") & (color["chrEnd"] == "#BEBEBE")
    color.drop(color[bad].index, inplace=True)

    color["strand"] = "+"
    color.loc[color["orient"] == "R", "strand"] = "-"

    color["rgb"] = color["color"].map(hex_to_rgb)
    color["score"] = 0
    color["tst"] = color["chrStart"]
    color["ten"] = color["chrEnd"]

    out = color[
        [
            "chr",
            "chrStart",
            "chrEnd",
            "Repeat",
            "score",
            "strand",
            "tst",
            "ten",
            "rgb",
        ]
    ]
    ##ct     st      en      name    score   strand  tst     ten     color
    out = out.rename(
        columns={
            "chr": "#ct",
            "chrStart": "st",
            "chrEnd": "en",
            "Repeat": "name",
            "rgb": "color",
        }
    )

    out.to_csv(args.outfile, sep="\t", index=False)
