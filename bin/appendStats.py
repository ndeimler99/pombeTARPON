#!/usr/bin/env python3 

import argparse
import pysam



def main(args):
    first_file = True
    with open(args.outFile, "w") as out_fh:
        for file in args.stats_files:
            label = file.strip().split("/")[-1].split(".")[0]
            with open(file, "r") as fh:
                linecount = 0
                for line in fh:
                    if linecount == 0 and not first_file:
                        linecount += 1
                        continue
                    elif linecount == 0 and first_file:
                        line=line.strip().split("\t")
                        line.append("Sample")
                        out_fh.write('\t'.join(line) + "\n")
                        first_file = False
                        continue
                    line = line.strip().split("\t")
                    line.append(label)
                    out_fh.write('\t'.join(line) + "\n")

def argparser():
    parser = argparse.ArgumentParser()
    parser.add_argument("--stats_files", nargs="+", required=True)
    parser.add_argument("--outFile", required=True)
    return parser

if __name__ == "__main__":
    args = argparser().parse_args()
    main(args)