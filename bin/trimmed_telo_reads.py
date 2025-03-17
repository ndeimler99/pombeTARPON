#!/usr/bin/env python3 

import argparse

def main(args):

    args.read_length = int(args.read_length)
    stats_dict = {}
    with open(args.stats_file, "r") as stats_fh:
        linecount = 0
        for line in stats_fh:
            if linecount == 0:
                linecount += 1
                continue
            line = line.strip().split()
            stats_dict[line[0]] = int(line[4])
    
    with open(args.no_telo, "w") as no_telo_fh, open(args.no_telo_subtelo, "w") as trimmed_fh:
        with open(args.input_file, "r") as input_fh:
            header = ''
            seq = ''
            for line in input_fh:
                if line.startswith(">"):
                    if header != "":
                        no_telo_fh.write("{}\n{}\n".format(header, seq[0:stats_dict[header.strip(">")]]))
                        if len(seq) >= args.read_length:
                            trimmed_fh.write("{}\n{}\n".format(header, seq[stats_dict[header.strip(">")]-args.read_length:stats_dict[header.strip(">")]]))
                    header = line.strip()
                    seq = ''
                else:
                    seq += line.strip()
            no_telo_fh.write("{}\n{}\n".format(header, seq[0:stats_dict[header.strip(">")]]))
            if len(seq) >= args.read_length:
                trimmed_fh.write("{}\n{}\n".format(header, seq[stats_dict[header.strip(">")]-args.read_length:stats_dict[header.strip(">")]]))


def argparser():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_file", required=True)
    parser.add_argument("--stats_file", required=True)
    parser.add_argument("--no_telo", required=True)
    parser.add_argument("--no_telo_subtelo", required=True)
    parser.add_argument("--read_length", required=True)
    return parser

if __name__ == "__main__":
    args = argparser().parse_args()
    main(args)