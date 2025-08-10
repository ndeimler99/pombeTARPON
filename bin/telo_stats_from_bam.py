#!/usr/bin/env python3 

import argparse
import pysam

def avg_quality(seq):
    return sum([i for i in seq])/len(seq)

def main(args):

    aln_file = pysam.AlignmentFile(args.input_file, "rb", check_sq=False)

    with open(args.stats_file, "w") as stats_fh:
        stats_fh.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format("read_id", "strand", "read_length", "read_quality", "telo_start", "telo_end", "telo_length", "telo_quality"))
        for aln in aln_file:
            # print(aln.query_name)
            # print(aln.get_tag("XT"))
            # print(aln.get_tag("XE"))
            # print(len(aln.query_sequence))
            # print(aln.query_sequence[aln.get_tag("XT"):aln.get_tag("XE")])
            # print(aln.query_qualities[aln.get_tag("XT"):aln.get_tag("XE")])
            stats_fh.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(aln.query_name, aln.get_tag("XS"), len(aln.query_sequence), 
                                                                        avg_quality(aln.query_qualities), aln.get_tag("XT"), aln.get_tag("XE"),
                                                                        aln.get_tag("XL"), avg_quality(aln.query_qualities[aln.get_tag("XT"):aln.get_tag("XE")])))


def argparser():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_file", required=True)
    parser.add_argument("--stats_file", required=True)
    return parser

if __name__ == "__main__":
    args = argparser().parse_args()
    main(args)