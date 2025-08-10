#!/usr/bin/env python3 

import argparse
import pysam


def rev_complement(seq):
    rev_dict = {'A':'T', 'T':'A', 'C':'G', 'G':'C'}
    return ''.join([rev_dict[i] for i in seq[::-1]])

def main(args):
    
    aln_file = pysam.AlignmentFile(args.input_file, "rb", check_sq=False)
    args.repeat_count = int(args.repeat_count)
    args.fragment_read_length = int(args.fragment_read_length)
    args.min_repeat_threshold = int(args.min_repeat_threshold)
    putative_fh = pysam.AlignmentFile(args.putative_file, "wb", template=aln_file)
    non_telo_fh = pysam.AlignmentFile(args.non_telomeric, "wb", template=aln_file)
    chimeric_fh = pysam.AlignmentFile(args.chimeric_file, "wb", template=aln_file)
    reverse_complemented_fh = pysam.AlignmentFile(args.reverse_complement_file, "wb", template=aln_file)


    for aln in aln_file:
        g_strand = False
        g_count = aln.query_sequence[-(args.fragment_read_length):].count(args.repeat)
        c_count = aln.query_sequence[0:args.fragment_read_length].count(rev_complement(args.repeat))

        if g_count > args.repeat_count or c_count > args.repeat_count:
            if g_count > args.repeat_count:
                aln.set_tag("XS", "G")
                g_strand=True
                putative_fh.write(aln)
            elif c_count > args.repeat_count:
                if not g_strand:
                    aln.set_tag("XS", "C")
                    putative_fh.write(aln)
                else:
                    aln.set_tag("XS", "GC")
                    putative_fh.write(aln)
            
            
            freq = g_count / (g_count + c_count) * 100
            if freq < 100 - args.min_repeat_threshold and freq > args.min_repeat_threshold:
                chimeric_fh.write(aln)
            elif freq <= args.min_repeat_threshold:
                aln.set_tag("XS", "C")
                q=aln.query_qualities
                aln.query_sequence = rev_complement(aln.query_sequence)
                aln.query_qualities = q[::-1]
                reverse_complemented_fh.write(aln)
            elif freq >= 100 - args.min_repeat_threshold:
                aln.set_tag("XS", "G")
                reverse_complemented_fh.write(aln)
            else:
                chimeric_fh.write(aln)

        else:
            non_telo_fh.write(aln)
            


def argparser():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_file", required=True)
    parser.add_argument("--repeat", required=True)
    parser.add_argument("--repeat_count", required=True)
    parser.add_argument("--fragment_read_length", required=True)
    parser.add_argument("--min_repeat_threshold", required=True)
    parser.add_argument("--putative_file", required=True)
    parser.add_argument("--non_telomeric", required=True)
    parser.add_argument("--chimeric_file", required=True)
    parser.add_argument("--reverse_complement_file", required=True)
    return parser

if __name__ == "__main__":
    args = argparser().parse_args()
    main(args)