#!/usr/bin/env python3 

import argparse
import pysam
import regex

# def get_ratio(seq, repeat, errors):
#     matches = list(regex.finditer(r'(%s){e<=%s}' % (repeat, errors), seq, overlapped=True))
#     summation = 0
#     for match in matches:
#         summation += match.span()[1] - match.span()[0]
#     return summation

def main(args):

    args.telo_end_repeat_errors = int(args.telo_end_repeat_errors)
    args.telo_end_repeat_count = int(args.telo_end_repeat_count)
    args.telo_start_canonical_errors = int(args.telo_start_canonical_errors)
    args.telo_start_repeat_count = int(args.telo_start_repeat_count)
    args.telo_start_repeat_errors = int(args.telo_start_repeat_errors)

    if args.repeat == "GGTTAC":
        args.repeat = 'GGTTAC((AC)|A){0,1}G{0,7}?'

    aln_file = pysam.AlignmentFile(args.input_file, "rb", check_sq=False)
    coords_file = pysam.AlignmentFile(args.coords_file, "wb", template=aln_file)
    no_coords_file = pysam.AlignmentFile(args.no_coords_file, "wb", template=aln_file)

    telo_seq_freq = []
    for aln in aln_file:
        matches = list(regex.finditer(r'(%s){e<=%s}' % (args.repeat * args.telo_end_repeat_count,args.telo_end_repeat_errors), aln.query_sequence, overlapped=True))
        
        telo_end = None
        telo_start = None
        if len(matches) > 0:
            telo_end = matches[-1].span()[1]
            aln.set_tag("XE", telo_end)

        if telo_end is not None:
            matches = list(regex.finditer(r'(%s){e<=%s}' % (args.repeat * args.telo_start_repeat_count, args.telo_start_repeat_errors), aln.query_sequence, overlapped=True))
            if len(matches) > 0:
                for match in matches:
                    if match.span()[0] > telo_end:
                        continue
                    telo_seq = aln.query_sequence[match.span()[0]:telo_end]
                    freq = telo_seq.count("GGTTAC")*6 / len(telo_seq) * 100
                    if freq >= 50:
                        telo_start = match.span()[0]
                        aln.set_tag("XT", telo_start)
                        break
        
        
            matches = list(regex.finditer(r'(%s){e<=%s}' % (args.telo_start_canonical, args.telo_start_canonical_errors), aln.query_sequence, overlapped=True))
            if len(matches) > 0:
                for match in matches:
                    if match.span()[0] > telo_end:
                        continue
                    if telo_start is None:
                        telo_seq = aln.query_sequence[match.span()[0]:telo_end] 
                        #freq = get_ratio(telo_seq, args.repeat, args.telo_start_repeat_errors) / len(telo_seq) * 100
                        freq = telo_seq.count("GGTTAC")*6 / len(telo_seq) * 100
                        if freq >= 50:
                            telo_start = match.span()[0]
                            aln.set_tag("XT", telo_start)
                            break
                    elif match.span()[0] < telo_start:
                        telo_seq = aln.query_sequence[match.span()[0]:telo_end] 
                        #freq = get_ratio(telo_seq, args.repeat, args.telo_start_repeat_errors) / len(telo_seq) * 100
                        freq = telo_seq.count("GGTTAC")*6 / len(telo_seq) * 100
                        if freq >= 50:
                            telo_start = match.span()[0]
                            aln.set_tag("XT", telo_start)
                            break

        if telo_start is None or telo_end is None:
            no_coords_file.write(aln)
        else:
            aln.set_tag("XL", telo_end-telo_start)
            coords_file.write(aln)

def argparser():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_file", required=True)
    parser.add_argument("--repeat", required=True)
    parser.add_argument("--telo_end_repeat_count")
    parser.add_argument("--telo_end_repeat_errors", required=True)
    parser.add_argument("--coords_file", required=True)
    parser.add_argument("--no_coords_file", required=True)
    parser.add_argument("--telo_start_canonical", required=True)
    parser.add_argument("--telo_start_canonical_errors", required=True)
    parser.add_argument("--telo_start_repeat_count", required=True)
    parser.add_argument("--telo_start_repeat_errors", required=True)
    return parser

if __name__ == "__main__":
    args = argparser().parse_args()
    main(args)