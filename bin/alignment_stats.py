#!/usr/bin/env python3 

import argparse
import pysam



def main(args):

    aln_file = pysam.AlignmentFile(args.alignment_file, "rb", check_sq=False)

    hypothetical = {}
    for chrom in aln_file.header.to_dict()['SQ']:
        if chrom["SN"] in ["AB325691", "MTR"]:
            continue
        hypothetical[chrom["SN"]] = chrom["LN"]

    alignment_counts = {}
    for aln in aln_file:
        if aln.query_sequence is None or aln.reference_name is None:
            continue
        if aln.reference_name in ["AB325691", "MTR"]:
            continue
        if aln.reference_name in alignment_counts:
            alignment_counts[aln.reference_name]["count"] += 1
            alignment_counts[aln.reference_name]["nucl_count"] += len(aln.query_sequence)
        else:
            alignment_counts[aln.reference_name] = {"count": 1, "nucl_count": len(aln.query_sequence)}
    
    alignment_results = {}
    for chrom in alignment_counts:
        alignment_results[chrom] = {"perc_reads": alignment_counts[chrom]["count"] / sum([alignment_counts[chrom]["count"] for chrom in alignment_counts]) * 100,
                                    "perc_nucl": alignment_counts[chrom]["nucl_count"] / sum([alignment_counts[chrom]["nucl_count"] for chrom in alignment_counts])*100}

    
    with open(args.stats_file, "w") as stats_fh:
        stats_fh.write("Sample\tChrom\tPerc_Reads\tPerc_Nucl\n")
        for chrom in alignment_results:
            stats_fh.write("{}\t{}\t{}\t{}\n".format(args.sample, chrom, alignment_results[chrom]["perc_reads"], alignment_results[chrom]["perc_nucl"]))

   

def argparser():
    parser = argparse.ArgumentParser()
    parser.add_argument("--alignment_file", required=True)
    parser.add_argument("--stats_file", required=True)
    parser.add_argument("--sample", required=True)
    return parser

if __name__ == "__main__":
    args = argparser().parse_args()
    main(args)