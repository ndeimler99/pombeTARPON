#!/usr/bin/env python3 

import argparse
import pysam

def create_fasta_dict(fh):
    fasta_dict = {}
    with open(fh, 'r') as fasta:
        header = ''
        seq = ''
        stats = ''
        for line in fasta:
            if line.startswith(">"):
                if seq != '':
                    fasta_dict[header] = seq
                seq = '' 
                header = line.split()[0].strip()
            else:
                seq += line.strip()
        fasta_dict[header] = seq
    return fasta_dict


def get_terminal_lengths(seq_dict):
    terminal_pos = []
    for read in seq_dict:
        terminal_pos.append(len(seq_dict[read].rstrip("-"))-1)
    return terminal_pos

def get_count(i, terminal_pos):
#     returns the number of reads with length greater or equal to i
    return (sum([1 for read in terminal_pos if read >= i]))

def get_max_nucl(nucl_dict, read_count):

    nucl = max(nucl_dict, key=nucl_dict.get)
    nucl_value = max(nucl_dict.values())
    
    if nucl_value > read_count:
        max_val = 0
        for nuc in nucl_dict:
            if nucl_dict[nuc] > max_val and nuc != '-':
                max_val = nucl_dict[nuc]
                nucl = nuc
                nucl_value = nucl_dict[nuc]
    return nucl, nucl_value


def get_consensus_seq(fasta_dict, min_perc):
    terminal_pos = get_terminal_lengths(fasta_dict)
    consensus_seq = []
    coverage = []
    consensus_cov = []

    for i in range(0, max(terminal_pos)):
        nucl_dict = {'A':0, 'T':0, 'C':0, 'G':0, '-':0}
        read_count = get_count(i, terminal_pos)
        for read in fasta_dict:
            nucl = fasta_dict[read][i].upper()
            nucl_dict[nucl] += 1
        nucl, nucl_value = get_max_nucl(nucl_dict, read_count)
        if nucl_value / read_count * 100 >= min_perc and nucl != '-':
            consensus_seq.append(nucl)
            coverage.append(read_count)
            consensus_cov.append(nucl_value)

    return consensus_seq, coverage, consensus_cov

def main(args):

    args.seq_percentage = int(args.seq_percentage)
    with open(args.output_fasta, "w") as out_fasta, open(args.output_stats, "w") as out_stats:
        out_stats.write("Consensus\tPos\tRead_Count\tConsensus_Cov\n")
        for file in args.consensus_files:
            label = file.strip().split(".")[1]
            seq_dict = create_fasta_dict(file)
            consensus_seq, coverage, consensus_cov = get_consensus_seq(seq_dict, args.seq_percentage)
            out_fasta.write(">{}\n{}\n".format(label, ''.join(consensus_seq)))
            for i in range(0,len(coverage)):
                out_stats.write("{}\t{}\t{}\t{}\n".format(label, i, coverage[i], consensus_cov[i]))

def argparser():
    parser = argparse.ArgumentParser()
    parser.add_argument("--consensus_files", nargs="+", required=True)
    parser.add_argument("--output_fasta", required=True)
    parser.add_argument("--output_stats", required=True)
    parser.add_argument("--seq_percentage", required=True)
    return parser

if __name__ == "__main__":
    args = argparser().parse_args()
    main(args)