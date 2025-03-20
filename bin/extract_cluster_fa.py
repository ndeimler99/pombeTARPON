#!/usr/bin/env python3 

import argparse

def rev_complement(seq):
    nucl_dict = {"A":"T", "T":"A","C":"G", "G":"C"}
    return ''.join([nucl_dict[i] for i seq[::-1]])
    
def meshclust_clusters(mesh_file):
    cluster_dict = {}
    with open(mesh_file, 'r') as fh:
        for line in fh:
            line = line.strip().split()
            if line == []:
                continue
            if line[0] in cluster_dict:
                cluster_dict[line[0]].append(line[1].strip())
            else:
                cluster_dict[line[0]] = [line[1].strip()]
    return cluster_dict

def get_fasta_dict(fasta_file):
    fasta_dict = {}
    with open(fasta_file, "r") as fasta_fh:
        header = ""
        seq = ""
        for line in fasta_fh:
            if line.startswith(">"):
                if header != "":
                    fasta_dict[header] = seq
                header = line.strip()
                seq = ""
            else:
                seq += line.strip()
        fasta_dict[header] = seq
    return fasta_dict

def main(args):

    args.minimum_cluster_size = int(args.minimum_cluster_size)

    seq_dict = get_fasta_dict(args.input_file)
    cluster_dict = meshclust_clusters(args.cluster_file)

    for cluster in cluster_dict:
        if len(cluster_dict[cluster]) > args.minimum_cluster_size:
            with open("{}.cluster_{}.fa".format(args.prefix, cluster), "w") as cluster_fh:
                for read in cluster_dict[cluster]:
                    cluster_fh.write("{}\n{}\n".format(read, rev_complement(seq_dict[read])))

def argparser():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_file", required=True)
    parser.add_argument("--cluster_file", required=True)
    parser.add_argument("--prefix", required=True)
    parser.add_argument("--minimum_cluster_size", required=True)
    return parser

if __name__ == "__main__":
    args = argparser().parse_args()
    main(args)