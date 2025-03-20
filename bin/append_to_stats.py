#!/usr/bin/env python3 

import argparse

def meshclust_clusters(mesh_file,min_size):
    cluster_dict = {}
    with open(mesh_file, 'r') as fh:
        for line in fh:
            line = line.strip().split()
            if line == []:
                continue
            if line[0] in cluster_dict:
                cluster_dict[line[0]].append(line[1].strip(">"))
            else:
                cluster_dict[line[0]] = [line[1].strip(">")]

    read_dict = {}
    for cluster in cluster_dict:
        if len(cluster_dict[cluster]) >= min_size:
            for read in cluster_dict[cluster]:
                read_dict[read] = cluster

    return read_dict

def get_fasta_dict(fasta_file):
    return_dict = {}
    with open(fasta_file, "r") as fh:
        header =''
        seq = ''
        for line in fh:
            if line.startswith(">"):
                if header != "":
                    return_dict[header.strip(">")] = seq
                header = line.strip().split()[0]
                seq = ""
            else:
                seq += line.strip()
        return_dict[header.strip(">")] = seq
    return return_dict

# def get_seq_composition(telomere):
#     long = regex.finditer("GGTTACAC{e<=1}")

def main(args):
    
    args.minimum_cluster_size = int(args.minimum_cluster_size)
    cluster_dict = meshclust_clusters(args.cluster_file, args.minimum_cluster_size)
    fasta_dict = get_fasta_dict(args.fasta_file)

    with open(args.stats_file, "r") as stats_fh, open(args.new_stats_file, "w") as stats_out:
        linecount = 0
        for line in stats_fh:
            line = line.strip().split()
            if linecount == 0:
                linecount += 1
                line.extend(["Cluster"])
                # line.extend(["Cluster","GGTTAC","GGTTACA","GGTTACAC"])
                stats_out.write("{}\n".format('\t'.join(line)))
                continue
            
            if line[0] in cluster_dict:
                line.append(cluster_dict[line[0]])
            else:
                line.append("Not_Clustered")
            
            # freqs = get_seq_composition(fasta_dict[read][int(line[4]):])
            stats_out.write("{}\n".format('\t'.join(line)))
            

def argparser():
    parser = argparse.ArgumentParser()
    parser.add_argument("--stats_file", required=True)
    parser.add_argument("--cluster_file", required=True)
    parser.add_argument("--minimum_cluster_size", required=True)
    parser.add_argument("--new_stats_file", required=True)
    parser.add_argument("--fasta_file", required=True)
    parser.add_argument("--repeat", required=True)
    return parser

if __name__ == "__main__":
    args = argparser().parse_args()
    main(args)