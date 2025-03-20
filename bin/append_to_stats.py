#!/usr/bin/env python3 

import argparse
import cairo
import numpy as np

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


def main(args):
    
    args.minimum_cluster_size = int(args.minimum_cluster_size)
    cluster_dict = meshclust_clusters(args.cluster_file, args.minimum_cluster_size)
    with open(args.stats_file, "r") as stats_fh, open(args.new_stats_file, "w") as stats_out:
        linecount = 0
        for line in stats_fh:
            line = line.strip().split()
            if linecount == 0:
                linecount += 1
                line.append("Cluster")
            elif line[0] in cluster_dict:
                line.append(cluster_dict[line[0]])
            else:
                line.append("Not_Clustered")
            stats_out.write("{}\n".format('\t'.join(line)))
            

def argparser():
    parser = argparse.ArgumentParser()
    parser.add_argument("--stats_file", required=True)
    parser.add_argument("--cluster_file", required=True)
    parser.add_argument("--minimum_cluster_size", required=True)
    parser.add_argument("--new_stats_file", required=True)
    return parser

if __name__ == "__main__":
    args = argparser().parse_args()
    main(args)