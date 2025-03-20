#!/usr/bin/env Rscript

library(ggplot2)
library(dplyr)


# script funtions to generate summary stats of entire sequencing run/experiment
args = commandArgs(trailingOnly=TRUE)
data_table <- read.table(args[1], sep="\t", header=TRUE)
