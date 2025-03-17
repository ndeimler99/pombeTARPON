#!/usr/bin/env Rscript

library(ggplot2)
library(dplyr)


# script funtions to generate summary stats of entire sequencing run/experiment
args = commandArgs(trailingOnly=TRUE)
data_table <- read.table(args[0], sep="\t", header=TRUE)
data_table <- reshape2::melt(data_table)
data_table <- data_table[data_table$variable=="Perc_Nucl",]
aln_plot <- ggplot(data=data_table) +
  geom_bar(mapping=aes(x=Sample,y=value,fill=Chrom), stat="identity") +
  scale_fill_manual(breaks=c("I", "II", "III", "MT"), values=c("#D81B60", "#1E88E5", "#FFC107", "#004D40")) +
  theme_minimal() +
  ylab("Percentage of Reads") +
  theme(axis.title.x = element_blank(),
        axis.title.y=element_text(size=20),
        axis.legend.title = element_text(size=15),
        axis.text=element_text(size=15))

ggsave("alignment_distribution.pdf", device="pdf", width=4, height=12, plot=aln_plot)