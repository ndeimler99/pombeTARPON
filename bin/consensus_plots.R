#!/usr/bin/env Rscript

library(ggplot2)
library(dplyr)


# script funtions to generate summary stats of entire sequencing run/experiment
args = commandArgs(trailingOnly=TRUE)
data_table <- read.table(args[1], sep="\t", header=TRUE)
for (cluster in unique(data_table$Consensus)){
  tmp_df <- data_table[data_table$Consensus == cluster,]
  plt <- ggplot(data=tmp_df) +
    geom_line(mapping=aes(x=Pos, y=Read_Count), size=2, color="#D81B60", label="All Reads") +
    geom_line(mapping=aes(x=Pos, y=Consensus_Cov), size=2, color="#1E88E5", label="Reads into Consensus") +
    xlab("Distance from Telomere") + ylab("Number of Reads") +
    theme_minimal() +
    theme(axis.title=element_text(size=20),
          axis.text=element_text(size=15),
          legend.text = element_text(size=15))
  ggsave(paste(cluster, "_coverage_plot.pdf", sep=""), device="pdf", width=12, height=8, plot=plt)
}
# data_table <- data_table[data_table$variable=="Perc_Nucl",]
# aln_plot <- ggplot(data=data_table) +
#   geom_bar(mapping=aes(x=Sample,y=value,fill=Chrom), stat="identity") +
#   scale_fill_manual(breaks=c("I", "II", "III", "MT"), values=c("#D81B60", "#1E88E5", "#FFC107", "#004D40")) +
#   theme_minimal() +
#   ylab("Percentage of Reads") +
#   theme(axis.title.x = element_blank(),
#         axis.title.y=element_text(size=20),
#         legend.title = element_text(size=15),
#         axis.text=element_text(size=15))

# ggsave("alignment_distribution.pdf", device="pdf", width=4, height=12, plot=aln_plot)