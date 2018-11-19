#!/usr/bin/env Rscript

library(readr)
library(dplyr)
library(stringr)
library(ggplot2)

## bypass X11 rendering
## "device='png'" is workaround for bug in ggplot2 v3.0.0
options(bitmapType='cairo', device='png')

args = commandArgs(trailingOnly=TRUE)

if (length(args) != 5) {
    stop("Incorrect number of arguments!\nUsage:\n> plot_peak_cdf.R <covFile> <peakStart> <gene> <batch> <outFile>\n")
}

arg.covFile <- args[1]
arg.peakStart <- as.integer(args[2])
arg.gene <- args[3]
arg.batch <- args[4]
arg.outFile <- args[5]

df <- read_tsv(arg.covFile, col_names=c("chromosome", "position", "count")) %>%
    mutate(distance=abs(position - arg.peakStart)) %>%
    arrange(distance) %>%
    filter(distance < 1000) %>%
    mutate(cdf=cumsum(count))

cutoff <- sum(df$count)*0.95
top95 <- df %>% filter(cdf >= cutoff) %>% summarise(val=first(distance)) %>% as.numeric

g <- ggplot(df, aes(x=distance, y=cdf)) +
    geom_line() + ##geom_hline(yintercept=cutoff) +
    geom_vline(xintercept=top95, color='red', linetype='dashed') +
    labs(title=sprintf("Cumulative Read Coverage for %s in %s", arg.gene, arg.batch),
         x="Distance from Transcript 3' End (nucleotides)",
         y="Cumulative Read Count") +
    geom_text(label="95%", color='red', aes(x=top95, y=0), hjust=-0.5)

ggsave(arg.outFile, g)
