#!/usr/bin/env Rscript

#' Parses RSeQC's read_distribution.py output into a data.table object
#' then generates a bar plot for where tags mapped in genome.

suppressMessages({
    require(data.table)
    require(ggplot2)
    require(stringr)
})

## bypass X11 rendering 
options(bitmapType = 'cairo')

args = commandArgs(trailingOnly=TRUE)

if (length(args) != 3) {
    stop("Incorrect number of arguments!\nUsage:\n> read_dist_plot.R <label> <inFile> <outFile>\n")
}

arg.label <- args[1]
arg.inFile <- args[2]
arg.outFile <- args[3]

dir.create(dirname(arg.outFile), recursive = TRUE, showWarnings = FALSE)

parseRSeQCLog <- function (fileName) {
    total_tags <- as.numeric(str_extract(readLines(fileName)[5], "\\d+$"))
    dt <- fread(fileName, skip=7, nrows=10)
    dt <- dt[!(Group %in% c("TES_down_1kb", "TSS_up_1kb", "TES_down_10kb", "TSS_up_10kb")),]
    ig_tags <- total_tags - dt[,sum(Tag_count)]
    dt <- rbind(dt, list("Intergenic", NA, ig_tags, NA))
    dt[, Proportion := Tag_count/total_tags]
    dt$Group <- factor(
        c("CDS", "5' UTR", "3' UTR", "Intronic", "Extended 5' UTR", "Extended 3' UTR", "Intergenic"),
        levels = c("Extended 5' UTR", "5' UTR", "CDS", "Intronic", "3' UTR", "Extended 3' UTR", "Intergenic"))
    dt[,.(Group, Tag_count, Proportion)]
}

tagCounts <- parseRSeQCLog(arg.inFile)
tagCounts[, label := arg.label]

g <- ggplot(tagCounts) +
    geom_bar(aes(x = Group, y = Tag_count, fill = Group), stat = 'identity') +
    scale_fill_brewer(palette = 'Set2') +
    labs(title = sprintf("Tag Annotations in %s", arg.label),
         x = "Genomic Region", y = "Tag Count") +
    theme(legend.position = 'none')

ggsave(arg.outFile, g)
