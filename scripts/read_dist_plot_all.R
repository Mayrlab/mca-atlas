#!/usr/bin/env Rscript

#' Parses RSeQC's read_distribution.py output into a data.table object
#' then generates a mosaic plot for where tags mapped in genome for
#' all aligned Mouse Cell Atlas libraries.

suppressMessages({
    require(data.table)
    require(ggplot2)
    require(ggmosaic)
    require(stringr)
})

## bypass X11 rendering 
options(bitmapType = 'cairo')

args = commandArgs(trailingOnly=TRUE)

if (length(args) != 4) {
    stop("Incorrect number of arguments!\nUsage:\n> read_dist_plot_all.R <metadataFile> <qcDir> <selector> <outFile>\n")
}

arg.metadataFile <- args[1]
arg.qcDir <- args[2]
arg.selector <- args[3]
arg.outFile <- args[4]

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

qcFiles <- list.files(arg.qcDir, arg.selector, full.names = TRUE)
names(qcFiles) <- str_extract(qcFiles, "SRR\\d+")

gse.info <- fread(arg.metadataFile)

tagCounts.dt <- data.table()
for (qcFile in qcFiles) {
    tagCounts <- parseRSeQCLog(qcFile)
    srr <- str_extract(qcFile, "SRR\\d+")
    tissue <- gse.info[Run == srr, tissue]
    tagCounts[, Tissue := tissue]
    tagCounts[, SRR := srr]
    tagCounts.dt <- rbind(tagCounts.dt, tagCounts)
}

tagCounts.dt[, Label := sprintf("%s (%s)", Tissue, SRR)]

g <- ggplot(data = tagCounts.dt) +
    geom_mosaic(aes(x=product(Group, Label), weight=Tag_count, fill=factor(Group)), offset = 0.004) +
    labs(x = "Tissue (Library)", title = "Tag Annotations in Mouse Cell Atlas Libraries") +
    guides(fill=guide_legend(title = "Genomic Region", reverse = TRUE)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.4))

ggsave(arg.outFile, g, width = 20)
