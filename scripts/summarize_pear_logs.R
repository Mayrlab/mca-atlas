#!/usr/bin/env Rscript

#' Parses PEAR's standard output into a data.table object
#' then generates a mosaic plot for how read pairs successfully
#' merged for all aligned Mouse Cell Atlas libraries.

suppressMessages({
    require(data.table)
    require(ggplot2)
    require(ggmosaic)
    require(stringr)
})

## bypass X11 rendering 
options(bitmapType = 'cairo')

args = commandArgs(trailingOnly=TRUE)

if (length(args) != 3) {
    stop("Incorrect number of arguments!\nUsage:\n> summarize_pear_logs.R <metadataFile> <logsDir> <outFile>\n")
}

arg.metadataFile <- args[1]
arg.logsDir <- args[2]
arg.outFile <- args[3]

dir.create(dirname(arg.outFile), recursive = TRUE, showWarnings = FALSE)

parseSRR <- function (fileName) {
    na.omit(str_match(readLines(fileName), 'Assembled reads.*(SRR\\d+)')[,2])[1]
}

parseReadCounts <- function (fileName) {
    readCounts <- list(
        Assembled = na.omit(str_match(readLines(fileName), 'Assembled reads[^0-9]*([,0-9]+) /')[,2])[1],
        Unassembled = na.omit(str_match(readLines(fileName), 'Not assembled reads[^0-9]*([,0-9]+) /')[,2])[1]
    )
    readCounts <- lapply(readCounts, function (s) { as.numeric(gsub(",", "", s)) })
    as.data.table(readCounts)
}

logFiles <- list.files(arg.logsDir, ".stdout$", full.names = TRUE)

gse.info <- fread(arg.metadataFile)

readCounts.dt <- data.table()
for (logFile in logFiles) {
    readCounts <- parseReadCounts(logFile)
    srr <- parseSRR(logFile)
    tissue <- gse.info[Run == srr, source_name]
    readCounts[, Tissue := tissue]
    readCounts[, SRR := srr]
    readCounts.dt <- rbind(readCounts.dt, readCounts)
}

readCounts.dt[, Label := sprintf("%s (%s)", Tissue, SRR)]

readCounts.mt <- melt(readCounts.dt, id.vars = c('SRR', 'Label', 'Tissue'),
                      measure.vars=c('Assembled', 'Unassembled'),
                      variable.name = "Read_Type", value.name = "Reads")

g <- ggplot(data = readCounts.mt) +
    geom_mosaic(aes(x=product(Read_Type, Label), weight=Reads, fill=factor(Read_Type)), offset = 0.004) +
    labs(x = "Tissue (Library)", title = "Paired-End Read Merging Results in Mouse Cell Atlas Libraries") +
    guides(fill=guide_legend(title = "Read Type")) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.4))

ggsave(arg.outFile, g, width = 16)
