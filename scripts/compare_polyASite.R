#!/usr/bin/env Rscript

library(BSgenome.Mmusculus.UCSC.mm10)
library(rtracklayer)
library(ggplot2)
library(stringr)
library(GenomicRanges)
library(readr)
library(dplyr)
library(gUtils)

## bypass X11 rendering 
options(bitmapType = 'cairo')

mm10 <- BSgenome.Mmusculus.UCSC.mm10

args = commandArgs(trailingOnly=TRUE)

if (length(args) != 5) {
    stop("Incorrect number of arguments!\nUsage:\n> compare_polyASite.R <pasFile> <score> <window> <bedFile> <outFile>\n")
}

arg.pasFile <- args[1]
arg.score <- as.numeric(args[2])
arg.window <- as.numeric(args[3])
arg.bedFile <- args[4]
arg.outFile <- args[5]

polyASites.df <- read_tsv(arg.pasFile, na=c("NA", "."),
                          col_names = c("chrom", "start", "stop", "name",
                                        "score", "strand", "pA_signal", "gene"),
                          col_types="ciiciccc")

polyASites.gr <- makeGRangesFromDataFrame(polyASites.df, keep.extra.columns = TRUE,
                                          seqinfo = seqinfo(mm10))

bed.gr <- import.bed(arg.bedFile, genome="mm10")
names(bed.gr) <- bed.gr$name

top.gr <- (bed.gr + arg.window) %&&% (polyASites.gr %Q% (score >= arg.score))

export.bed(bed.gr[bed.gr$name %in% names(top.gr),], arg.outFile)
