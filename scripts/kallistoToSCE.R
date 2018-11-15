#!/usr/bin/env Rscript

library(scater)
library(readr)
library(dplyr)

## bypass X11 rendering
## "device='png'" is workaround for bug in ggplot2 v3.0.0
##options(bitmapType='cairo', device='png')

################################################################################
                                        # Load Argument Values
################################################################################
args = commandArgs(trailingOnly=TRUE)

if (length(args) != 4) {
    stop("Incorrect number of arguments!\nUsage:\n> kallistoToSCE.R <inputPath> <txToGeneFile> <txsOutFile> <genesOutFile>\n")
}

arg.inputPath <- args[1]
arg.txToGeneFile <- args[2]
arg.txsOutFile <- args[3]
arg.genesOutFile <- args[4]

################################################################################
                                        # Load Data
################################################################################
loadSCEFromDir <- function (sample.path) {
    sample.names <- list.dirs(path=sample.path, full.names=FALSE, recursive=FALSE)
    sample.dirs <- list.dirs(path=sample.path, full.names=TRUE, recursive=FALSE)
    readKallistoResults(samples=sample.names, directories=sample.dirs, read_h5=TRUE)
}

## Load kallisto results as SingleCellExperiment
sce <- loadSCEFromDir(arg.inputPath)

## Load transcript annotations
tx2gene <- read_tsv(arg.txToGeneFile, col_names=c("transcript_id", "gene_id", "gene_symbol", "chromosome"))
rowData(sce) <- merge(rowData(sce), tx2gene, by.x="feature_id", by.y="transcript_id", sort=FALSE)

## Note mitochondrial transcripts
mitoTXs <- which(rowData(sce)$chromosome == "chrM")

## QC Metrics
sce <- calculateQCMetrics(sce, feature_controls=list(mito=mitoTXs))

################################################################################
                                        # Gene Counts
################################################################################
sce.genes <- summariseExprsAcrossFeatures(sce, exprs_values="counts", summarise_by="gene_symbol")
names(rowData(sce.genes)) <- c("gene_symbol")

## WARNING: There are 12 gene symbols that each have 2 ENSMUSG entries
## For now, we are just taking whichever is first.
genes.tbl <- tx2gene %>%
    select(-transcript_id) %>%
    distinct(gene_symbol, .keep_all=TRUE)

## Eliminate Gene Symbol ambiguity
rowData(sce.genes) <- merge(rowData(sce.genes), genes.tbl, by="gene_symbol", sort=FALSE)

## Note mitochondrial genes
mitoGenes <- which(rowData(sce.genes)$chromosome == "chrM")

## QC Metrics
sce.genes <- calculateQCMetrics(sce.genes, feature_controls=list(mito=mitoGenes))

################################################################################
                                        # Export Data
################################################################################
saveRDS(sce, arg.txsOutFile)
saveRDS(sce.genes, arg.genesOutFile)
