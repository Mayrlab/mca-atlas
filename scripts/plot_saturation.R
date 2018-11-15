#!/usr/bin/env Rscript

library(SingleCellExperiment)
library(readr)
library(dplyr)
library(stringr)
library(ggplot2)

## bypass X11 rendering
## "device='png'" is workaround for bug in ggplot2 v3.0.0
options(bitmapType='cairo', device='png')

args = commandArgs(trailingOnly=TRUE)

if (length(args) != 5) {
    stop("Incorrect number of arguments!\nUsage:\n> plot_saturation.R <sceRDS> <annotFile> <label> <batch> <outFile>\n")
}

arg.sceRDS <- args[1]
arg.annotFile <- args[2]
arg.label <- args[3]
arg.batch <- args[4]
arg.outFile <- args[5]

## Load SCE
sce <- readRDS(arg.sceRDS)

## Load Cell Labels
annots <- read_csv(arg.annotFile)

## Retrieve cell cluster label
getCluster <- function (cell.name) {
    if (cell.name %in% annots$Cell.name) {
        cluster = annots %>%
            filter(Cell.name == cell.name) %>%
            select(ClusterID) %>% as.character
    } else {
        cluster = NA_character_
    }
    cluster
}

## Retrieve cell type label
getCellType <- function (cell.name) {
    if (cell.name %in% annots$Cell.name) {
        cellType = annots %>%
            filter(Cell.name == cell.name) %>%
            select(Cell.Anno) %>% as.character %>%
            str_remove(pattern='\\s*\\(.*\\)\\s*$')
    } else {
        cellType = "Unknown"
    }
    cellType
}

## Extract Tx/Gene Counts and Read Counts per Cell
df <- tibble(counts=colSums(counts(sce) > 0), reads=colSums(counts(sce)),
             cluster=sapply(colnames(sce), getCluster), cellType=sapply(colnames(sce), getCellType))

## Plot base
g <- ggplot(df, aes(x = reads, y = counts)) + geom_point(alpha=0.4) +
    labs(title=paste("Saturation Curve for", arg.label, "in Library", arg.batch),
         x="Mapped Reads", y=arg.label)
ggsave(arg.outFile, g, width=6, height=6)

## Cell Labels
g <- ggplot(df, aes(x = reads, y = counts, color = cellType)) + geom_point(alpha=0.4) +
    labs(title=paste("Saturation Curve for", arg.label, "in Library", arg.batch),
         x="Mapped Reads", y=arg.label, color="Cell Type")
ggsave(str_replace(arg.outFile, "png$", "cellTypes.png"), g, width=8, height=6)

## Filtered for labeled only
g <- df %>%
    filter(cellType != "Unknown") %>%
    ggplot(aes(x = reads, y = counts, color = cellType)) + geom_point(alpha=0.4) +
    labs(title=paste("Saturation Curve for", arg.label, "in Library", arg.batch, "(labeled only)"),
         x="Mapped Reads", y=arg.label, color="Cell Type")
ggsave(str_replace(arg.outFile, "png$", "labeledOnly.png"), g, width=8, height=6)

