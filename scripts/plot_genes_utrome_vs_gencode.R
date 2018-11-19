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
    stop("Incorrect number of arguments!\nUsage:\n> plot_genes_utrome_vs_gencode.R <utromeRDS> <gencodeRDS> <annotFile> <batch> <outFile>\n")
}

arg.utromeRDS  <- args[1]
arg.gencodeRDS <- args[2]
arg.annotFile  <- args[3]
arg.batch      <- args[4]
arg.outFile    <- args[5]

## Load SCEs
sce.utrome  <- readRDS(arg.utromeRDS)
sce.gencode <- readRDS(arg.gencodeRDS)

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

## determine common subsets
genes.common <- intersect(rownames(sce.utrome), rownames(sce.gencode))
cells.common <- intersect(colnames(sce.utrome), colnames(sce.gencode))

genes.topCounts <- genes.common[order(rowData(sce.utrome[genes.common,])$total_counts, decreasing=TRUE)]
genes.topPctCells <- genes.common[order(rowData(sce.utrome[genes.common,])$n_cells_counts, decreasing=TRUE)]

genes.top50ish <- union(genes.topCounts[1:50], genes.topPctCells[1:50])

clusters <- sapply(cells.common, getCluster)
cellTypes <- sapply(cells.common, getCellType)

for (gene in genes.top50ish) {
    df <- tibble(utrome=counts(sce.utrome)[gene, cells.common],
                 gencode=counts(sce.gencode)[gene, cells.common],
                 cluster=clusters, cellType=cellTypes)

    xylim <- c(0, round(max(df[,c('gencode', 'utrome')])*1.01))
    
    ## Plot base
    g <- ggplot(df, aes(x=gencode, y=utrome)) + geom_point(alpha=0.4) + geom_abline() +
        labs(title=paste0(gene, " Gene Counts in ", arg.batch), x="GENCODE (counts per cell)", y="UTRome (counts per cell)") +
        coord_fixed(xlim=xylim, ylim=xylim)
    ggsave(str_replace(arg.outFile, "png$", paste0(gene, ".png")), g, width=6, height=6)
    
    ## Cell Labels
    g <- ggplot(df, aes(x=gencode, y=utrome, color=cellType)) + geom_point(alpha=0.4) + geom_abline() +
        labs(title=paste0(gene, " Gene Counts in ", arg.batch), color="Cell Type",
             x="GENCODE (counts per cell)", y="UTRome (counts per cell)") +
        coord_fixed(xlim=xylim, ylim=xylim)
    ggsave(str_replace(arg.outFile, "png$", paste0(gene, ".cellTypes.png")), g, width=8, height=6)
    
    ## Labeled only 
    g <- df %>% filter(cellType != "Unknown") %>%
        ggplot(aes(x=gencode, y=utrome, color=cellType)) + geom_point(alpha=0.4) + geom_abline() +
        labs(title=paste0(gene, " Gene Counts in ", arg.batch), color="Cell Type",
             x="GENCODE (counts per cell)", y="UTRome (counts per cell)") +
        coord_fixed(xlim=xylim, ylim=xylim)
    ggsave(str_replace(arg.outFile, "png$", paste0(gene, ".labeled.png")), g, width=8, height=6)
}
