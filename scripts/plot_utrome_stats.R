#!/usr/bin/env Rscript

library(BSgenome.Mmusculus.UCSC.mm10)
library(rtracklayer)
library(GenomicRanges)
library(GenomicFeatures)
library(gUtils)
library(ggplot2)
library(stringr)

## bypass X11 rendering
## "device='png'" is workaround for bug in ggplot2 v3.0.0
options(bitmapType='cairo', device='png')

args = commandArgs(trailingOnly=TRUE)

if (length(args) != 4) {
    stop("Incorrect number of arguments!\nUsage:\n> plot_utrome_stats.R <cores> <UTRome GTF> <resultsPrefix> <resultsSuffix>\n")
}

arg.cores <- as.integer(args[1])
arg.utromeFile <- args[2]
arg.outPrefix <- args[3]
arg.outSuffix <- args[4]

#arg.cores <- 4L
#arg.utromeFile <- "data/gff/adult.utrome.e20.t100.f0.999.w300.gtf"
#arg.outPrefix <- "mock"
#arg.outSuffix <- "e20.t100.f0.999.w300"

## Load Files
cat("Loading input files...\n")
utrome.gr <- keepStandardChromosomes(import(arg.utromeFile, genome='mm10'))

## Simple Counts
counts.type <- table(utrome.gr$type)
cat(sprintf("Found %d genes, %d transcripts, and %d exons.\n",
            counts.type[['gene']], counts.type[['transcript']], counts.type[['exon']]))

## Transcripts Per Gene
counts.txsPerGene <- as.data.frame(table(table(mcols(utrome.gr)[utrome.gr$type == 'transcript', 'gene_id'])))
names(counts.txsPerGene) <- c("utrs", "genes")

g <- ggplot(counts.txsPerGene, aes(x=utrs, y=genes)) + geom_bar(stat = 'identity') +
    labs(title = "Transcripts Per Gene in UTRome", x = "Transcripts Counts", y = "Gene Counts")
ggsave(sprintf("qc/utrome/%s.utrome.txsPerGene.%s.png", arg.outPrefix, arg.outSuffix), g)

## Exons Per Transcript
counts.exsPerTxs <- as.data.frame(table(table(mcols(utrome.gr)[utrome.gr$type == 'exon', 'transcript_id'])))
names(counts.exsPerTxs) <- c("exons", "transcripts")

g <- ggplot(counts.exsPerTxs, aes(x=exons, y=transcripts)) + geom_bar(stat = 'identity') +
    labs(title = "Exons Per Transcript in UTRome", x = "Exon Counts", y = "Transcript Counts")
ggsave(sprintf("qc/utrome/%s.utrome.exonsPerTx.%s.png", arg.outPrefix, arg.outSuffix), g)

## Lengths
lengths.df <- data.frame(length = width(utrome.gr), type = utrome.gr$type)

g <- ggplot(lengths.df, aes(x = type, y = length, fill = type)) +
    geom_violin() + scale_y_log10() +
    labs(title = "Lengths in Truncated UTRome", x = "Annotation", y = "Length")
ggsave(sprintf("qc/utrome/%s.utrome.lengths.%s.png", arg.outPrefix, arg.outSuffix), g)

## Distances
txs.gr <- utrome.gr %Q% (type == 'transcript')
dists.dt <- as.data.table(distanceToNearest(gr.end(txs.gr, ignore.strand = FALSE)))
dists.dt[, chromosome := str_replace(as.character(seqnames(txs.gr[queryHits])), "chr", "")]

g <- ggplot(dists.dt, aes(x = distance + 1)) +
    stat_ecdf(geom = 'step') + scale_x_log10() +
    labs(title="Empirical CDF of Distances between Cleavage Sites",
         x="Distance (nt)", y="Proportion Transcripts (Proportion)")
ggsave(sprintf("qc/utrome/%s.utrome.dists.%s.png", arg.outPrefix, arg.outSuffix), g,
       height = 8, width = 8)

g <- ggplot(dists.dt, aes(x = distance + 1, group = chromosome, color = chromosome)) +
    stat_ecdf(geom = 'step') + scale_x_log10() +
    labs(title="Empirical CDF of Distances between Cleavage Sites",
         x="Distance (nt)", y="Transcripts (Proportion)", color="Chromosome")
ggsave(sprintf("qc/utrome/%s.utrome.distsByChr.%s.png", arg.outPrefix, arg.outSuffix), g,
       height = 8, width = 10)

## Cleavage Sites per Gene
sites.dt <- gr2dt(gr.end(utrome.gr %Q% (type == 'transcript'), ignore.strand = FALSE))

counts.sitesPerGene <- sites.dt[, length(unique(start)), by = gene_id]
names(counts.sitesPerGene) <- c("gene", "utrs")

g <- ggplot(counts.sitesPerGene, aes(utrs)) + geom_bar() +
    labs(title = "Cleavage Sites Per Gene in UTRome", x = "Cleavage Sites", y = "Genes")
ggsave(sprintf("qc/utrome/%s.utrome.sitesPerGene.%s.png", arg.outPrefix, arg.outSuffix), g)

g <- ggplot(counts.sitesPerGene, aes(utrs)) + geom_bar() +
    labs(title = "Cleavage Sites Per Gene in UTRome", x = "Cleavage Sites", y = "Genes") +
    scale_y_log10()
ggsave(sprintf("qc/utrome/%s.utrome.sitesPerGene.log.%s.png", arg.outPrefix, arg.outSuffix), g)

## Transcripts Per Cleavage Site
counts.txsPerSite <- sites.dt[, length(transcript_id), by = .(gene_id, start)]
names(counts.txsPerSite) <- c('gene', 'start', 'txs')

g <- ggplot(counts.txsPerSite, aes(txs)) + geom_bar() +
    labs(title = "Transcripts Per Cleavage Site in UTRome", x = "Transcripts", y = "Cleavage Sites")
ggsave(sprintf("qc/utrome/%s.utrome.txsPerSite.%s.png", arg.outPrefix, arg.outSuffix), g)

g <- ggplot(counts.txsPerSite, aes(txs)) + geom_bar() +
    labs(title = "Transcripts Per Cleavage Site in UTRome", x = "Transcripts", y = "Cleavage Sites") +
    scale_y_log10()
ggsave(sprintf("qc/utrome/%s.utrome.txsPerSite.log.%s.png", arg.outPrefix, arg.outSuffix), g)
