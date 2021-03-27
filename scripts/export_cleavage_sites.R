#!/usr/bin/env Rscript

library(data.table)
library(rtracklayer)
library(GenomicRanges)
library(GenomicFeatures)
library(gUtils)
library(BSgenome.Mmusculus.UCSC.mm10)
library(ggplot2)
library(ggmosaic)
library(readr)
library(dplyr)

## bypass X11 rendering 
options(bitmapType='cairo')

if (interactive()) {
    args <- c("1",
              "data/bed/cleavage-sites/adult.cleavage.e3.t200.bed.gz",
              "data/bed/cleavage-sites/adult.cleavage.e3.t200.f0.0.bed.gz",
              "/data/mayrc/db/mm10/gencode.vM21.annotation.mRNA_ends_found.gff3.gz",
              "/data/mayrc/db/mm10/polyASite.clusters.mm10.bed.gz",
              "1000",
              "5000",
              "adult",
              "e3.t200.f0.0")
} else {
    args = commandArgs(trailingOnly=TRUE)

    if (length(args) != 9) {
        stop("Incorrect number of arguments!\nUsage:\n> export_cleavage_sites.R <cores> <fullInFile> <filteredInFile> <annotFile> <polyASiteFile> <extendUpstream> <extendDownstream> <resultsPrefix> <resultsSuffix>\n")
    }
}

arg.cores <- as.integer(args[1])
arg.fullInFile <- args[2]
arg.filteredInFile <- args[3]
arg.annotFile <- args[4]
arg.pasFile <- args[5]
arg.extUpstream <- as.integer(args[6])
arg.extDownstream <- as.integer(args[7])
arg.outPrefix <- args[8]
arg.outSuffix <- args[9]


## Load References
mm10 <- BSgenome.Mmusculus.UCSC.mm10

annotation.gr <- keepStandardChromosomes(import(arg.annotFile, genome='mm10'),
                                         pruning.mode="coarse")

polyASites.df <- read_tsv(arg.pasFile, na=c("NA", "."),
                          col_names = c("chrom", "start", "stop", "name",
                                        "score", "strand", "pA_signal", "gene"),
                          col_types="ciiciccc")

polyASites.gr <- makeGRangesFromDataFrame(polyASites.df, keep.extra.columns = TRUE,
                                          seqinfo = seqinfo(mm10))

## Load data
sites.gr <- import(arg.filteredInFile, genome='mm10')
names(sites.gr) <- sites.gr$name
sites.gr <- keepStandardChromosomes(sites.gr, pruning.mode="coarse")

sites.all.gr <- import(arg.fullInFile, genome='mm10')
names(sites.all.gr) <- sites.all.gr$name
sites.all.gr <- keepStandardChromosomes(sites.all.gr, pruning.mode="coarse")


## Validated Sites
validated.sites.gr <- (sites.all.gr + 20) %&&%
    gr.end(annotation.gr %Q% (type == 'transcript' &
                              transcript_type == 'protein_coding'),
           ignore.strand = FALSE)

unvalidated.sites.gr <- sites.all.gr %Q% !(name %in% names(validated.sites.gr))

export.bed(validated.sites.gr, sprintf("data/bed/cleavage-sites/%s.validated.%s.bed.gz", arg.outPrefix, arg.outSuffix))

## Supported Sites
supported.sites.gr <- (unvalidated.sites.gr + 20) %&&% (polyASites.gr %Q% (score >= 3))

unsupported.sites.gr <- unvalidated.sites.gr %Q% !(name %in% names(supported.sites.gr))

export.bed(supported.sites.gr, sprintf("data/bed/cleavage-sites/%s.supported.%s.bed.gz", arg.outPrefix, arg.outSuffix))

## Likely Sites
likely.sites.gr <- sites.gr %Q% (name %in% names(unsupported.sites.gr))

unlikely.sites.gr <- unsupported.sites.gr %Q% !(name %in% names(sites.gr))

export.bed(likely.sites.gr, sprintf("data/bed/cleavage-sites/%s.likely.%s.bed.gz", arg.outPrefix, arg.outSuffix))

export.bed(unlikely.sites.gr, sprintf("data/bed/cleavage-sites/%s.unlikely.%s.bed.gz", arg.outPrefix, arg.outSuffix))

## Annotate Sites with Transcriptomic Regions
gene_names.extUTR3 <- function (sites.gr, annotation.gr, ext.width=5000L, mc.cores=1) {
    gns <- gr.val(
        query = sites.gr,
        target = flank(gr.end(annotation.gr %Q% (type == 'three_prime_UTR' &
                                                 transcript_type == 'protein_coding'),
                              ignore.strand = FALSE),
                       width = ext.width, start = FALSE),
        val = c("gene_name"), ignore.strand = FALSE,
        mc.cores = mc.cores)$gene_name
    gns <- as.character(gns)
    ifelse(gns %in% c("", ".", " ", "NA"), NA_character_, gns)
}

gene_names.extUTR5 <- function (sites.gr, annotation.gr, ext.width=5000L, mc.cores=1) {
    gns <- gr.val(
        query = sites.gr,
        target = flank(gr.start(annotation.gr %Q% (type == 'five_prime_UTR' &
                                                   transcript_type == 'protein_coding'),
                                ignore.strand = FALSE),
                       width = ext.width, start = TRUE),
        val = c("gene_name"), ignore.strand = FALSE,
        mc.cores = mc.cores)$gene_name
    gns <- as.character(gns)
    ifelse(gns %in% c("", ".", " ", "NA"), NA_character_, gns)
}

gene_name.label <- function (annotated.gr) {
    regions <- annotated.gr$type
    gnl <- rep_len(NA_character_, length(regions))
    for (i in 1:length(regions)) {
        if (regions[i] != "intergenic") {
            if (regions[i] == "intron") {
                gnl[i] <- elementMetadata(annotated.gr)[i, "transcript"]
            } else {
                gnl[i] <- elementMetadata(annotated.gr)[i, regions[i]]
            }
        }
    }
    gnl
}

annotate.sites <- function (sites.gr, annotation.gr,
                            upstream=5000L, downstream=5000L,
                            mc.cores=1) {
    gr_cols <- c("five_prime_UTR", "exon", "intron", "three_prime_UTR",
                 "intergenic", "extended_five_prime_UTR", "extended_three_prime_UTR")
    annot <- gr.val(sites.gr,
                    annotation.gr %Q% (!is.na(transcript_type) & transcript_type == "protein_coding"),
                    val=c("gene_name"), by=c("type"),
                    by.prefix = NULL, ignore.strand=FALSE,
                    mc.cores=mc.cores)
    annot$extended_three_prime_UTR <- gene_names.extUTR3(sites.gr, annotation.gr, downstream, mc.cores)
    annot$extended_five_prime_UTR <- gene_names.extUTR5(sites.gr, annotation.gr, upstream, mc.cores)
    for (k in gr_cols) {
        if (!(k %in% colnames(elementMetadata(annot)))) {
            elementMetadata(annot)[,k] <- NA_character_
        }
    }
    annot$type <- ifelse(
        is.na(annot$three_prime_UTR), no="three_prime_UTR",
        yes=ifelse(
            is.na(annot$five_prime_UTR), no="five_prime_UTR",
            yes=ifelse(
                is.na(annot$exon), no="exon",
                yes=ifelse(
                    is.na(annot$transcript), no="intron",
                    yes=ifelse(
                        is.na(annot$extended_five_prime_UTR), no="extended_five_prime_UTR",
                        yes=ifelse(
                            is.na(annot$extended_three_prime_UTR), no="extended_three_prime_UTR",
                            yes="intergenic")
                    )
                )
            )
        )
    )
    annot$gene <- gene_name.label(annot)
    annot
}

## Annotate Subsets
validated.sites.gr <- annotate.sites(validated.sites.gr, annotation.gr,
                                     upstream = arg.extUpstream, downstream = arg.extDownstream,
                                     mc.cores = arg.cores)

supported.sites.gr <- annotate.sites(supported.sites.gr, annotation.gr,
                                     upstream = arg.extUpstream, downstream = arg.extDownstream,
                                     mc.cores = arg.cores)

## TODO: Had to add a "+ 1" to give width to the sites, o.w., was returning no intersections
likely.sites.gr <- annotate.sites(likely.sites.gr + 1, annotation.gr,
                                  upstream = arg.extUpstream, downstream = arg.extDownstream,
                                  mc.cores = arg.cores)

if (length(unlikely.sites.gr) > 0) { ## Include possibility of passing all sites
    unlikely.sites.gr <- annotate.sites(unlikely.sites.gr + 1, annotation.gr,
                                        upstream = arg.extUpstream, downstream = arg.extDownstream,
                                        mc.cores = arg.cores)
}

## Tabulate Annotation Categories
annotation_counts.df <- data.frame(list(Class=character(0),
                                        Region=character(0),
                                        Count=numeric(0)))

t_v <- table(validated.sites.gr$type)
annotation_counts.df <- rbind(annotation_counts.df,
                              data.frame(Class="Validated\n(GENCODE)",
                                         Region=names(t_v),
                                         Count=as.numeric(t_v)))

t_s <- table(supported.sites.gr$type)
annotation_counts.df <- rbind(annotation_counts.df,
                              data.frame(Class="Supported\n(PolyASite)",
                                         Region=names(t_s),
                                         Count=as.numeric(t_s)))

t_l <- table(likely.sites.gr$type)
annotation_counts.df <- rbind(annotation_counts.df,
                              data.frame(Class="Likely\n(cleanUpdTSeq)",
                                         Region=names(t_l),
                                         Count=as.numeric(t_l)))

if (length(unlikely.sites.gr) > 0) {
    t_u <- table(unlikely.sites.gr$type)
    annotation_counts.df <- rbind(annotation_counts.df,
                                  data.frame(Class="Unlikely\n(cleanUpdTSeq)",
                                             Region=names(t_u),
                                             Count=as.numeric(t_u)))
}

annotation_counts.df <- annotation_counts.df %>%
    mutate(Region = recode_factor(Region,
                                  extended_five_prime_UTR = "5' UTR (-5 Kb)",
                                  five_prime_UTR = "5' UTR",
                                  exon = "CDS",
                                  intron = "Intronic",
                                  three_prime_UTR = "3' UTR",
                                  extended_three_prime_UTR = "3' UTR (+5 Kb)",
                                  intergenic = "Intergenic"))



## Plot Annotation Categories
g <- ggplot(data = annotation_counts.df) +
    geom_mosaic(aes(x=product(Region, Class), weight=Count, fill=factor(Region)), offset = 0.004) +
    labs(x = "Evidence Group", title = "Cleavage Site Annotations in Mouse Cell Atlas") +
    guides(fill=guide_legend(title = "Region Annotation", reverse = TRUE)) +
    scale_fill_brewer(palette = "Set3")

ggsave(sprintf("qc/cleavage-sites/%s.annotations.%s.png", arg.outPrefix, arg.outSuffix), g,
       width = 12)

## Plot Scores By Evidence
g <- rbind(data.frame(Evidence="Validated\n(GENCODE)", Score=validated.sites.gr$score),
           data.frame(Evidence="Supported\n(PolyASite)", Score=supported.sites.gr$score),
           data.frame(Evidence="Likely\n(cleanUpdTSeq)", Score=likely.sites.gr$score)) %>%
    { 
        if (length(unlikely.sites.gr) > 0) {
            rbind(., data.frame(Evidence="Unlikely\n(cleanUpdTSeq)", Score=unlikely.sites.gr$score))
        } else { . } 
    } %>%
    ggplot() +
    geom_violin(aes(Evidence, Score, fill=Evidence), draw_quantiles = c(0.025, 0.5, 0.975)) +
    scale_y_log10() + theme(legend.position="none") +
    ggtitle("Score Distributions for Mouse Cell Atlas Cleavage Sites By Evidence Level")

ggsave(sprintf("qc/cleavage-sites/%s.site_scores.%s.png", arg.outPrefix, arg.outSuffix), g,
       width = 10)

## Export terminal exon sites
export.bed(grbind(likely.sites.gr %Q% (type == 'three_prime_UTR'),
                  supported.sites.gr %Q% (type == 'three_prime_UTR')),
           sprintf("data/bed/cleavage-sites/%s.utr3.%s.bed.gz", arg.outPrefix, arg.outSuffix))

## Export Extended 3' UTR sites
export.bed(grbind(likely.sites.gr %Q% (type == 'extended_three_prime_UTR'),
                  supported.sites.gr %Q% (type == 'extended_three_prime_UTR')),
           sprintf("data/bed/cleavage-sites/%s.extutr3.%s.bed.gz", arg.outPrefix, arg.outSuffix))
