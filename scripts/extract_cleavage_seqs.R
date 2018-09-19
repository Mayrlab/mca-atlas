#!/usr/bin/env Rscript

library(BSgenome.Mmusculus.UCSC.mm10)
library(rtracklayer)
library(GenomicRanges)
library(GenomicFeatures)
library(gUtils)
library(stringr)
library(pbmcapply)

## bypass X11 rendering 
options(bitmapType='cairo', device='png')

## Load References
mm10 <- BSgenome.Mmusculus.UCSC.mm10

fileStr <- "data/bed/cleavage-sites/adult.%s.e%s.t%s.f0.999.bed.gz"

PASs <- c("AATAAA", "ATTAAA", "TTTAAA", "AAGAAA", "TATAAA", "AATATA", "AATGAA", "AGTAAA", "AATACA", "CATAAA", "GATAAA", "ACTAAA", "AATAGA", "AAAAAA")

grToPAS <- function (label, epsilon, threshold, dist=40) {
    file <- sprintf(fileStr, label, epsilon, threshold)
    gr <- keepStandardChromosomes(import(file, genome='mm10'), pruning.mode='coarse')
    ends.gr <- gr.end(gr, width=dist, clip=FALSE, ignore.strand=FALSE)
    ends.seqs <- as.character(getSeq(mm10, ends.gr))
    detected.PASs <- sapply(PASs, function (pas) {
        sapply(ends.seqs, str_detect, pattern=pas)
    })
    counts.PASs <- apply(detected.PASs, 2, sum)
    df <- data.frame(t(unlist(counts.PASs)))
    df$label <- label
    df$epsilon <- epsilon
    df$threshold <- threshold
    df$sum.any <- sum(apply(detected.PASs, 1, any))
    df$total <- nrow(detected.PASs)
    df
}

covars <- expand.grid(c("validated", "utr3", "extutr3"), c(3,5,20), c(50,100,300,1000))
names(covars) <- c("label", "epsilon", "threshold")

pas.df <- do.call(rbind, pbmcmapply(grToPAS, label=covars$label, epsilon=covars$epsilon, threshold=covars$threshold,
                                    SIMPLIFY=FALSE, USE.NAMES=FALSE, mc.cores=4))
outFile <- "qc/cleavage-sites/adult.PAS.counts.tsv"
write.table(pas.df, outFile, row.names=FALSE, col.names=TRUE, quote=FALSE, sep='\t')
