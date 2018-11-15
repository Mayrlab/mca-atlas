#!/usr/bin/env Rscript

library(BSgenome.Mmusculus.UCSC.mm10)
library(rtracklayer)
library(GenomicRanges)
library(GenomicFeatures)
library(gUtils)
library(stringr)
library(pbmcapply)
library(ggplot2)

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

covars <- expand.grid(c("validated", "utr3", "extutr3"), c(3,5,10,15,20), c(50,100,200,300,400,1000))
names(covars) <- c("label", "epsilon", "threshold")

print("Counting polyadenylation signals...")
pas.df <- do.call(rbind, pbmcmapply(grToPAS, label=covars$label, epsilon=covars$epsilon, threshold=covars$threshold,
                                    SIMPLIFY=FALSE, USE.NAMES=FALSE, mc.cores=4))

print("Saving results table...")
outFile <- "qc/cleavage-sites/adult.PAS.counts.tsv"
write.table(pas.df, outFile, row.names=FALSE, col.names=TRUE, quote=FALSE, sep='\t')

print("Plotting results...")
g <- ggplot(data=subset(pas.df, label == "validated"), aes(group=epsilon, color=factor(epsilon))) +
    geom_line(aes(x=threshold, y=sum.any/total)) + coord_cartesian(ylim=c(0,1)) +
    labs(title="Proportion of Validated Sites Containing Known PASs", x="Threshold (Read Counts)", y="Proportion") +
    guides(color=guide_legend(title="Window Width (nt)")) + theme(legend.justification=c(1,0), legend.position=c(0.98,0.02))
ggsave("qc/cleavage-sites/adult.pas.ratios.validated.png", g)

g <- ggplot(data=subset(pas.df, label == "validated"), aes(group=epsilon, color=factor(epsilon))) +
    geom_line(aes(x=threshold, y=sum.any)) +
    labs(title="Number of Validated Sites Containing Known PASs", x="Threshold (Read Counts)", y="Cleavage Sites") +
    guides(color=guide_legend(title="Window Width (nt)")) + theme(legend.justification=c(1,0), legend.position=c(0.98,0.02))
ggsave("qc/cleavage-sites/adult.pas.counts.validated.png", g)

g <- ggplot(data=subset(pas.df, label == "utr3"), aes(group=epsilon, color=factor(epsilon))) +
    geom_line(aes(x=threshold, y=sum.any/total)) + coord_cartesian(ylim=c(0,1)) +
    labs(title="Proportion of 3' UTR Sites Containing Known PASs", x="Threshold (Read Counts)", y="Proportion") +
    guides(color=guide_legend(title="Window Width (nt)")) + theme(legend.justification=c(1,0), legend.position=c(0.98,0.02))
ggsave("qc/cleavage-sites/adult.pas.ratios.utr3.png", g)

g <- ggplot(data=subset(pas.df, label == "utr3"), aes(group=epsilon, color=factor(epsilon))) +
    geom_line(aes(x=threshold, y=sum.any)) +
    labs(title="Number of 3' UTR Sites Containing Known PASs", x="Threshold (Read Counts)", y="Cleavage Sites") +
    guides(color=guide_legend(title="Window Width (nt)")) + theme(legend.justification=c(1,0), legend.position=c(0.98,0.02))
ggsave("qc/cleavage-sites/adult.pas.counts.utr3.png", g)

g <- ggplot(data=subset(pas.df, label == "extutr3"), aes(group=epsilon, color=factor(epsilon))) +
    geom_line(aes(x=threshold, y=sum.any/total)) + coord_cartesian(ylim=c(0,1)) +
    labs(title="Proportion of Extended 3' UTR Sites Containing Known PASs", x="Threshold (Read Counts)", y="Proportion") +
    guides(color=guide_legend(title="Window Width (nt)")) + theme(legend.justification=c(1,0), legend.position=c(0.98,0.02))
ggsave("qc/cleavage-sites/adult.pas.ratios.extutr3.png", g)

g <- ggplot(data=subset(pas.df, label == "extutr3"), aes(group=epsilon, color=factor(epsilon))) +
    geom_line(aes(x=threshold, y=sum.any)) +
    labs(title="Number of Extended 3' UTR Sites Containing Known PASs", x="Threshold (Read Counts)", y="Cleavage Sites") +
    guides(color=guide_legend(title="Window Width (nt)")) + theme(legend.justification=c(1,0), legend.position=c(0.98,0.02))
ggsave("qc/cleavage-sites/adult.pas.counts.extutr3.png", g)

## AAUAAA
g <- ggplot(data=subset(pas.df, label == "validated"), aes(group=epsilon, color=factor(epsilon))) +
    geom_line(aes(x=threshold, y=AATAAA/total)) + coord_cartesian(ylim=c(0,1)) +
    labs(title="Proportion of Validated Sites Containing AAUAAA", x="Threshold (Read Counts)", y="Proportion") +
    guides(color=guide_legend(title="Window Width (nt)")) + theme(legend.justification=c(1,0), legend.position=c(0.98,0.02))
ggsave("qc/cleavage-sites/adult.pas.AAUAAA.ratios.validated.png", g)

g <- ggplot(data=subset(pas.df, label == "validated"), aes(group=epsilon, color=factor(epsilon))) +
    geom_line(aes(x=threshold, y=AATAAA)) +
    labs(title="Number of Validated Sites Containing AAUAAA", x="Threshold (Read Counts)", y="Cleavage Sites") +
    guides(color=guide_legend(title="Window Width (nt)")) + theme(legend.justification=c(1,0), legend.position=c(0.98,0.02))
ggsave("qc/cleavage-sites/adult.pas.AAUAAA.counts.validated.png", g)

g <- ggplot(data=subset(pas.df, label == "utr3"), aes(group=epsilon, color=factor(epsilon))) +
    geom_line(aes(x=threshold, y=AATAAA/total)) + coord_cartesian(ylim=c(0,1)) +
    labs(title="Proportion of 3' UTR Sites Containing AAUAAA", x="Threshold (Read Counts)", y="Proportion") +
    guides(color=guide_legend(title="Window Width (nt)")) + theme(legend.justification=c(1,0), legend.position=c(0.98,0.02))
ggsave("qc/cleavage-sites/adult.pas.AAUAAA.ratios.utr3.png", g)

g <- ggplot(data=subset(pas.df, label == "utr3"), aes(group=epsilon, color=factor(epsilon))) +
    geom_line(aes(x=threshold, y=AATAAA)) +
    labs(title="Number of 3' UTR Sites Containing AAUAAA", x="Threshold (Read Counts)", y="Cleavage Sites") +
    guides(color=guide_legend(title="Window Width (nt)")) + theme(legend.justification=c(1,0), legend.position=c(0.98,0.02))
ggsave("qc/cleavage-sites/adult.pas.AAUAAA.counts.utr3.png", g)

g <- ggplot(data=subset(pas.df, label == "extutr3"), aes(group=epsilon, color=factor(epsilon))) +
    geom_line(aes(x=threshold, y=AATAAA/total)) + coord_cartesian(ylim=c(0,1)) +
    labs(title="Proportion of Extended 3' UTR Sites Containing AAUAAA", x="Threshold (Read Counts)", y="Proportion") +
    guides(color=guide_legend(title="Window Width (nt)")) + theme(legend.justification=c(1,0), legend.position=c(0.98,0.02))
ggsave("qc/cleavage-sites/adult.pas.AAUAAA.ratios.extutr3.png", g)

g <- ggplot(data=subset(pas.df, label == "extutr3"), aes(group=epsilon, color=factor(epsilon))) +
    geom_line(aes(x=threshold, y=AATAAA)) +
    labs(title="Number of Extended 3' UTR Sites Containing AAUAAA", x="Threshold (Read Counts)", y="Cleavage Sites") +
    guides(color=guide_legend(title="Window Width (nt)")) + theme(legend.justification=c(1,0), legend.position=c(0.98,0.02))
ggsave("qc/cleavage-sites/adult.pas.AAUAAA.counts.extutr3.png", g)

