#!/usr/bin/env Rscript

library(BSgenome.Mmusculus.UCSC.mm10)
library(cleanUpdTSeq)
library(rtracklayer)
library(ggplot2)
library(stringr)

## bypass X11 rendering 
options(bitmapType = 'cairo')

args = commandArgs(trailingOnly=TRUE)

if (length(args) != 4) {
    stop("Incorrect number of arguments!\nUsage:\n> classify_cleavage_sites.R <inFile> <likelihood> <resultsFile> <plotFile>\n")
}

arg.inFile <- args[1]
arg.likelihood <- args[2]
arg.resultsFile <- args[3]
arg.plotFile <- args[4]

formatString <- str_replace(arg.inFile, "bed.gz$", "f%s.bed.gz")

peaks.gr <- import(arg.inFile, genome='mm10')

## Need row names or else cleanUpdTSeq fails
names(peaks.gr) <- elementMetadata(peaks.gr)$name

## Having issue with chr4_GL456350_random, so just leave out for now
peaks.gr <- keepStandardChromosomes(peaks.gr, pruning.mode="coarse")

peaks.features <- buildFeatureVector(peaks.gr, BSgenomeName=Mmusculus, sampleType='unknown',
                                     upstream=40, downstream=30, wordSize=6, method="NaiveBayes",
                                     alphabet=c("ACGT"), replaceNAdistance=30, ZeroBasedIndex=0,
                                     fetchSeq=TRUE)

## load default classifier
data(classifier)

res <- predictTestSet(testSet.NaiveBayes=peaks.features, outputFile=NULL,
                      classifier=classifier, assignmentCutoff=0.5)

passing <- res[res$`prob True` > as.numeric(arg.likelihood), 'PeakName']
gr.filtered <- peaks.gr[peaks.gr$name %in% passing, ]
export.bed(object=gr.filtered, con=sprintf(formatString, arg.likelihood))

gz.out <- gzfile(arg.resultsFile, "w")
write.table(res, gz.out, sep='\t', row.names=FALSE, quote=FALSE)

g <- ggplot(data=res, aes(x=`prob True`)) + geom_histogram() +
    labs(title="Posterior Probabilities for True Poly-A Cleavage Site",
         x="Posterior Probability", y="Genomic Sites")
ggsave(arg.plotFile, g)
