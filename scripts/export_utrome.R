#!/usr/bin/env Rscript

library(BSgenome.Mmusculus.UCSC.mm10)
library(rtracklayer)
library(GenomicRanges)
library(GenomicFeatures)
library(gUtils)

## bypass X11 rendering 
options(bitmapType='cairo')

args = commandArgs(trailingOnly=TRUE)

if (length(args) != 7) {
    stop("Incorrect number of arguments!\nUsage:\n> export_utrome.R <cores> <upstreamUTRFile> <downstreamUTRFile> <annotFile> <txLength> <resultsPrefix> <resultsSuffix>\n")
}

arg.cores <- as.integer(args[1])
arg.upstreamUTRFile <- args[2]
arg.downstreamUTRFile <- args[3]
arg.annotFile <- args[4]
arg.txLength <- as.integer(args[5])
arg.outPrefix <- args[6]
arg.outSuffix <- args[7]

## Load References
mm10 <- BSgenome.Mmusculus.UCSC.mm10

## Load Files
cat("Loading input files...\n")
annotation.gr <- keepStandardChromosomes(import(arg.annotFile, genome='mm10'))
upstream.annotation.gr <- keepStandardChromosomes(import(arg.upstreamUTRFile, genome='mm10'))
extended.annotation.gr <- keepStandardChromosomes(import(arg.downstreamUTRFile, genome='mm10'))
augmented.gr <- grbind(annotation.gr, upstream.annotation.gr, extended.annotation.gr)

## Adjust to max transcript size
genes.gr <- augmented.gr %Q% (type == 'gene')
txs.gr <- augmented.gr %Q% (type == 'transcript' & transcript_type == 'protein_coding')
exons.gr <- augmented.gr %Q% (type == 'exon' & transcript_type == 'protein_coding')

## TODO: Figure out where the duplicated exons are coming from
cat("Deduplicating exons...\n")
exons.gr <- exons.gr[!duplicated(exons.gr$ID)]

cat("Intersecting exons with parent transcripts...\n")
children.txs.dt <- gr.findoverlaps(
    exons.gr, txs.gr,
    scol = c("ID"), qcol = c("Parent"), ignore.strand = FALSE,
    return.type = 'data.table', mc.cores = arg.cores)

children.txs.dt <- children.txs.dt[ID == Parent]
children.txs.dt[, width := end - start]
children.txs.dt <- children.txs.dt[order(seqnames, start), ]

clipTranscript <- function (dt) {
    tx.strand <- dt$strand[[1]]
    stopifnot(all(dt[order(start), start] == dt[order(end), start]))
    stopifnot(all(dt[order(start), start] == dt$start))
    if (tx.strand == '+') {
        ## end is fixed
        tx.end <- max(dt$end)

        ## compute cumulative transcript width from end including current exon
        tx.width <- rev(dt[order(-end), cumsum(width)])

        ## select the exon that brings the length over max
        exon.pos <- max(which(tx.width >= arg.txLength), 0)

        if (exon.pos == 0) {
            ## if none of the exons did, then return existing min
            tx.start <- min(dt$start)
        } else {
            ## otherwise, clip so full length is target
            tx.start <- dt$end[[exon.pos]] - (arg.txLength - dt[-(1:exon.pos), sum(width)])
        }
    } else if (tx.strand == '-') {
        ## start is fixed
        tx.start <- min(dt$start)

        ## compute cumulative transcript width from end including current exon
        tx.width <- dt[order(start), cumsum(width)]

        ## select the exon that brings the length over target length
        exon.pos <- min(which(tx.width >= arg.txLength), Inf)

        ## if none of the exons did, then return current end
        if (exon.pos == Inf) {
            tx.end <- max(dt$end)
        } else if (exon.pos == 1) {
            ## if the first exon is already > target length, then clip from there
            tx.end <- dt$start[[exon.pos]] + arg.txLength
        } else {
            ## otherwise, subtract the widths of other exons
            tx.end <- dt$start[[exon.pos]] + (arg.txLength - dt[1:(exon.pos - 1), sum(width)])
        }
    } else {
        warning(sprintf("Strand should be '-' or '+', but found '%s'", tx.strand))
    }
    list(tx.start = tx.start, tx.end = tx.end)
}

cat(sprintf("Truncating transcripts to a max of %d nucleotides...\n", arg.txLength))
txs.bounds.dt <- children.txs.dt[, clipTranscript(.SD), by = subject.id]

utrs.gr <- txs.gr[txs.bounds.dt$subject.id,]
start(utrs.gr) <- txs.bounds.dt$tx.start
end(utrs.gr) <- txs.bounds.dt$tx.end

cat("Deduplicating truncated transcripts...\n")
utrs.gr <- utrs.gr[!duplicated(utrs.gr),]

## Update child entries
cat("Intersecting truncated transcripts with child exons...\n")
children.utrs.dt <- gr.findoverlaps(
    utrs.gr, exons.gr,
    qcol = c("ID"), scol = c("Parent"), ignore.strand = FALSE,
    return.type = 'data.table', mc.cores = arg.cores)
children.utrs.dt <- children.utrs.dt[ID == Parent]

clipChild <- function (tx.start, tx.end, subject.id) {
    cat("\tExtracting exons...", '\n')
    child.gr <- exons.gr[subject.id, ]
    cat("\tUpdating start positions...", '\n')
    start(child.gr) <- pmax.int(start(child.gr), tx.start)
    cat("\tUpdating end positions...", '\n')
    end(child.gr) <- pmin.int(end(child.gr), tx.end)
    child.gr
}

cat("Adjusting exons to truncation bounds...\n")
exons.truncated.gr <- children.utrs.dt[, clipChild(start, end, subject.id)]

## Update genes based on new transcripts
cat("Intersecting genes with truncated transcripts...\n")
genes.txs.dt <- gr.findoverlaps(
    utrs.gr, genes.gr,
    scol = c("ID"), qcol = c("Parent"), ignore.strand = FALSE,
    return.type = 'data.table', mc.cores = arg.cores)
genes.txs.dt <- genes.txs.dt[ID == Parent]

cat("Adjusting gene bounds to match transcript ranges...\n")
## TODO: Improve speed of this method
gene.positions <- genes.txs.dt[, .(start = min(start), end = max(end)), by = subject.id]

genes.truncated.gr <- genes.gr[gene.positions$subject.id, ]
start(genes.truncated.gr) <- gene.positions$start
end(genes.truncated.gr) <- gene.positions$end

## Combine hierarchy of elements
utrome.gr <- grbind(genes.truncated.gr, utrs.gr, exons.truncated.gr)

cat("Creating TxDb object from UTRome annotation...\n")
## TODO: Figure out why some transcripts get dropped!
## TODO: Add more metadata
utrome.txdb <- makeTxDbFromGRanges(utrome.gr, taxonomyId = 10090,
                                   metadata = data.frame(name = c("Genome"), value = c("mm10")))

## Extract exons
cat("Extracting exons from TxDb...\n")
utrome.exons <- exons(utrome.txdb, columns = c('gene_id', 'tx_name', 'exon_name'))
mcols(utrome.exons)$type <- 'exon'
mcols(utrome.exons)$source <- 'GENCODE.vM17_MouseCellAtlas'
mcols(utrome.exons)$transcript_id <- as(mcols(utrome.exons)$tx_name, "CharacterList")
mcols(utrome.exons)$exon_id <- as(mcols(utrome.exons)$exon_name, "CharacterList")
mcols(utrome.exons)$tx_id <- NULL
mcols(utrome.exons)$tx_name <- NULL

## Extract transcripts
cat("Extracting transcripts from TxDb...\n")
utrome.txs <- transcripts(utrome.txdb, columns = c('gene_id', 'tx_id', 'tx_name'))
mcols(utrome.txs)$type <- 'transcript'
mcols(utrome.txs)$source <- 'GENCODE.vM17_MouseCellAtlas'
mcols(utrome.txs)$transcript_id <- as(mcols(utrome.txs)$tx_name, "CharacterList")
mcols(utrome.txs)$transcript_name <- mcols(utrome.txs)$tx_name
mcols(utrome.txs)$tx_id <- NULL
mcols(utrome.txs)$tx_name <- NULL

## Extract genes
cat("Extracting genes from TxDb...\n")
utrome.genes <- genes(utrome.txdb, columns = c('gene_id'))
mcols(utrome.genes)$type <- 'gene'
mcols(utrome.genes)$source <- 'GENCODE.vM17_MouseCellAtlas'


#########################
##  Export Annotation  ##
#########################
cat("Exporting UTRome GTF...\n")
utrome.gtf.file <- sprintf("data/gff/%s.utrome.%s.gtf", arg.outPrefix, arg.outSuffix)
export(utrome.genes, utrome.gtf.file, format = "GTF")
export(utrome.txs, utrome.gtf.file, format = "GTF", append = TRUE)
export(utrome.exons %Q% (strand == '+'), utrome.gtf.file, format = "GTF", append = TRUE)
export(sort(utrome.exons %Q% (strand == '-'), decreasing = TRUE), utrome.gtf.file, format = "GTF", append = TRUE)


########################
##  Export Sequences  ##
########################
cat("Exporting UTRome FASTA...\n")
utrome.seq.file <- sprintf("data/gff/%s.utrome.%s.fasta", arg.outPrefix, arg.outSuffix)
utrome.seq <- extractTranscriptSeqs(mm10, utrome.txdb, use.names=TRUE)
writeXStringSet(utrome.seq, utrome.seq.file, format = "fasta")

cat("All operations complete.")
