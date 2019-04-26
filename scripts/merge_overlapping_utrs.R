#!/usr/bin/env Rscript

library(rtracklayer)
library(GenomicRanges)
library(gUtils)
library(stringr)
library(dplyr)
library(data.table)
library(readr)

if (interactive()) {
    arg.gtf <- "/data/mayrc/data/mca/gff/adult.utrome.e3.t200.f0.999.w500.gtf"
    arg.overlap <- 200
    arg.outFile <- "/scratch/fanslerm/adult.utrome.e3.t200.f0.999.w500.merge.tsv"
}

utrs.gr <- import(arg.gtf)

txs.gr <- utrs.gr[utrs.gr$type == 'transcript',]

txs.ends.gr <- gr.end(txs.gr, ignore.strand=FALSE)


## test.ends.gr <- txs.ends.gr[seqnames(txs.ends.gr) == "chr18",]

## dist.map <- distanceToNearest(test.ends.gr, ignore.strand=FALSE)

## mismatches <- mcols(test.

## clear.idx <- which(mcols(dist.map)$distance >= 200)

## df.out <- data.frame(tx_from=mcols(txs.ends.gr)[clear.idx, 'transcript_id'],
##                      tx_to=mcols(txs.ends.gr)[clear.idx, 'transcript_id'])

## dist.map <- dist.map[-clear.idx,]

## df.hits <- tibble(tx_query=mcols(test.ends.gr)[queryHits(dist.map), 'transcript_id'],
##                   tx_subject=mcols(test.ends.gr)[subjectHits(dist.map), 'transcript_id'])


## df.test <- tibble(tx=test.ends.gr$transcript_id,
##                   gene=test.ends.gr$gene_id,
##                   end=end(test.ends.gr))

## extend downstream and get all intersections on the same strand
overlaps.gr <- gr.findoverlaps(gr.start(txs.ends.gr, width=arg.overlap, clip=FALSE, ignore.strand=FALSE),
                               txs.ends.gr, ignore.strand=FALSE)

## remove self-intersections
overlaps.uniq.gr <- overlaps.gr[mcols(overlaps.gr)$query.id != mcols(overlaps.gr)$subject.id, ]

## for cases where two txs end at the same place, choose the one that came second in the order
overlaps.uniq.gr <- overlaps.uniq.gr[ifelse(strand(overlaps.uniq.gr) == '+',
                                            overlaps.uniq.gr$query.id < overlaps.uniq.gr$subject.id,
                                            overlaps.uniq.gr$query.id > overlaps.uniq.gr$subject.id),]

## Let's work with data.table, since it'll treat all columns identically; also faster
## If multiple intersections, keep only the downstream most tx for the ending tx
remapped.dt <- gr2dt(overlaps.uniq.gr)[,.(subject.id=if(strand == '+') {max(subject.id)} else {min(subject.id)}), by=c('query.id')]

setkey(remapped.dt, 'query.id')

## some chaining may occur, e.g.,
##   tx.1 -> tx.2
##   tx.2 -> tx.3
##   tx.3 -> tx.3
## In this case, since tx.2 will be merged with tx.3, then we also
## need to merge tx.1 to tx.3.  If we think of this as a transition
## matrix, we basically want to apply it repeatedly until we reach
## equilibrium. It's a little tricky with data.table, but it works.
##
## while there are still ending txs that are getting mapped elsewhere
while (nrow(remapped.dt[subject.id %in% query.id]) > 0) {
    ## identify the starting txs
    unmerged.query.ids <- remapped.dt[subject.id %in% query.id, query.id]

    ## get their current ending txs
    subject.ids.old <- remapped.dt[J(unmerged.query.ids), subject.id]

    ## figure out their next ending txs
    subject.ids.new <- remapped.dt[J(subject.ids.old), subject.id]

    ## get row numbers for all the maps we're updating
    unmerged.which <- remapped.dt[subject.id %in% query.id, which=TRUE]

    ## for each row, update the ending txs
    set(remapped.dt, i=unmerged.which, j='subject.id', value=subject.ids.new)
}

## for those that need collapsing map to downstream transcripts
## otherwise, map to itself
df.txmap <- bind_rows(
    tibble(tx_in=mcols(txs.ends.gr)[remapped.dt$query.id, 'transcript_id'],
           tx_out=mcols(txs.ends.gr)[remapped.dt$subject.id, 'transcript_id']),
    tibble(tx_in=mcols(txs.ends.gr)[-remapped.dt$query.id, 'transcript_id'],
           tx_out=mcols(txs.ends.gr)[-remapped.dt$query.id, 'transcript_id'])) %>%
    arrange(-desc(tx_in)) %>%
    mutate(gene_symbol=str_match(tx_out, "(^.*)\\.\\d+$")[,2])

write_tsv(df.txmap, arg.outFile)
