library(Biostrings)
library(tidyverse)

## bypass X11 rendering
## "device='png'" is workaround for bug in ggplot2 v3.0.0
options(bitmapType='cairo', device='png')

utrome <- readDNAStringSet("data/gff/adult.utrome.e3.t200.f0.999.w500.fasta")
utrome <- utrome[width(utrome) >= 50]
utrome <- subseq(utrome, -50, -1)

writeXStringSet(utrome, "data/fasta/adult.utrome.ends.e3.t200.f0.999.w500.fasta")

txid2txname <- read_tsv("data/gff/adult.utrome.txid2txname.e3.t200.f0.999.w500.tsv",
                        col_names=c("id", "name", "symbol", "chr", "utr"))

gencodeUTR <- txid2txname %>% filter(utr == 'GENCODE') %>% pull(id)
upstreamUTR <- txid2txname %>% filter(utr == 'upstream') %>% pull(id)
downstreamUTR <- txid2txname %>% filter(utr == 'downstream') %>% pull(id)

utrome.seqs <- as.character(utrome)

PASs <- c("AATAAA", "ATTAAA", "TTTAAA", "AAGAAA", "TATAAA", "AATATA", "AATGAA", "AGTAAA", "AATACA", "CATAAA", "GATAAA", "ACTAAA", "AATAGA")#, "AAAAAA")

utrome.detected.PASs <- sapply(PASs, function (pas) {
    sapply(utrome.seqs, str_detect, pattern = pas)
})

##mean(apply(utrome.detected.PASs, 1, any))
##apply(utrome.detected.PASs, 2, mean)

any.gencode = mean(apply(utrome.detected.PASs[rownames(utrome.detected.PASs) %in% gencodeUTR,], 1, any))
each.gencode = apply(utrome.detected.PASs[rownames(utrome.detected.PASs) %in% gencodeUTR,], 2, mean)
any.upstream = mean(apply(utrome.detected.PASs[rownames(utrome.detected.PASs) %in% upstreamUTR,], 1, any))
each.upstream = apply(utrome.detected.PASs[rownames(utrome.detected.PASs) %in% upstreamUTR,], 2, mean)
any.downstream = mean(apply(utrome.detected.PASs[rownames(utrome.detected.PASs) %in% downstreamUTR,], 1, any))
each.downstream = apply(utrome.detected.PASs[rownames(utrome.detected.PASs) %in% downstreamUTR,], 2, mean)

## Bin Tian, 2005 data for mouse
each.tian = c(0.5916, 0.1611, 0.0108, 0.0215, 0.0379, 0.0171, 0.0090, 0.0328, 0.0165, 0.0180, 0.0116, 0.0064, 0.0036)

df <- as.data.frame(rbind(each.gencode, each.upstream, each.downstream, each.tian))
df <- cbind(df, list(Source=factor(c("GENCODE", "Upstream", "Downstream", "Tian"),
                                   levels=c("GENCODE", "Upstream", "Downstream", "Tian"),
                                   ordered=TRUE)))

g <- df %>% gather(PAS, Proportion, -Source) %>%
    mutate(PAS = factor(PAS, levels=PASs, ordered=TRUE)) %>%
    ggplot(aes(x=PAS, y=Proportion, fill=Source)) +
    geom_col(position='dodge')

ggsave("qc/PAS_frequencies.utrome.e3.t200.f0.999.w500.png", g, width=12, height=6)
