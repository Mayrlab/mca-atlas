library(Biostrings)
library(stringr)

utrome <- readDNAStringSet("data/gff/adult.utrome.e20.t100.f0.999.w300.fasta")
utrome <- utrome[width(utrome) >= 40]
utrome <- subseq(utrome, -40, -1)

writeXStringSet(utrome, "data/fasta/adult.utrome.ends.e20.t100.f0.999.w300.fasta")


utrome.seqs <- as.character(utrome)

PASs <- c("AATAAA", "ATTAAA", "TTTAAA", "AAGAAA", "TATAAA", "AATATA", "AATGAA", "AGTAAA", "AATACA", "CATAAA", "GATAAA", "ACTAAA", "AATAGA", "AAAAAA")

utrome.detected.PASs <- sapply(PASs, function (pas) {
    sapply(utrome.seqs, str_detect, pattern = pas)
})

mean(apply(utrome.detected.PASs, 1, any))

apply(utrome.detected.PASs, 2, mean)

