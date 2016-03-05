## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = T,eval = T)

## ----collapse=T----------------------------------------------------------
# Load the GenomicRanges Library .. here is suppress messages for a cleaner document
suppressPackageStartupMessages(
  library(GenomicRanges)
  )

melPeak_Rep1 <- read.delim("data/MacsPeaks/mycmelrep1_peaks.xls",sep="\t",comment.char = "#")
melPeak_Rep2 <- read.delim("data/MacsPeaks/mycmelrep2_peaks.xls",sep="\t",comment.char = "#")

melRep1_GR <- GRanges(
                  seqnames=melPeak_Rep1[,"chr"],
                  IRanges(melPeak_Rep1[,"start"],
                  melPeak_Rep1[,"end"]
                  )
                )

mcols(melRep1_GR) <- melPeak_Rep1[,c("abs_summit", "fold_enrichment")]

melRep1_GR

melRep2_GR <- GRanges(
                  seqnames=melPeak_Rep2[,"chr"],
                  IRanges(melPeak_Rep2[,"start"],
                  melPeak_Rep2[,"end"]
                  )
                )

mcols(melRep2_GR) <- melPeak_Rep2[,c("abs_summit", "fold_enrichment")]

melRep2_GR


## ---- warnings=F,collapse=T----------------------------------------------
# Number of peaks

length(melRep1_GR)

length(melRep2_GR)

# Number of peaks on chromosome 4

# Using table on logical vector

table(seqnames(melRep1_GR) %in% "4")

table(seqnames(melRep2_GR) %in% "4")

# Indexing and recounting

length(melRep1_GR[seqnames(melRep1_GR) %in% "4"])

length(melRep2_GR[seqnames(melRep2_GR) %in% "4"])


# Number of peaks on chromosome 4 and with 5 fold enrichment above input

# Using table on logical vector

table(seqnames(melRep1_GR) %in% "4" & melRep1_GR$fold_enrichment > 5)

table(seqnames(melRep2_GR) %in% "4" & melRep2_GR$fold_enrichment > 5)

# Indexing and recounting

length(melRep1_GR[seqnames(melRep1_GR) %in% "4" & melRep1_GR$fold_enrichment > 5])

length(melRep2_GR[seqnames(melRep2_GR) %in% "4" & melRep2_GR$fold_enrichment > 5])


## ---- warnings=F,collapse=T----------------------------------------------

# Using table

table(melRep1_GR %over% melRep2_GR)

# Using index and recounting
length(melRep1_GR[melRep1_GR %over% melRep2_GR])


# Using table

table(!melRep1_GR %over% melRep2_GR)

# Using index and recounting
length(melRep1_GR[!melRep1_GR %over% melRep2_GR])

commonMelPeaks <- melRep1_GR[melRep1_GR %over% melRep2_GR]

## ---- warnings=F,collapse=T----------------------------------------------

# Using table

melRep1_GRSummits <- melRep1_GR

start(melRep1_GRSummits) <- end(melRep1_GRSummits) <- melRep1_GR$abs_summit

melRep1_GRSummits

