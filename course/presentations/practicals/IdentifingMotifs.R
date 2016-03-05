## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = F,eval = F)

## ----collapse=T,echo=F---------------------------------------------------
# Load the GenomicRanges Library .. here is suppress messages for a cleaner document
suppressPackageStartupMessages(
  library(GenomicRanges)
)
suppressPackageStartupMessages(
  library(BSgenome)
)
suppressPackageStartupMessages(
  library(BSgenome.Mmusculus.UCSC.mm9)
)


## ----collapse=T----------------------------------------------------------
  library(GenomicRanges)
  library(BSgenome)
  library(BSgenome.Mmusculus.UCSC.mm9)


## ----collapse=T----------------------------------------------------------

melPeak_Rep1 <- read.delim("data/MacsPeaks/mycmelrep1_peaks.xls",sep="\t",comment.char = "#")
melPeak_Rep2 <- read.delim("data/MacsPeaks/mycmelrep2_peaks.xls",sep="\t",comment.char = "#")

melRep1_GR <- GRanges(
                  seqnames=melPeak_Rep1[,"chr"],
                  IRanges(melPeak_Rep1[,"start"],
                  melPeak_Rep1[,"end"]
                  )
                )

mcols(melRep1_GR) <- melPeak_Rep1[,c("abs_summit", "fold_enrichment")]


melRep2_GR <- GRanges(
                  seqnames=melPeak_Rep2[,"chr"],
                  IRanges(melPeak_Rep2[,"start"],
                  melPeak_Rep2[,"end"]
                  )
                )

mcols(melRep2_GR) <- melPeak_Rep2[,c("abs_summit", "fold_enrichment")]

melRep1_GR_InRep1AndRep2 <- melRep1_GR[melRep1_GR %over% melRep2_GR]
melRep1_GR_InRep1AndRep2

## ----collapse=T----------------------------------------------------------
melRep1_GR_InRep1AndRep2 <- melRep1_GR_InRep1AndRep2[ order(melRep1_GR_InRep1AndRep2$fold_enrichment,decreasing=T)
                         ]
top500_melRep1_Rep1AndRep2 <- melRep1_GR_InRep1AndRep2[1:500,]
top500_melRep1_1And2_Resized <- resize(top500_melRep1_Rep1AndRep2,200,fix = "center")
top500_melRep1_1And2_Resized[1:3,]

## ----collapse=T----------------------------------------------------------
genome <- BSgenome.Mmusculus.UCSC.mm9

seqlevelsStyle(top500_melRep1_1And2_Resized) <- "UCSC"

top500_melRep1_1And2_Resized_Seq <- getSeq(genome,top500_melRep1_1And2_Resized)
names(top500_melRep1_1And2_Resized_Seq) <- paste0("peak_",seqnames(top500_melRep1_1And2_Resized),"_",
                                         start(top500_melRep1_1And2_Resized),
                                         "-",
                                         end(top500_melRep1_1And2_Resized))

top500_melRep1_1And2_Resized_Seq[1:4,]

writeXStringSet(top500_melRep1_1And2_Resized_Seq,file="top500_melRep1_1And2.fa")

## ----collapse=T----------------------------------------------------------

start(top500_melRep1_Rep1AndRep2) <- end(top500_melRep1_Rep1AndRep2) <- top500_melRep1_Rep1AndRep2$abs_summit
top500_melRep1_1And2_Resized <- resize(top500_melRep1_Rep1AndRep2,200,fix = "center")


genome <- BSgenome.Mmusculus.UCSC.mm9

seqlevelsStyle(top500_melRep1_1And2_Resized) <- "UCSC"

top500_melRep1_1And2_Resized_Seq <- getSeq(genome,top500_melRep1_1And2_Resized)
names(top500_melRep1_1And2_Resized_Seq) <- paste0("peak_",seqnames(top500_melRep1_1And2_Resized),"_",
                                         start(top500_melRep1_1And2_Resized),
                                         "-",
                                         end(top500_melRep1_1And2_Resized))

writeXStringSet(top500_melRep1_1And2_Resized_Seq,file="abssummit_top500_melRep1_1And2.fa")

