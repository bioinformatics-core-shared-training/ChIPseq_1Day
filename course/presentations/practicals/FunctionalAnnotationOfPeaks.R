## ----collapse=T,echo=F---------------------------------------------------
knitr::opts_chunk$set(echo = F,eval = F)

## ----collapse=T,echo=F---------------------------------------------------
# Load the GenomicRanges Library .. here is suppress messages for a cleaner document
suppressPackageStartupMessages(
  library(GenomicRanges)
)
suppressPackageStartupMessages(
  library(ChIPseeker)
)
suppressPackageStartupMessages(
  library(org.Mm.eg.db)
)
suppressPackageStartupMessages(
  library(TxDb.Mmusculus.UCSC.mm9.knownGene)
)
suppressPackageStartupMessages(
  library(GenomeInfoDb)
)


## ----collapse=T,echo=T---------------------------------------------------
library(GenomicRanges)
library(ChIPseeker)
library(org.Mm.eg.db)
library(TxDb.Mmusculus.UCSC.mm9.knownGene)

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


## ----collapse=T,echo=T---------------------------------------------------
library(GenomeInfoDb)

seqlevelsStyle(melRep1_GR) <- "UCSC"
seqlevelsStyle(melRep2_GR) <- "UCSC"

peakAnno_MelRep1 <- annotatePeak(melRep1_GR, tssRegion=c(-1000, 1000), 
                         TxDb=TxDb.Mmusculus.UCSC.mm9.knownGene, annoDb="org.Mm.eg.db")
peakAnno_MelRep2 <- annotatePeak(melRep2_GR, tssRegion=c(-1000, 1000), 
                         TxDb=TxDb.Mmusculus.UCSC.mm9.knownGene, annoDb="org.Mm.eg.db")

plotAnnoBar(peakAnno_MelRep1)
plotAnnoBar(peakAnno_MelRep2)


plotDistToTSS(peakAnno_MelRep1)
plotDistToTSS(peakAnno_MelRep2)


## ----collapse=T,echo=T---------------------------------------------------
Mel_Rep1_and_Rep2 <- list(peakAnno_MelRep1,peakAnno_MelRep2)
names(Mel_Rep1_and_Rep2) <- c("Rep1","Rep2")
plotAnnoBar(Mel_Rep1_and_Rep2)

## ----collapse=T,echo=T---------------------------------------------------
peakAnno_MelRep1_GR <- as.GRanges(peakAnno_MelRep1)
peakAnno_MelRep2_GR <- as.GRanges(peakAnno_MelRep2)

MelRep1_genesWithPeakInTSS <- unique(peakAnno_MelRep1_GR[peakAnno_MelRep1_GR$annotation == "Promoter",]$SYMBOL)

MelRep1_genesWithPeakInTSS <- unique(peakAnno_MelRep1_GR[peakAnno_MelRep1_GR$annotation == "Promoter",]$SYMBOL)

## ----collapse=T,echo=T---------------------------------------------------
melRep1_GR_InRep1AndRep2 <- melRep1_GR[melRep1_GR %over% melRep2_GR]
melRep1_GR_InRep1AndRep2

## ----collapse=T,echo=F---------------------------------------------------

suppressPackageStartupMessages(
  library(rGREAT)
)


## ----collapse=T,echo=T---------------------------------------------------
library(rGREAT)

## ----collapse=T,echo=T---------------------------------------------------
# Should already be in UCSC format..but just incase
seqlevelsStyle(melRep1_GR_InRep1AndRep2) <- "UCSC"

great_Job <- submitGreatJob(melRep1_GR_InRep1AndRep2,species="mm9")
availableCategories(great_Job)

great_ResultTable = getEnrichmentTables(great_Job,category=
                          "Pathway Data")
names(great_ResultTable)
msigdb_great_ResultTable <- great_ResultTable[["MSigDB Pathway"]]
msigdb_great_ResultTable[1:4,]

