## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = T,eval=T)

## ----cache=T-------------------------------------------------------------
suppressPackageStartupMessages(
  library(GenomicRanges)
)
suppressPackageStartupMessages(
  library(rGREAT)
)
suppressPackageStartupMessages(
  library(AnnotationHub)
)

library(AnnotationHub)
library(GenomicRanges)
library(rGREAT)

ah = AnnotationHub()

query(ah, c("Ebf", "Mus Musculus"))

ebf_proB <- ah[["AH27957"]]

ebf_B <- ah[["AH27958"]]



## ----cache=T-------------------------------------------------------------
peaksGR_List <- c(ebf_proB,ebf_B)
peaksGR_flat <- reduce(peaksGR_List)
common_Ebf_Peaks <- peaksGR_flat[peaksGR_flat %over% ebf_proB & peaksGR_flat %over% ebf_B ]

## ----cache=T-------------------------------------------------------------

seqlevelsStyle(common_Ebf_Peaks) <- "UCSC"

great_Job <- submitGreatJob(common_Ebf_Peaks,species="mm10")
availableCategories(great_Job)

great_ResultTable = getEnrichmentTables(great_Job,category=
                          "Pathway Data")
names(great_ResultTable)
msigdb_great_ResultTable <- great_ResultTable[["MSigDB Pathway"]]


