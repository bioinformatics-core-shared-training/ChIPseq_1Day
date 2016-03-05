## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = F,eval=F)

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


## ----cache=T,eval=F------------------------------------------------------
## seqlevelsStyle(common_Ebf_Peaks) <- "UCSC"
## common_Ebf_Peaks_Resized <- resize(common_Ebf_Peaks,200,fix = "center")
## 
## 
## genome <- BSgenome.Mmusculus.UCSC.mm10
## 
## seqlevelsStyle(common_Ebf_Peaks) <- "UCSC"
## 
## common_Ebf_Peaks_Seq <- getSeq(genome,common_Ebf_Peaks)
## 
## names(common_Ebf_Peaks_Seq) <- paste0("peak_",seqnames(common_Ebf_Peaks),"_",
##                                          start(common_Ebf_Peaks),
##                                          "-",
##                                          end(common_Ebf_Peaks))
## 
## writeXStringSet(common_Ebf_Peaks_Seq,file="common_Ebf_Peaks.fa")

