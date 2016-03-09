## ----eval=F--------------------------------------------------------------
## setwd("/PathToMyDownload/ChIPseq_short/course")
## # e.g. setwd("~/Downloads/ChIPseq_short/course")

## ----eval=F--------------------------------------------------------------
## QCresult <- ChIPQCsample(reads="/pathTo/myChIPreads.bam",
##                          peaks="/pathTo/myChIPpeaks.bed",
##                          genome="mm9",
##                          blacklist = "/pathTo/mm9_Blacklist.bed")

## ----eval=F,echo=T-------------------------------------------------------
## library(ChIPQC)
## load("data/robjects/ChIPQCwithPeaks.RData")
## QCmetrics(res)

## ----eval=T,echo=F-------------------------------------------------------
library(ChIPQC)
load("data/robjects/ChIPQCwithPeaks.RData")
knitr:::kable(QCmetrics(res))

## ----fig.height=5, fig.width=15------------------------------------------
frip(res)
plotFrip(res)

## ----fig.height=5, fig.width=15------------------------------------------
ccplot <- plotCC(res)
ccplot$layers <- ccplot$layers[1]
ccplot

## ----fig.height=5, fig.width=15------------------------------------------
plotSSD(res)+xlim(0,14)

## ----fig.height=5, fig.width=15------------------------------------------
plotSSD(res)+xlim(0.2,0.4)

## ----eval=T,echo=F,  warning=FALSE,collapse=T----------------------------
macsPeaksFiles <- dir("data/MacsPeaks/", full.names=T)
macsPeaks_DF <- read.delim(macsPeaksFiles[1],sep="",comment="#")
knitr:::kable(macsPeaks_DF[1:3,])

## ----eval=T,echo=F,  warning=FALSE,collapse=T----------------------------
macsPeaks_DF <- read.delim(macsPeaksFiles[1],comment="",stringsAsFactors = F)
macsPeaks_DF[1:6,]

## ----eval=T,echo=T, eval=T, echo=T, warning=FALSE------------------------
peakfile <- "data/MacsPeaks/mycch12rep1_peaks.xls"
macsPeaks_DF <- read.delim(peakfile,comment.char="#")
macsPeaks_DF[1:2,]

## ----eval=T,echo=T, eval=T, echo=T, warning=FALSE------------------------
library(GenomicRanges)
macsPeaks_GR <- GRanges(
 seqnames=macsPeaks_DF[,"chr"],
 IRanges(macsPeaks_DF[,"start"],
         macsPeaks_DF[,"end"]
 )
)
macsPeaks_GR

## ----eval=T,echo=T, eval=T, echo=T, warning=FALSE------------------------
seqnames(macsPeaks_GR)
ranges(macsPeaks_GR)

## ----eval=T,echo=T, eval=T, echo=T, warning=FALSE------------------------
mcols(macsPeaks_GR) <- macsPeaks_DF[,c("abs_summit", "fold_enrichment")]
macsPeaks_GR

## ----eval=T,echo=T, eval=T, echo=T, warning=FALSE------------------------
macsPeaks_GR[1:3]
# or macsPeaks_GR[1:3,]

## ----eval=T,echo=T, eval=T, echo=T, warning=FALSE------------------------
macsPeaks_GR[seqnames(macsPeaks_GR) %in% "1"]

## ----eval=T,echo=T, eval=T, echo=T, warning=FALSE------------------------
names(macsPeaks_GR) <- macsPeaks_DF[,"name"]
macsPeaks_GR["mycch12rep1_peak_33496"]

## ----echo=F--------------------------------------------------------------
library(ChIPQC)
library(DESeq2)

## ----eval=T,echo=T, eval=T, echo=T, warning=FALSE------------------------
peakfile <- "data/MacsPeaks/mycch12rep1_peaks.xls"
singlePeakSet <- ChIPQC:::GetGRanges(peakfile, sep="\t", simple=F)
singlePeakSet[1:2,]

## ----commonpeaks_1, eval=T, echo=T, warning=FALSE------------------------
firstPeakSet <- ChIPQC:::GetGRanges("data/MacsPeaks//mycch12rep1_peaks.xls", sep="\t", simple=F)
secondPeakSet <- ChIPQC:::GetGRanges("data/MacsPeaks//mycch12rep2_peaks.xls", sep="\t", simple=F)

## ----commonpeaks_2, eval=T, echo=T, warning=FALSE------------------------
OnlyfirstPeakSet <- firstPeakSet[!firstPeakSet %over% secondPeakSet]
firstANDsecondPeakSets <- firstPeakSet[firstPeakSet %over% secondPeakSet]
length(OnlyfirstPeakSet)
length(firstANDsecondPeakSets)

## ------------------------------------------------------------------------
foldEnrichment <- firstPeakSet$fold_enrichment
# or foldEnrichment <- firstPeakSet[,"fold_enrichment"]
foldEnrichment[1:10]

## ----boxplotOfFE, eval=T, echo=T, fig.height=5, fig.width=9, warning=FALSE----
FirstOnly_FE <- log2(OnlyfirstPeakSet$fold_enrichment)
FirstAndSecond_FE <- log2(firstANDsecondPeakSets$fold_enrichment)

boxplot(FirstOnly_FE,
        FirstAndSecond_FE,
        names=c("Only_in_First","Common_to_first_second"),
        ylab="log2 Fold_Enrichment")

## ----eval=T,echo=T, eval=T, echo=T, warning=FALSE------------------------
firstANDsecondPeakSets <- firstPeakSet[firstPeakSet %over% secondPeakSet]
secondANDfirstPeakSets <- secondPeakSet[secondPeakSet %over% firstPeakSet]

length(firstANDsecondPeakSets)
length(secondANDfirstPeakSets)

## ----eval=T,echo=T, eval=T, echo=T, warning=FALSE------------------------
allPeaks <- c(firstPeakSet,secondPeakSet)
allPeaksReduced <- reduce(allPeaks)
length(allPeaks)
length(allPeaksReduced)

## ----eval=T,echo=T, eval=T, echo=T, warning=FALSE------------------------
commonPeaks <- allPeaksReduced[allPeaksReduced %over% firstPeakSet 
                               & allPeaksReduced %over% secondPeakSet]
length(commonPeaks)

## ----eval=F--------------------------------------------------------------
## source("https://bioconductor.org/biocLite.R")
## biocLite("rGREAT")

## ----eval=T,echo=T, eval=T, echo=T, warning=FALSE,tidy=T-----------------
library(TxDb.Mmusculus.UCSC.mm9.knownGene)
library(org.Mm.eg.db)
library(GenomeInfoDb)
library(ChIPseeker)


## ----eval=T,echo=T, eval=T, echo=T, warning=FALSE,tidy=T-----------------
commonPeaks[1:2,]
seqlevelsStyle(commonPeaks) <- "UCSC"
commonPeaks[1:2,]

## ----eval=T,echo=T, message=F,messages=F, eval=T, echo=T, warning=FALSE,tidy=T----
peakAnno <- annotatePeak(commonPeaks, tssRegion=c(-1000, 1000), 
                         TxDb=TxDb.Mmusculus.UCSC.mm9.knownGene, annoDb="org.Mm.eg.db")

## ----eval=T,echo=T, message=F,messages=F, eval=T, echo=T, warning=FALSE,tidy=T----
class(peakAnno)
peakAnno

## ----eval=T,echo=T, eval=T, echo=T, warning=FALSE,tidy=T-----------------
peakAnno_GR <- as.GRanges(peakAnno)
peakAnno_GR[1:3,]

## ---- eval=T, echo=T, fig.height=5, fig.width=15, warning=FALSE, tidy=T----
plotAnnoBar(peakAnno)

## ----eval=T,echo=T, eval=T, echo=T, warning=FALSE,fig.height=5, fig.width=15,tidy=T----
plotDistToTSS(peakAnno)

## ---- eval=T, echo=T, fig.height=10, fig.width=25, warning=FALSE, tidy=T----
upsetplot(peakAnno, vennpie=F)

## ----eval=T,echo=T, eval=T, echo=T, warning=FALSE,tidy=T-----------------
genesWithPeakInTSS <- unique(peakAnno_GR[peakAnno_GR$annotation == "Promoter",]$geneId)

genesWithPeakInTSS[1:2]

## ----eval=T,echo=T, eval=T, echo=T, warning=FALSE,tidy=T-----------------
allGenes <- unique(
  unlist(
    keys(TxDb.Mmusculus.UCSC.mm9.knownGene, "GENEID")
    )
  )

length(allGenes)

## ----eval=T,echo=T, eval=T, echo=T, warning=FALSE,tidy=T-----------------
allGenesForGOseq <- as.integer(allGenes %in% genesWithPeakInTSS)
names(allGenesForGOseq) <- allGenes
allGenesForGOseq[1:3]

## ------------------------------------------------------------------------
library(KEGG.db)
library(goseq)
xx <- as.list(KEGGPATHID2NAME)
temp <- cbind(names(xx),unlist(xx))
addKeggTogoseq <- function(JX,temp){
  for(l in 1:nrow(JX)){
    if(JX[l,1] %in% temp[,1]){
      JX[l,"term"] <- temp[temp[,1] %in% JX[l,1],2]
      JX[l,"ontology"] <- "KEGG"
    }
    
  }
  return(JX)
}

## ----eval=T,echo=T, eval=T, echo=T, warning=FALSE,tidy=T-----------------
library(goseq)

pwf=nullp(allGenesForGOseq,"mm9","knownGene",plot.fit=FALSE)

Kegg_MycPeaks <- goseq(pwf,"mm9","knownGene",test.cats=c("KEGG"),method="Hypergeometric")

Kegg_MycPeaks <- addKeggTogoseq(Kegg_MycPeaks,temp)

Kegg_MycPeaks[1:4,]


## ----eval=T,echo=T, eval=T, echo=T, warning=FALSE,tidy=T-----------------
library(rGREAT)
seqlevelsStyle(commonPeaks) <- "UCSC"

## ----eval=T,echo=T, eval=T, echo=T,messages=F,message=F,warning=FALSE,tidy=T----
great_Job <- submitGreatJob(commonPeaks,species="mm9")
availableCategories(great_Job)

## ----eval=T,echo=T, eval=T, echo=T, warning=FALSE,tidy=T-----------------
great_ResultTable = getEnrichmentTables(great_Job,category=
                          "Regulatory Motifs")
names(great_ResultTable)
great_ResultTable[["MSigDB Predicted Promoter Motifs"]][1:4,]


## ---- echo=TRUE,collapse=F-----------------------------------------------

library(BSgenome)
library(BSgenome.Mmusculus.UCSC.mm9)
genome <- BSgenome.Mmusculus.UCSC.mm9
seqlevelsStyle(commonPeaks) <- "UCSC"

## ---- echo=TRUE,collapse=F-----------------------------------------------
commonPeaks <- resize(commonPeaks,300,fix="center")
commonPeaks[1:4,]

## ---- echo=TRUE,collapse=F-----------------------------------------------
commonPeaksSequences <- getSeq(genome,GRanges(commonPeaks))
names(commonPeaksSequences) <- paste0("peak_",seqnames(commonPeaks),"_",
                                         start(commonPeaks),
                                         "-",
                                         end(commonPeaks))

commonPeaksSequences[1:2,]

## ---- echo=TRUE,collapse=F-----------------------------------------------
writeXStringSet(commonPeaksSequences,file="consensusPeaks.fa")


## ----eval=T,echo=F,cache=T-----------------------------------------------
library(AnnotationHub)
ah = AnnotationHub()

## ----eval=F,echo=T,cache=T-----------------------------------------------
## library(AnnotationHub)
## ah = AnnotationHub()
## rowResults <- display(ah)

## ----eval=T,echo=T,cache=T-----------------------------------------------
query(ah, c("Myc","BED", "Mus Musculus"))
cmycAnnoHub <- ah[["AH28051"]]
cmycAnnoHub[1:3,]

## ----eval=F,echo=T,cache=T-----------------------------------------------
## library(rtracklayer)
## export.bed(commonPeaks,con = "consensusPeaksForIGV.bed")
## 

## ----eval=F,echo=T,cache=T-----------------------------------------------
## # Something like this
## install.packages( "C:/Users/MYUSERNAME/Desktop/BSgenome.Mmusculus.UCSC.mm10_1.4.0.tar.gz", lib="C:/Users/MYUSERNAME/Desktop/Rlibs", repos = NULL, type = "source"
##                   )
## 
## library(BSgenome.Mmusculus.UCSC.mm10, lib.loc="C:/Users/tcarroll/MYUSERNAME/Rlibs"
##         )

## ---- echo=TRUE,collapse=F-----------------------------------------------
listOfPeaks <- GRangesList(lapply(macsPeaksFiles,
                                  function(x)ChIPQC:::GetGRanges(x,sep="\t",simplify=T)
                                  )
                           )
flattenedPeaks <- reduce(unlist(listOfPeaks))

## ---- echo=TRUE----------------------------------------------------------
matOfOverlaps <- sapply(listOfPeaks,function(x)
                          as.integer(flattenedPeaks %over% x)
                        )

colnames(matOfOverlaps) <- basename(gsub("_peaks\\.xls","",macsPeaksFiles))


mcols(flattenedPeaks) <- matOfOverlaps

flattenedPeaks[1:2,]

## ---- echo=TRUE----------------------------------------------------------
library(limma)
vennCounts(as.data.frame(mcols(flattenedPeaks)))

## ---- echo=T,fig.width=10, fig.height=5----------------------------------
vennDiagram(as.data.frame(elementMetadata(flattenedPeaks)))

## ---- echo=TRUE----------------------------------------------------------

mych12Peaks <- flattenedPeaks[elementMetadata(flattenedPeaks)$mycch12rep1 + elementMetadata(flattenedPeaks)$mycch12rep2 == 2]
mycMelPeaks <- flattenedPeaks[elementMetadata(flattenedPeaks)$mycmelrep1 +  elementMetadata(flattenedPeaks)$mycmelrep2 == 2]

mych12Peaks[1:3,]

## ---- echo=TRUE----------------------------------------------------------

mycMelPeaks_Only <- flattenedPeaks[elementMetadata(flattenedPeaks)$mycmelrep1 + elementMetadata(flattenedPeaks)$mycmelrep2 == 2 &
elementMetadata(flattenedPeaks)$mycch12rep1 + 
                 elementMetadata(flattenedPeaks)$mycch12rep2 == 0]

mycMelPeaks_Only[1,]

## ---- echo=TRUE----------------------------------------------------------

highConfidence_Only <- flattenedPeaks[elementMetadata(flattenedPeaks)$mycmelrep1 +                  elementMetadata(flattenedPeaks)$mycmelrep2 == 2 |
elementMetadata(flattenedPeaks)$mycch12rep1 + 
                 elementMetadata(flattenedPeaks)$mycch12rep2 == 2]


## ---- echo=TRUE----------------------------------------------------------
boxplot(width(highConfidence_Only))
abline(h=400,col="red")

## ---- echo=TRUE----------------------------------------------------------
PeaksToCount <- resize(highConfidence_Only,width = 400,fix = "center")
PeaksToCount[1:2,]

## ---- echo=TRUE----------------------------------------------------------
load("data/robjects/MycCounts.Rdata")
countTable[1:3,]

## ---- echo=TRUE----------------------------------------------------------
library("DESeq2")

colData <- data.frame(SampleName=colnames(countTable[,-c(3,6)]),
                      CellLine=c("ch12","ch12","mel","mel"))

dds <- DESeqDataSetFromMatrix(countData = countTable[,-c(3,6)],
                              colData = colData,
                              design = ~ CellLine,
                              rowRanges=PeaksToCount)

dds <- DESeq(dds)

## ---- echo=TRUE----------------------------------------------------------

test_cellline <- results(dds, contrast=c("CellLine","ch12","mel"),
                         format="GRanges")

UpinMel <- test_cellline[test_cellline$padj < 0.05 & !is.na(test_cellline$padj) 
                         & test_cellline$log2FoldChange > 0]

UpinMel

