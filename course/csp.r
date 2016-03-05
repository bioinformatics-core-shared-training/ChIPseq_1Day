testC <- readChar("~/Downloads/tommy/short_ChIPseq/chipseq_Intro_2.md", file.info("~/Downloads/tommy/short_ChIPseq/chipseq_Intro_2.md")$size)
testCC <- strsplit(testC,"\n")
presBreaks <- grep("^=*=$",testCC[[1]])
testCC[[1]][grep("^=*=$",testCC[[1]])-1] <- paste0("##",
                                                   testCC[[1]][grep("^=*=$",testCC[[1]])-1])
testCC[[1]][grep("^=*=$",testCC[[1]])] <- "\n"
lapply(testCC[[1]], cat, "\n", file="~/Downloads/tommy/short_ChIPseq/chipseq_Intro_2_SinglePage.md", append=TRUE)
getwd()
