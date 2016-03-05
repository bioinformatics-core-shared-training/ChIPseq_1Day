# Integration with Rmarkdown
<!--
%\VignetteEngine{knitr::rmarkdown}
%\VignetteIndexEntry{Supplementary Materials}
-->


This is an example of integrating a Tracktables report into a RMarkdown document.




```r
library(tracktables)
options(markdown.HTML.header = system.file("misc", "datatables.html", package = "knitr"))
oldFileLocations <- system.file("extdata",package="tracktables")

dir.create(file.path(getwd(),"IGVDirectory"),
           showWarnings = FALSE,recursive = TRUE)
file.copy(oldFileLocations,
          file.path(getwd(),"IGVDirectory"),
          recursive = TRUE)
```

```
## [1] TRUE
```

```r
fileLocations <- file.path(getwd(),"IGVDirectory","extdata")
```

Next the samplesheet of metadata and filesheet of locations is created.

```r
bigwigs <- dir(fileLocations,pattern="*.bw",full.names=TRUE)
intervals <- dir(fileLocations,pattern="*.bed",full.names=TRUE)
bigWigMat <- cbind(gsub("_Example.bw","",basename(bigwigs)),
                   bigwigs)
intervalsMat <- cbind(gsub("_Peaks.bed","",basename(intervals)),
                      intervals)

FileSheet <- merge(bigWigMat,intervalsMat,all=TRUE)
FileSheet <- as.matrix(cbind(FileSheet,NA))
colnames(FileSheet) <- c("SampleName","bigwig","interval","bam")

SampleSheet <- cbind(as.vector(FileSheet[,"SampleName"]),
                     c("EBF","H3K4me3","H3K9ac","RNAPol2"),
                     c("ProB","ProB","ProB","ProB"))
colnames(SampleSheet) <- c("SampleName","Antibody","Species")
```
The tracktables report is created from a call to \texttt{maketracktable}. By default all paths are created relative the directory specified by \texttt{basedirectory}.


```r
HTMLreport <- maketracktable(fileSheet=FileSheet,
                               SampleSheet=SampleSheet,
                               filename="IGVEx3.html",
                               basedirectory=getwd(),
                               genome="mm9")
```

```
## tracktables uses the Datatables javascript libraries.
##             For information on Datatables see http://datatables.net/
```

```r
cat(HTMLreport)
```



<link rel="stylesheet" type="text/css" href="https://cdn.rawgit.com/ThomasCarroll/tracktables-Data/master/js/jquery.datatables.css">
<link rel="stylesheet" type="text/css" href="https://cdn.rawgit.com/ThomasCarroll/tracktables-Data/master/js/tracktables.css">
<script type="text/javascript" language="javascript" src="https://cdn.rawgit.com/ThomasCarroll/tracktables-Data/master/js/jquery.min.js"></script><script type="text/javascript" language="javascript" src="https://cdn.rawgit.com/ThomasCarroll/tracktables-Data/master/js/datatables.js"></script><script type="text/javascript" language="javascript">var loc = window.location.pathname;
var dir = loc.substring(0, loc.lastIndexOf('/'));
var igvtable = [["EBF","EBF","ProB","<a class='table' href='http://localhost:60151/load?file=".concat(dir.concat("/EBFigv.xml&merge=true'".concat(">EBF</a>"))),"<a class='table' href='EBFGI.html'>Intervals</a>"]
,["H3K4me3","H3K4me3","ProB","<a class='table' href='http://localhost:60151/load?file=".concat(dir.concat("/H3K4me3igv.xml&merge=true'".concat(">H3K4me3</a>"))),'No Intervals']
,["H3K9ac","H3K9ac","ProB","<a class='table' href='http://localhost:60151/load?file=".concat(dir.concat("/H3K9acigv.xml&merge=true'".concat(">H3K9ac</a>"))),'No Intervals']
,["RNAPol2","RNAPol2","ProB","<a class='table' href='http://localhost:60151/load?file=".concat(dir.concat("/RNAPol2igv.xml&merge=true'".concat(">RNAPol2</a>"))),'No Intervals']
];
$(document).ready(function() {
    $('#demo').html( '<table cellpadding="0" cellspacing="0" border="0" class="display" id="example"></table>' );
    $('#example').dataTable( {
    "data": igvtable,
columns:[{"title":"SampleName"},
{"title":"Antibody"},
{"title":"Species"},
{"title":"IGV_Link"},
{"title":"Intervals"}]
} );
} );
</script>

<section><div id="tttext"><ul><li>To take advantage of the integration with IGV, <b>IGV must be already running </b>on your machine or can be launched now from this <a class="main" href="http://www.broadinstitute.org/igv/projects/current/igv.php">webstart</a>.</li>
</ul></div>
<div id="demo"></div></section>

