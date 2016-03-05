# ChIP-seq (short course)  


##The Course

This course introduces some of the fundamental concepts and tools of ChIP-seq analysis in R/Bioconductor. This course is a shortened version of our full course and covers downstream analysis of ChIP-seq following peak calling and aquisition of QC results.

The course consists of 7 sections.
* Evaluation of ChIP-seq specific QC metrics
* Working with peaks.
* Functional annotation of peaks.
* Motif identification.
* External/Public datasets.
* Exporting data.
* Complex overlaps.
* Simple differential ChIP-seq.

Each section is presented as both HTMl and Rpres markdown ( to allow for intergration of the presentation in the RStudio enviroment itself).  Exercises and answer sheets are included after all subsections to practice techniques and provide future reference examples. 

Course material and exercises are available to view as rendered HTML slides or single page HTML at [http://mrccsc.github.io/ChIPseq_short/](http://mrccsc.github.io/ChIPseq_short).  
All material is available to download under GPL v3 license.

For  information on other courses run by our team see our [github IO page](http://mrccsc.github.io/training/).


## The Team
This course was created and conducted by the MRC Clinical Sciences Centre Bioinformatics Team at Imperial College London, Hammersmith Hospital.  
For more information on the team see our [github IO page](http://mrccsc.github.io/).


This course is free for MRC CSC and Imperial staff and students. If you would like to attend a future course contact thomas.carroll@imperial.ac.uk.

## Setting up.


#### Install R.

R can be installed from the R-project website.  
R 3.1.0 or higher is required for this course.

http://www.r-project.org/

#### Install RStudio.

RStudio can be installed from the R-project website. 

http://www.rstudio.com/

#### Install required packages.

Having downloaded R and RStudio, some additional packages are required (rmarkdown and ggplot2).  
To install these,
* First launch RStudio
* Install the packages in the R console
<pre>
source("https://bioconductor.org/biocLite.R")
biocLite("GenomicRanges")
biocLite("ChIPQC")
biocLite("BSgenome.Mmusculus.UCSC.mm9")
biocLite("ChIPseeker")
biocLite("goseq")
biocLite("rGREAT")
biocLite("AnnotationHub")
</pre>

#### Download the material
The material can either be downloaded as a [zip](https://github.com/mrccsc/ChIPseq_short/archive/master.zip)
<pre>
wget https://github.com/mrccsc/ChIPseq_short/archive/master.zip ./
</pre>
or checked out from our Github repository
https://github.com/mrccsc/ChIPseq_short/
