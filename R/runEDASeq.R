setwd("~/Projects/edaseq/R/")
require("shiny")
require("ShortRead")
require("Rsamtools")
require("DESeq")
require("IRanges")
source("./AllClasses.R"); source("./AllGenerics.R"); source("./functions.R")
source("./methods-BamFileList.R"); source("./methods-FastqFileList.R")
source("./methods-Others.R"); source("./methods-SeqExpressionSet.R")

filelist <- c("../../Masters Thesis/SRR772.bam","../../Masters Thesis/SRR773.bam","../../Masters Thesis/SRR774.bam")
filelist <- c("../../umi/SRR1043197.fastq","../../umi/SRR1043198.fastq","../../umi/SRR1043199.fastq")

runEDASEQ <- function(filelist){
  source("./ui.R"); source("./server.R")
  names(filelist) <- gsub("\\.fastq.*", "", basename(filelist))
  f <- BamFileList(filelist)
  runApp(list(ui=ui.edaseq(f),server=server.edaseq(f)))
}
runEDASEQ(filelist)


files <- list.files(file.path(system.file(package = "leeBamViews"),"bam"), pattern = "bam$", full.names = TRUE)
names(files) <- gsub("\\.bam", "", basename(files))
gt <- gsub(".*/", "", files)
gt <- gsub("_.*", "", gt)
lane <- gsub(".*(.)$", "\\1", gt)
geno <- gsub(".$", "", gt)
pd <- DataFrame(geno=geno, lane=lane, row.names=paste(geno,lane,sep="."))
bfs <- BamFileList(files)
elementMetadata(bfs) <- pd
bfs


