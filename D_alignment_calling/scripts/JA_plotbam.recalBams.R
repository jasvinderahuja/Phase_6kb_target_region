#!/usr/bin/env Rscript
library(tidyverse)


args<-commandArgs(TRUE)
bamsfolder <- args[1]
plotfile <- args[2]
fasta_file <- args[3]

attr_file <- str_replace(fasta_file, ".fasta", "attrib.txt")
#BiocManager::install("Gviz")

getBams <- function ( bamdirs ) {
  bais <- list.files(path=bamsfolder, pattern="fin.recal.bam.bai", full.names=T, recursive=FALSE)
  bams<-gsub(".bai", "", bais)
  return(bams)
}
getdirname <- function(pathtofile) {
  dirnsplit <- strsplit(dirname(pathtofile), "/")
  dirname1 <- dirnsplit[[1]][length(dirnsplit[[1]])]
  return(dirname1)
}

getdirnames <- function(paths) {
  return(lapply(paths, getdirname))
}

plotbamd<- function (sam){
  require(data.table)
  require(Gviz)
  require(rtracklayer)
  require(GenomicFeatures)
  require(Biostrings)
  require(Rsamtools)
  require(readr)
  require(reshape2)
  require(ggplot2)
  require(stringr)
  require(dplyr)
  require(tidyr)
  require(VariantAnnotation)
  
  options(ucscChromosomeNames=FALSE)
  S4921<-readDNAStringSet(fasta_file)
  STE50toFUS1_S4921 <- read_delim(attr_file, 
                                  "\t", escape_double = FALSE, trim_ws = TRUE)
  sTrack<-SequenceTrack(S4921, chromosome="STE50toFUS1_S4921", name="Seq")
  axisTrack<- GenomeAxisTrack(range=IRanges(start = STE50toFUS1_S4921$start, end=STE50toFUS1_S4921$end), labelPos = "below")
  gens<-subset(STE50toFUS1_S4921, groups=="his4::URA3rev-tel1-ARG4")
  feats<-subset(STE50toFUS1_S4921,groups=="markers")
  mTrack <- AnnotationTrack(start = feats$start, width = feats$width, chromosome = feats$`# seqname`, strand = feats$strand,name = "Pri,Ec, DSB", labelPos = "below",stacked="dense",  shape="box", collapse=TRUE,col = NULL,)
  gTrack <- AnnotationTrack(start = gens$start, width = gens$width, chromosome = gens$`# seqname`, strand = gens$strand, group = gens$genes, feature = gens$genes, name = "h::UtA",featureAnnotation = "feature", fontcolor.feature = "White", cex.feature = 0.4, his4 = "black", ura3 = "darkgreen",arg4 = "darkgreen", tel1 = "darkred", collapse=TRUE, stacked="dense",col = NULL)
  fl<-paste0(basename(sam))
  gTrack <- AnnotationTrack(start = gens$start, width = gens$width, chromosome = gens$`# seqname`, strand = gens$strand, group = gens$genes, feature = gens$genes, name = "h::UtA",featureAnnotation = "feature", fontcolor.feature = "White", cex.feature = 0.4, his4 = "black", ura3 = "darkgreen",arg4 = "darkgreen", tel1 = "darkred", collapse=TRUE, stacked="dense")
  dr <- getdirname(sam)
  assign = paste0(dr,".",fl)
  S1AlTrack<-AlignmentsTrack(sam, isPaired=TRUE,name = assign)
  return(c( ass = assign,
            Path = sam, 
            dir = dr, 
            file = fl, 
            axis = axisTrack, 
            seq  = sTrack, 
            gene = gTrack, 
            m = mTrack,
            Alignment = S1AlTrack)
        )
}
plottheBams <- function(i) {
  res <- tryCatch(bamplot<-plotbam(i), error = function(e) e)
  if (inherits(res, "error")) (assign = paste0(getdirname(i),".",basename(i), "error"))
  return(bamplot)
}

plotallbams <-function(){
  a<-getBams(bamsfolder)
  bamPlots<-lapply(a, plotbamd)
  saveRDS(bamPlots, file=paste0(bamsfolder, "bamPlots.rds", collapse = NULL))
}
plotBams <- plotallbams()
plotBams <- readRDS(paste0(bamsfolder, "bamPlots.rds", collapse = NULL))
filetibble<-cbind(unlist(lapply(plotBams,"[",3)), 
                  unlist(lapply(plotBams,"[",4)), 
                  unlist(lapply(plotBams,"[",2)), 
                  unlist(lapply(plotBams,"[",1)))
colnames(filetibble) <- c("dir", "file", "path", "assign")
filetibble<-dplyr::as_tibble(filetibble)

finditem <- function(plotBams, dr) {
  for (i in 1:length(plotBams)) {
    if(plotBams[[i]]$Path == dr) { b = i; break }
  }
  return(b)
}
plotbam <- function(dr, st, en) {
      i <- finditem(plotBams, dr)
      plotTracks(c(plotBams[[i]]$axis, 
                  plotBams[[i]]$seq, 
                  plotBams[[i]]$gene, 
                  plotBams[[i]]$m, 
                  plotBams[[i]]$Alignment), 
                sizes=c(1,1,2,1, 40), 
                chromosome="STE50toFUS1_S4921", 
                from=st, to=en, 
                stacked="dense", 
                cex.feature = 0.5, 
                showOverplotting = TRUE
                )
}

filetoget <- function(ngs) {
  fl <- filetibble %>%
  dplyr::filter(dir == ngs) %>%
  dplyr::select(file)
  return(fl)
}
a<-getBams(bamsfolder)
pdf(file=plotfile, width=11, height = 8, pointsize= 1/300)
for(i in a) {
      message(paste(i))
      res <- tryCatch(p<-plotbam(i, 1500, 10000), error = function(e) e)
      if(inherits(res, "error")) next
      print(p)
}
dev.off()








