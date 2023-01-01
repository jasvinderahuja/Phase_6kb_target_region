#!/usr/bin/env Rscript
## cmd: JA_make_pre.sam2tsv.R {input} {outputfile}
## input: the output from jvarkit:sam2tsv
## output: an RDS(R data file formatted for insertions)

##setwd("/Volumes/data/projects/PacBio/Pacbio5_exo1ndPII/PB5_tetrads_aligned/sam2tsv/")
## s2t_in <- "PCR_552.sam2tsv.out"
## outputfile <- "PCR_552.sam2tsv.pre.RDS"
now <- Sys.time()
suppressPackageStartupMessages(library(zoo))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(readr))

message(paste0("1. started on ", now))
args<-commandArgs(TRUE)
outputfile<-paste0(args[2])
s2t_in <- args[1]

read.sam2tsv <- function(arg){
        message(arg)
        a<-readr::read_tsv(arg,
            col_names = TRUE, 
            col_types=cols(
###Read-Name	Flag	MAPQ	CHROM	READ-POS0	READ-BASE	READ-QUAL	REF-POS1	REF-BASE	CIGAR-OP
                      `#Read-Name` = col_character(),
                      Flag = col_double(),
                      CHROM = col_character(),
                      `READ-POS0` = col_double(),
                      `READ-BASE` = col_character(),
                      MAPQ = col_character(),
                      `REF-POS1` = col_double(),
                      `REF-BASE` = col_character(),
                      `CIGAR-OP` = col_character(),
                      `READ-QUAL` = col_character()
                    ),
            na = c("", "NA", "."),
            quote=""
            ) %>% 
            dplyr::select(`#READ_NAME`=`#Read-Name`, 
                          FLAG = Flag, 
                          CHROM, 
                          READ_POS=`READ-POS0`, 
                          BASE=`READ-BASE`, 
                          QUAL=MAPQ, 
                          REF_POS=`REF-POS1`, 
                          REF=`REF-BASE`, 
                          OP=`CIGAR-OP`,
                          READ_QUAL = `READ-QUAL`) %>%
            filter(!is.na(READ_POS) & !is.na(BASE) & !is.na(REF)) %>%
            filter(OP!="S") %>%
            mutate(REF_POS=na.locf(REF_POS)) 
        
          message(paste0("a done, with rows = ", nrow(a)))
        return(a)
}


proc.a2 <- read.sam2tsv(s2t_in) 

message(paste0("2. read the ", s2t_in, " ", now))
glimpse(proc.a2)
## tail(proc.a2)

message(paste0("unique reads = ", length(unique(proc.a2$`#READ_NAME`))))

saveRDS(proc.a2, outputfile)

gc()
