#!/usr/bin/env Rscript
## cmd: JA_create.cut.haplotypes.R {inputfile} {outputfile} {samplesheet}
## inputfile: processed RDS with positions formatted for insertions
## outputfile: R list containing abundance_data, haplotype_data, abundance_plot, 
##                               profile_plot, & identifiers = PCR, tetrad, spore)
## samplesheet: needs to have PCR (unique identifier), sample, tetrad, spore  

## OTHER FILES THIS NEEDS -- look at MakeStandard_PB2.R for how to make them
## pos_gps: position groups eg when 8901:8910 = 8905 then that information
## Standard_PB2: this file tells that if position X has genotype Y then it is SNP/WT
## could not confirm 8216:8221, 8223:8229, 8559:8564 or rep(8219, 6), rep(8227,7), rep(8561, 6) as valid positions
## also see MakeStandard_PB2.R

pos_gps <- "/data/ahujajs/helpers/SNP_table_pMJ1101/JA_pos_groups.RDS"
stnds <- "/data/ahujajs/helpers/SNP_table_pMJ1101/JA_standards_PB2.RDS"


### TEST ###
# begin_forlder = "/Volumes/data/projects/Pacbio/Pacbio3D/"

# pos_gps <- "/Volumes/data/helpers/SNP_table_pMJ1101/JA_pos_groups.RDS"
# stnds <- "/Volumes/data/helpers/SNP_table_pMJ1101/JA_standards_PB2.RDS"
# input_folder <- paste0(begin_forlder, "PB3D/bams/")
# inp <- list.files(paste0(begin_forlder, "PB3D/sam2tsv/"), pattern = "pre.RDS")

# max_n_reads.cut.to <- 20
# Pacbio_samplesheet_name <- paste0(begin_forlder, "samplesheet/Pacbio3D_Plate25_28A_samplesheet_genetics.xlsx")
# inp <- list.files(paste0(begin_forlder, "PB3D/sam2tsv/"), pattern = "pre.RDS")
# inputfile<- paste0(begin_forlder, "PB3D/sam2tsv/", inp[3])

###

read_pos <- readRDS(pos_gps)
Standard_PB2 <- readRDS(stnds)

##


dtime <- function(){now <- Sys.time(); return(now)}
message(paste0("begin ", dtime()))
suppressPackageStartupMessages(library(zoo))
suppressPackageStartupMessages(library(tidyverse))

suppressPackageStartupMessages(library(readxl))
## Pacbio_samplesheet <- read_excel("../SampleSheet/Pacbio2_samples_with_error.xlsx")


make_string <- function(x){return(toString(na.omit(unique(x))))}


#message(paste0("started on ", now))
args<-commandArgs(TRUE)
top_n_Reads <- as.numeric(args[4])
Pacbio_samplesheet_name <- paste0(args[3])
outputfile<-paste0(args[2])
inputfile<-paste0(args[1])

message (paste0("most frequent ",  top_n_Reads," haplotypes. "))

Pacbio_samplesheet <- read_excel(Pacbio_samplesheet_name) 

message(paste0("got position grouping data, standards and samplesheet from ", pos_gps, ", ", Pacbio_samplesheet_name, " & ", stnds, " respectively."))


read_RDS <- function(file_RDS){
  ## extracting the file_id and messages for troubleshooting
  file_idn <- str_extract(file_RDS, "PCR_(.+)sam2tsv") %>% 
    str_replace_all("sam2tsv|pre|PCR_|\\.", "") %>%
    as.numeric()
  message(paste0(file_RDS, " file_id= ", file_idn))
  message(paste0("from ", dtime()))
  d_RDS <- readRDS(file_RDS) %>%
    ## filter out the softclipped bases, usually on ends and not aligned to particular spot
    ## fill in the empty REF_POS, with previous REF_POS, it means insertions
    ## left join the read_pos tibble with this new table to id scorable positions (by filtering for them)
    filter(OP!="S") %>%
    mutate(REF_POS=na.locf(REF_POS)) %>%
    left_join(read_pos, by="REF_POS") %>% 
    filter(!is.na(score_pos)) %>% 
    ## now combine the multiple positions that make up one position like 2580:2593 is 2585
    ## note: it also combines intertion into one string
    ungroup() %>%
    arrange(`#READ_NAME`, score_pos, REF_POS) %>%
    group_by(`#READ_NAME`, score_pos) %>%
    summarise(bases = paste0(BASE, sep="", collapse="")) %>%
    ## put the file_id
    mutate(PCR = file_idn)

  return(d_RDS)
          
  invisible(gc())
}

## inputfile <- "../PB2/sam2tsv/PCR_413.sam2tsv.pre.RDS"

proc.a <- read_RDS(inputfile) %>% 
              left_join(Standard_PB2, by = c("score_pos", "bases"))

message(paste0("read the RDS at ", dtime()))

sample_details = subset(Pacbio_samplesheet, as.numeric(PCR) == unique(as.numeric(proc.a$PCR)))
message("file loaded")

invisible(gc())

## data for abundance plot ##
for.abundance.plot <- proc.a %>%
                        group_by(PCR, score_pos, Parent) %>% 
                        tally(name="nReads_abundance")

message("file tally")
invisible(gc())

proc.a.spread <- proc.a %>% 
  select(PCR,`#READ_NAME`, score_pos, Parent) %>%
  spread(score_pos, Parent) %>%
  group_by_at(vars(-c(`#READ_NAME`))) %>%
  summarize(nReads=n()) %>%
  dplyr::filter(nReads > top_n_Reads) %>% 
  dplyr::slice_max(nReads, n = top_n_Reads, with_ties = FALSE)

## proc.a.spread <- subset(proc.a.spr, nReads > top_n_Reads)


### NOW I WANT TO REMOVE reads where any marker is NA
### and set the haplotype with max rows as one,  and so on in decreasing order

proc.a.spread.back <- proc.a.spread   %>% 
  dplyr::arrange(desc(nReads)) %>%
  rowid_to_column(var = "row_n") %>% 
  gather("Position", "Parent", -c("PCR", "nReads", "row_n"))

proc.a.spread.back.NA.rows <- proc.a.spread.back %>%
  ungroup() %>%
  filter(is.na(Parent)) %>%
  select(row_n) %>%
  distinct()

proc.a.spread.back_1 <- proc.a.spread.back %>% 
                          filter(!row_n %in% proc.a.spread.back.NA.rows$row_n) 

begin_haplotypes_n <- sort(unique(proc.a.spread.back_1$nReads))
  
proc.a.spread.back_FINAL <- proc.a.spread.back_1 
##                          filter(as.numeric(n) > max_n_reads.cut.to)

if(length(unique(proc.a.spread.back_FINAL$row_n))>1){
  proc.a.spread.back_FINAL$row_n <- factor(proc.a.spread.back_FINAL$row_n, levels = c(max(proc.a.spread.back_FINAL$row_n):min(proc.a.spread.back_FINAL$row_n)))
} else {
  proc.a.spread.back_FINAL$row_n <- factor(proc.a.spread.back_FINAL$row_n)
}
invisible(gc())

tot_reads <- length(unique(proc.a$`#READ_NAME`))

tetrad_det <- paste0("Sample = ", sample_details$Sample, " \t tetrad = ", sample_details$tetrad, " \t spore = ", sample_details$spore)
primers_etc <- paste0(" \n barcode = ", sample_details$fwd_barcode_5to3, " & ", sample_details$rev_barcode_5to3)
PCR_reads_details <- paste0("Total reads =", tot_reads, 
                            "\t PCR=", unique(proc.a.spread.back_FINAL$PCR), 
                            "\n haplotypes I began with n[supporting reads] = ", paste0(begin_haplotypes_n %>% unique(), collapse=", "),
                            "\n most frequent ", paste0(top_n_Reads), " haplotypes."
                            )



message(paste0("reads=", tot_reads, " sample details ", tetrad_det))

col_ors <- c("S" = "red", "W" = "blue", "black", "orange", "gray")

profile_plot_cut <- ggplot(proc.a.spread.back_FINAL, 
                       aes(x=as.numeric(Position), y=as.factor(row_n), group=as.factor(row_n), color=Parent))+
                    geom_point()+
                    geom_line()+
                    scale_color_manual(values=col_ors, aesthetics = c("fill", "colour"), limits = c("S", "W"))+
                    ylab("Reads/haplotype")+
                    labs(title = tetrad_det,
                            subtitle = PCR_reads_details,
                            caption =  primers_etc)+
                    expand_limits(y=0)+
                    geom_text(aes(label=nReads, x=2000, y=as.factor(row_n)))
  

read_abundance_plot_cut <- ggplot(for.abundance.plot, aes(y=nReads_abundance, x=score_pos, group=Parent, color=Parent))+
  geom_point()+
  expand_limits(y=0)+
  scale_color_manual(values=col_ors, aesthetics = c("fill", "color"), limits = c("S", "W"))+
  ylab("reads per base call")+
  labs(title = tetrad_det,
       subtitle = PCR_reads_details,
       caption =  primers_etc)



results.out <- list(abundance_data = for.abundance.plot,
                    plot_profile_data = proc.a.spread.back_FINAL,
                    abundance_plot = read_abundance_plot_cut,
                    profile_plot_cut = profile_plot_cut,
                    PCR = sample_details$PCR,
                    tetrad = sample_details$tetrad, 
                    spore = sample_details$spore
                    )

message(paste0("done at ", dtime()))
message(outputfile)
saveRDS(results.out, outputfile)

pdf_name <- paste0(str_replace(outputfile, ".RDS", ".pdf"))
pdf(pdf_name, paper = "letter", width=8, height=10.5) 
profile_plot_cut
dev.off()

