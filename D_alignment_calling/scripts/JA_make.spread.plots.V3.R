#!/usr/bin/env Rscript
## cmd: Rscript JA_make_spread.plots.R {input_folder} {samplesheet_with_genetics} {outputfile}
## describe
## input_folder: with files having suffix -- "plotData.cut10.RDS"
## samplesheet_with_genetics: needs have PCR, tetrad, spore, Sample, Nat_P, Hyg_P
## outputfile: the pdf file with absolute path
##

message("make.spread.plots.V3.R")
now <- Sys.time()
args<-commandArgs(TRUE)

input_folder <- args[1]
samplesheet_genetics <- args[2]
outputfile <- args[3]


### TEST ###
## begin_forlder = "/Volumes/data/projects/Pacbio/Pacbio2_exo1nd/"
## input_folder <- paste0(begin_forlder, "PB2_alignment/bams/")
## samplesheet_with_genetics <- paste0(begin_forlder, "SampleSheet/Pacbio2_samples_with_error.xlsx")
## outputfile <- paste0(begin_forlder, "PB2_alignment/PB2_alignment_plots.CUT.V2.try.pdf")

#### 

message(paste0("input_folder = ", input_folder,
               " \n samplesheet = ", samplesheet_genetics,
               " \n outputfile = ", outputfile))


suppressPackageStartupMessages(library(zoo))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(ggfortify))
suppressPackageStartupMessages(library(kableExtra))
suppressPackageStartupMessages(library(gridExtra))
suppressPackageStartupMessages(library(readxl))

make_string <- function(x){return(toString(na.omit(unique(x))))}



file_idn <- function(file_RDS){
  file_id <- str_extract(file_RDS, "PCR_(.+)plotData") %>% 
              str_replace_all("plotData|pre|PCR_|\\.", "") %>%
              as.numeric()
  return(file_id)
}



samplesheet_with_genetics <- read_excel(samplesheet_genetics) %>%
                                mutate(PCR = as.numeric(PCR))

data.files <- list.files(input_folder, pattern="plotData.CUT.RDS") %>% 
                    as_tibble() %>%
                    rename(file_ids = value) %>%
                    mutate(PCR = as.numeric(file_idn(file_ids)))

PCR_ids_all <- (data.files %>% arrange(PCR))$file_ids %>% unique() 

sample_genetics_plots <- left_join(samplesheet_with_genetics, data.files, by = c("PCR")) %>%
                              filter(!is.na(file_ids)) %>%
                              arrange(tetrad,spore)
                    


tetrads.exo1 <- sample_genetics_plots$tetrad %>% unique() %>% sort()

getRDS.plot <- function(id_PCR){
  d_plot_RDS <- (readRDS(paste0(input_folder, id_PCR)))$profile_plot_cut+geom_point(size=1)
  return(d_plot_RDS)
}

### MAKE PLOT PER PCR ###
#update_geom_default("point", list(size=1))
theme_set(theme_grey(base_size=6))
# plotfile <- str_replace(input_file, ".RDS", ".pdf") %>% str_replace_all("sam2tsv", "plots")
pdf(paste0(str_replace(outputfile, "plots", "spore_plots")), paper = "letter", width=7.5, height=10.5) 
for(i in PCR_ids_all){
  getRDS.plot(i) %>% print()
}
dev.off()


tetrad_plots <- function(tetrad_no){
#  message(tetrad_no)
  sample_genetics_plots_tetrad <- sample_genetics_plots %>%
                        filter(tetrad == tetrad_no) %>%
                        arrange(spore) %>%
                        select(PCR, Sample, file_ids) %>%
                        distinct()
  d_plots <- lapply(unique(sample_genetics_plots_tetrad$file_ids), getRDS.plot)
#  d_plots[[9]] <-tbl
  return(d_plots)
}

#update_geom_default("point", list(size=1))
theme_set(theme_grey(base_size=6))
plot_details <- paste0("Haplotype Plot (line connects haplotype) \t ", input_folder)
message(paste0(tetrads.exo1, sep=", "))
library(gridExtra)
# plotfile <- str_replace(input_file, ".RDS", ".pdf") %>% str_replace_all("sam2tsv", "plots")
pdf(paste0(outputfile), paper = "letter", width=8, height=10.5) 
    for(i in tetrads.exo1){
      message(i)
      grid.arrange(grobs=tetrad_plots(i), ncol=1, top = paste0("Tetrad = ", i), bottom=plot_details, gp=gpar(fontsize=8))
    }
dev.off()

getRDS.plotData <- function(id_PCR){
  d_plot_RDS <- (readRDS(paste0(input_folder, id_PCR)))$plot_profile_data
  return(d_plot_RDS)
}

all_data <- bind_rows(lapply(unique(sample_genetics_plots$file_ids), getRDS.plotData))
message(nrow(all_data))
saveRDS(all_data, paste0(str_replace(outputfile, ".pdf", ".RDS")))
