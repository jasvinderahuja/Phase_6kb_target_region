#!/usr/bin/env Rscript
### syntax == Rscript Make_fasta_folders_seqkit.R PB5_samplesheet_39_fixed.xlsx PB5.fastq Fwd_barcode_scan.tsv Rev_barcode_scan.tsv
### jasvinderahuja@gmail.com
### date 2/16/2021

args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args) < 4) {
  stop(paste0('argument 1 = excel samplesheet with "PCR", "fwd_barcode_5to3" and "rev_barcode_5to3". \n', 
              'I give every PCR a unique number that helps to demultiplex that number goes in the PCR column. \n',
              'argument 2 = the input.fasta or fastq file. \n',
              'argument 3 = Fwd barcode result from seqkit. \n',
              'argument 4 = Rev barcode result from seqkit. \n'
  ), call.=FALSE)
} 

suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(foreach))
suppressPackageStartupMessages(library(itertools))
suppressPackageStartupMessages(library(doMC))
suppressPackageStartupMessages(library(readxl))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(stringr))

samplesheet.file <- args[1]
fasta_run <- args[2]
barcodeID.fwd <- args[3]
barcodeID.rev <- args[4]
out_fastq <- paste0(str_replace(basename(samplesheet.file), ".xlsx", "_binned_fastqs"))
outputfile2 <- paste0(samplesheet.file, ".with.samples.RDS")

cat(paste0(
  "\n",
  "samplesheet used is ", samplesheet.file, "\n",
  "the sequence file used is - ",  fasta_run, "\n",
  "forward barcodes details are from - ", barcodeID.fwd, "\n",
  "reverse barcodes details are from - ", barcodeID.rev, "\n",
  "new fastq files will be in ", out_fastq, "\n",
  "the barcode details are in ", outputfile2, "\n",
  "\n"
  ))


#source("http://bioconductor.org/biocLite.R")
## codes2.samplesheet %>% filter(!is.na(demult_id))
#biocLite("Biostrings")
library(seqinr)


detectBatchCPUs <- function() {
  ncores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK"))
  if (is.na(ncores)) {
    ncores <- as.integer(Sys.getenv("SLURM_JOB_CPUS_PER_NODE"))
  }
  if (is.na(ncores)) {
    return(4) # for helix
  }
  return(ncores)
}

n.cores<-detectBatchCPUs()

no_w <- function(){Sys.time()}
message(paste0("Make_fasta_folders started on ", no_w(), ". \n check in ~30 mins."))

make_string <- function(x){return(toString(na.omit(unique(x))))}


samples_seqkit_fwd <- read_tsv(barcodeID.fwd)
samples_seqkit_rev <- read_tsv(barcodeID.rev)

Pacbio5_SampleSheet.orig <- read_excel(samplesheet.file) %>%
  arrange(Col_barcode, Row_barcode)
cat(paste0("samplesheet is ", samplesheet.file))
## glimpse(Pacbio5_SampleSheet.orig)


cat("\n  ## SHEET WITH READS WHERE TWO BARCODES WERE DETECTED ## \n")
seqkit_barcodes.all <- full_join(
  samples_seqkit_fwd %>% 
    select(seqID, Col_barcode = patternName, fwd_barcode_5to3 = matched) %>%
    mutate(fwd_barcode_5to3 = str_to_upper(fwd_barcode_5to3)),
  samples_seqkit_rev %>% 
    select(seqID, Row_barcode = patternName, rev_barcode_5to3 = matched) %>%
    mutate(rev_barcode_5to3 = str_to_upper(rev_barcode_5to3)),
  by="seqID") %>% 
  distinct() 

saveRDS(seqkit_barcodes.all, "seqkit_barcodes.all.RDS")
cat("\n written as seqkit_barcodes.all.RDS \n")

seqkit_barcodes <- seqkit_barcodes.all %>%
  filter(!is.na(Col_barcode)) %>%
  filter(!is.na(Row_barcode)) %>%
  arrange(Col_barcode, Row_barcode)


cat("\n")
cat("\n")
cat("### going to merge BarcodeID and samplesheet #### \n")
cat("\n")
cat("\n")
cat("seqkit_barcodes \n")
glimpse(seqkit_barcodes)
cat("Pacbio5_SampleSheet.orig \n")
glimpse(Pacbio5_SampleSheet.orig)

cat("### done merging BarcodeID and samplesheet #### \n")

cat("\n")
cat("\n")
cat("\n")

samplesheet_with_2codes <- left_join(seqkit_barcodes, Pacbio5_SampleSheet.orig, 
                                    by=c("Col_barcode", "Row_barcode", "fwd_barcode_5to3", "rev_barcode_5to3")) 

cat("samples with 2 barcodes. \n")
glimpse(samplesheet_with_2codes)
cat("\n written as samplesheet_with_both_codes.RDS or samplesheet_with_codes.tsv \n")
## saveRDS(samplesheet_with_2codes, "samplesheet_with_both_codes.RDS")
write_delim(samplesheet_with_2codes, file="samplesheet_with_codes.tsv", delim="\t")

codes2.samplesheet <- samplesheet_with_2codes  ## relay to next part

message("now binning...   \n")
############

suppressPackageStartupMessages(library(ShortRead))
if (str_detect(fasta_run, "fastq")) {
  sp <- SolexaPath(system.file('extdata', package='ShortRead'))
  pac2 <- readFastq(analysisPath(sp), pattern=basename(fasta_run), dirPath = dirname(fasta_run))
#  pac2 <- readQualityScaledDNAStringSet(fasta_run)
      frmt="fastq"
      system(paste0("mkdir ", out_fastq))
        } else if (str_detect(fasta_run, "fasta")) {
            pac2 <- readDNAStringSet(fasta_run, format="fasta")
            frmt = "fasta"
            system("mkdir ../fasta")
              } else if (!str_detect(fasta_run, "fasta|fastq")) {
                  stop("suffix fasta or fastq required on input file")
              }

message(paste0(fasta_run, " / pac2 has rows =", length(pac2), " -- ", no_w()))

sample_seq_relevant <- codes2.samplesheet %>%
  filter(!is.na(PCR)) 

samples_seq_relevant.nReads <- sample_seq_relevant %>% 
                                  group_by(PCR) %>%
                                    tally() %>%
                                      filter(n>10)
samples_seq_relevant.g10 <- left_join(sample_seq_relevant, samples_seq_relevant.nReads, by="PCR") %>%
                              filter(!is.na(n))

Samples_unique_g10 <- unique(samples_seq_relevant.g10$PCR) %>% sort() ## unique samples with reads greater than 10.
cat("samples = ", length(Samples_unique_g10))

registerDoMC(n.cores-1)
message(paste0("variable = n.cores has = ", n.cores, " cores registered. Use",(n.cores-1),  " cores.", no_w()))


#source("http://bioconductor.org/biocLite.R")
#biocLite("Biostrings")
library(Biostrings)
library(seqinr)
library(tidyverse)  
library(readxl)


### Progress bar
n=0
n_jobs <- length(Samples_unique_g10)
cat(paste0("samples = ", n_jobs))

pb <- txtProgressBar(min=0, initial=0, max = n_jobs, style = 3)
progress <- function(ja) setTxtProgressBar(pb, ja)
opts <- list(progress = progress)
message(paste0("going parallel with binning fasta", no_w()))

## sink file to track progress
writeLines(c(""), "log_Make_fastq.txt") 

foreach(i = Samples_unique_g10, .options.snow = opts) %dopar% ({
  ## track progress
  sink("log_Make_fastq.txt", append=TRUE) 
  cat(paste("Starting iteration ",i," of ",n_jobs, "rows.  At ", no_w(),"\n")) 
  sink()
  # trimws(i)
  suppressPackageStartupMessages(library(seqinr))
  suppressPackageStartupMessages(library(Biostrings))
  pac2.sample.a <- pac2[id(pac2) %in% subset(sample_seq_relevant, PCR==i)$seqID]
  file_name_JA <- paste0(out_fastq, "/PCR_", i, ".", frmt, collapse="", sep="")
  writeFastq(pac2.sample.a, file_name_JA, compress=FALSE)
  setTxtProgressBar(pb, i)
  invisible(gc())
})
