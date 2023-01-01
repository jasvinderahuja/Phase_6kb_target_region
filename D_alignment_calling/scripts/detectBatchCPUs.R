#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(foreach))
suppressPackageStartupMessages(library(doMC))
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
registerDoMC(n.cores-1)
message(paste0("variable = n.cores has = ", n.cores, " cores registered. Use",(n.cores-1),  " cores"))
### example that works
# library(doParallel)
# source("../scripts/detectBatchCPUs.R")
# source("/data/ahujajs/scripts/detectBatchCPUs.R")
# feach <- foreach(i = 1:nrow(positionsTOcorrect), .combine = rbind) %dopar% {
#   require(tidyverse)
#   require(stringr)
#   pos <-  if (positionsTOcorrect$checkF[i]==positionsTOcorrect$checkTo[i]) {
#     positionsTOcorrect$checkF[i]
#   } else {
#     (positionsTOcorrect$checkF[i]:positionsTOcorrect$checkTo[i])
#   }
#   Calls.parts <- Calls.all.parts %>% 
#     filter(snp_start %in% pos) %>% 
#     rowwise() %>% 
#     mutate (nt_det1 = ifelse (str_detect(seq, positionsTOcorrect$checkS[i]), positionsTOcorrect$S[i], NA), 
#             nt_det2 = ifelse (str_detect(seq, positionsTOcorrect$checkW[i]), positionsTOcorrect$W[i], NA),
#             snp_start = positionsTOcorrect$from[i]) %>%
#     distinct()
# }