#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(zoo))
suppressPackageStartupMessages(library(tidyverse))


stds <- readRDS("PBStandard.RDS")
PB1_Standard <- stds %>% 
                  arrange(Parent, REF_POS) %>% 
                  select(Parent, REF_POS, BASE) %>% 
                  group_by(Parent) %>% 
                  spread(REF_POS, BASE)

dtime <- function(){now <- Sys.time()}
#message(paste0("started on ", now))
#
read.sam2tsv <- function(arg){
          a <- read_tsv(arg, 
              col_names = TRUE, 
              col_types=cols(
              `#READ_NAME` = col_character(),
              FLAG = col_double(),
              CHROM = col_character(),
              READ_POS = col_double(),
              BASE = col_character(),
              QUAL = col_character(),
              REF_POS = col_double(),
              REF = col_character(),
              OP = col_character()
            ),
            na = c("", "NA", ".")
      )
    return(a)
}



# Red = S4955
# Blue = S5085

# S4955 <- "../PB1/sam2tsv/S4955.sam2tsv.out"
# S4955 <- read.sam2tsv(S4955)
# saveRDS(S4955, "S4955.RDS")
# S4955xS5085 <- "../PB1/sam2tsv/S4955xS5085.sam2tsv.out"
# S4955xS5085 <- read.sam2tsv(S4955xS5085)
# saveRDS(S4955xS5085, "S4955xS5085.RDS")
# S5085 <- "../PB1/sam2tsv/S5085.sam2tsv.out"
# S5085 <- read.sam2tsv(S5085)
# saveRDS(S5085, "S5085.RDS")


positions <- c(2580:2593, 2968, 2969, 2970, 
               3180, 3230:3233, 3807:3814, 3847, 4176, 4185, 
               4497, 4593, 4656, 4719, 4788, 4889, 5103, 5205, 5318, 5417, 5529, 5630, 5725, 
               5826, 5835, 5941, 6041, 6149, 6281, 6368, 6455, 6531, 6533, 6617, 6734, 6842, 
               6914, 7040, 7127, 7199, 7271, 7346, 7401, 7404:7406, 8078, 8090:8102,
               8200:8207, 8592)
pos_groups <- c(rep(2585, 14), 2968, 2969, 2970, 
                3180, rep(3230, 4), rep(3808, 8), 3847, 4176, 4185, 
                4497, 4593, 4656, 4719, 4788, 4889, 5103, 5205, 5318, 5417, 5529, 5630, 5725, 
                5826, 5835, 5941, 6041, 6149, 6281, 6368, 6455, 6531, 6533, 6617, 6734, 6842, 
                6914, 7040, 7127, 7199, 7271, 7346, 7401, rep(7405,3), 8078, rep(8095, 13), 
                rep(8202, 8), 8592)
saveRDS(read_pos, "pos_groups.RDS")

## could not confirm 8216:8221, 8223:8229, 8559:8564 or rep(8219, 6), rep(8227,7), rep(8561, 6) as valid positions

read_pos <- list(REF_POS=positions, score_pos=pos_groups) %>% as_tibble()

read_RDS <- function(file_RDS="S4955.RDS"){
      ## extracting the file_id and messages for troubleshooting
        file_idn <- str_extract(file_RDS, "pre(.+)sam2tsv") %>% 
                    str_replace_all("sam2tsv|pre|PCR_|\\.", "") %>%
                      as.numeric()
        message(paste0(file_RDS))
        message(dtime())
        message(paste0(" file_id= ", file_idn))
        pre.RDS <- readRDS(file_RDS) %>%
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
          summarise(BASES = paste0(BASE, sep="", collapse="")) %>%
      ## tally and get frequency of base occurance at all positions
          group_by(score_pos, BASES) %>% 
          tally() %>% 
          mutate(freq=n/sum(n)) %>% 
          filter(freq>0.1) %>% 
       ## put the file_id
          mutate(file_id = file_idn)
      
        return(pre.RDS)
        gc()
}

## read in the sam2tsv files I saved as RDS for ~ quicker loading
blue.pres <- list.files(".", "Blue.pre")
red.pres <- list.files(".", "Red.pre")

blues <- bind_rows(lapply(blue.pres, read_RDS))
reds <- bind_rows(lapply(red.pres, read_RDS))

message("#1")

blues_gp <- blues %>% 
              select(score_pos, file_id, BASES) %>%
              mutate(file_id = paste0("blue_", file_id)) %>%
              group_by(score_pos, file_id) %>%
              summarize(BASES = paste0(BASES, collapse="_")) %>%
              group_by(score_pos) %>%
              spread(file_id, BASES)
message("#2")
reds_gp <- reds %>% 
              select(score_pos, file_id, BASES) %>%
              mutate(file_id = paste0("red_", file_id)) %>%
              group_by(score_pos, file_id) %>%
              summarize(BASES = paste0(BASES, collapse="_")) %>%
              group_by(score_pos) %>%
              spread(file_id, BASES)
message("#3")
red_blue <- left_join(reds_gp, blues_gp, by = "score_pos")
alarm()
########### Now combine the reds and blues and create a standard #######

read_together <- function(file_RDS="S4955.RDS"){

  message(file_RDS)
  message(dtime())
  dpre.RDS <- readRDS(file_RDS) %>%
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
    summarise(BASES = paste0(BASE, sep="", collapse="")) %>%
    ## tally and get frequency of base occurance at all positions
    group_by(score_pos, BASES) %>% 
    tally() %>% 
    mutate(freq=n/sum(n)) %>% 
    filter(freq>0.1)
  return(dpre.RDS)
  gc()
}

alarm()
gc()
bind_rows(lapply(blue.pres, readRDS)) %>% saveRDS("blue.pre.sam2tsv.RDS")
gc()
message("#4")
bind_rows(lapply(red.pres, readRDS)) %>% saveRDS("red.pre.sam2tsv.RDS")
gc()
message("#5")
blues_together <- read_together("blue.pre.sam2tsv.RDS")
message("#6")
reds_together <- read_together("red.pre.sam2tsv.RDS")
message("#7")

saveRDS(list(blues = blues_together, reds = reds_together), "blue_red_allDATA.RDS")

blue_bases <- blues_together %>% 
  select(score_pos, BASES) %>%
  group_by(score_pos) %>%
  summarize(bases = paste0(unique(BASES), collapse="_")) %>%
  mutate(Parent = "W")
message("#8")


red_bases <- reds_together %>% 
  select(score_pos, BASES) %>%
  group_by(score_pos) %>%
  summarize(bases = paste0(unique(BASES), collapse="_")) %>%
  mutate(Parent = "S")
message("#9")


compare_SW <- left_join(blue_bases %>% select(score_pos, bases), red_bases %>% select(score_pos, bases), by = "score_pos")
message("#10")


standards_PB2 <- bind_rows(blue_bases, red_bases)
saveRDS(standards_PB2, "standards_PB2.RDS")
message("#11")

standards.results <-  list(positions_scored = read_pos,
                            compare_parents_individually = red_blue,
                            compare_parents = compare_SW,
                            final_standard = standards_PB2,
                            old_PB1_standard = PB1_Standard
                            )
saveRDS(standards.results, "standards.results.RDS")
suppressPackageStartupMessages(library(openxlsx))
write.xlsx(standards.results, "standards.results.xlsx")
message("DONE")
alarm()
# standards <- bind_rows(blues %>% mutate(P="W"), reds %>% mutate(P="S"))
# saveRDS(standards, "PBStandard_need_work.RDS")
## saveRDS(standards, "PBStandard.RDS")