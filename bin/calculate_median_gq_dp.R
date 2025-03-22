# Set env. ------------------------------------------------
library(tidyverse)
library(dplyr)
library(tidylog)
library(data.table)

# Read input ---------------------------------------------- [ Read arguments ]

input_file <- read_tsv("aux/raw_sample_variant_gq_dp.tsv",
                       col_types = c("cccccccc"),
                       col_names = c("chr","pos","ref","alt","sample","gq","dp")) %>% 
  mutate(variant_id = paste(chr,pos,ref,alt, sep = "_")) 

# Retrieve and calculate metrics -------------------------- [ By variant ]

stats_obj <- input_file %>% 
  mutate(gq = case_when(gq == "." ~ NA, TRUE ~ gq),
         dp = case_when(dp == "." ~ NA, TRUE ~ dp)) %>% 
  mutate(gq = as.numeric(gq),
         dp = as.numeric(dp)) %>% 
  group_by(variant_id) %>% 
  summarise(medianGQ = median(gq, na.rm = T),
            medianDP = median(dp, na.rm = T)) %>% 
  ungroup() %>% 
  mutate(medianGQ = case_when(is.na(medianGQ) ~ 0, TRUE ~ medianGQ),
         medianDP = case_when(is.na(medianDP) ~ 0, TRUE ~ medianDP))

out_file <- input_file %>% 
  distinct(chr, pos, ref, alt, variant_id) %>% 
  rename(`#CHROM` = chr,
         POS = pos,
         REF = ref,
         ALT = alt) %>% 
  left_join(stats_obj, 
            by = "variant_id")

# Write output -------------------------------------------- 

out_file %>% 
  select(-variant_id) %>% 
  write_tsv("aux/gq_dp_annotation_file.txt", col_names = T)

# Answers to exercise questions --------------------------

input_file %>% 
  mutate(gq = case_when(gq == "." ~ NA, TRUE ~ gq),
         dp = case_when(dp == "." ~ NA, TRUE ~ dp)) %>% 
  mutate(gq = as.numeric(gq),
         dp = as.numeric(dp)) %>% 
  summarise(minGQ = min(gq, na.rm = T),
            maxGQ = max(gq, na.rm = T),
            minDP = min(dp, na.rm = T),
            maxDP = max(dp, na.rm = T))

out_file %>% 
  summarise(minMedGQ = min(medianGQ, na.rm = T),
            maxMedGQ = max(medianGQ, na.rm = T),
            minMedDP = min(medianDP, na.rm = T),
            maxMedDP = max(medianDP, na.rm = T))


