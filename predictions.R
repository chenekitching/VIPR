#!/usr/bin/env Rscript

# Load libraries ----------------------------------------------------------

library(dplyr)
library(readr)
library(modelr)
library(tidymodels)

# CommandArgs -------------------------------------------------------------
args = commandArgs(trailingOnly=TRUE)

# Read files --------------------------------------------------------------
mod <- readRDS(args[1])
input <- read_tsv(args[2]) %>%
  rowid_to_column("ID")

# Functions ---------------------------------------------------------------

get_annotations <- function(annovar_file, af_colname){
  
  #Remove all rows with LRG
  rm_lrg <- annovar_file %>%
    dplyr::filter(!grepl('LRG', Chr))
  
  #Select features needed for ML/EDA
  annovar_features <- rm_lrg %>%
    select(ID,
          Func.refGene,
           ExonicFunc.refGene, 
           af_colname,
           FATHMM_score,
           GenoCanyon_score,
           LRT_score,
           CADD_phred,
           `GERP++_RS`,
           `M-CAP_score`,
           MetaLR_score,
           MetaSVM_score, 
           MutationTaster_score,
           MutationAssessor_score,
           PROVEAN_score, 
           SIFT_score,
           SiPhy_29way_logOdds, 
           `fathmm-MKL_coding_score`,
           integrated_fitCons_score, 
           integrated_confidence_value, 
           phastCons100way_vertebrate,
           phastCons30way_mammalian,
           phyloP100way_vertebrate,
           phyloP30way_mammalian, 
           Otherinfo6) 
  
  
  #Replace "." with NA
  annovar_features[annovar_features == "."] <- NA
  annovar_features[annovar_features == ""] <- NA
  
  #Numeric features
  features_add <- annovar_features %>% 
    mutate_at(c(af_colname,
                'FATHMM_score',
                'GERP++_RS',
                'GenoCanyon_score',
                'LRT_score',
                "CADD_phred",
                'M-CAP_score',
                'MetaLR_score', 
                'MetaSVM_score', 
                'MutationTaster_score',
                'MutationAssessor_score',
                'PROVEAN_score', 
                'SIFT_score', 
                'SiPhy_29way_logOdds',
                'fathmm-MKL_coding_score',
                'integrated_fitCons_score', 
                'integrated_confidence_value',
                'phastCons100way_vertebrate',
                'phastCons30way_mammalian', 
                'phyloP100way_vertebrate',
                'phyloP30way_mammalian'),
              as.numeric)
  

}
rm_missing_rows <- function(df){
  
  nas <- df %>% 
    mutate(n_NAs = rowSums(is.na(.)))
  rm_missing <- subset(nas, n_NAs <= 18) %>%
    select(-n_NAs) %>%
    mutate(AF_NA = ifelse(is.na(AF), 1, 0)) %>%
    mutate(AF = replace_na(AF, 0)) 
  
}


# Preprocessing -----------------------------------------------------------

input_feat <- get_annotations(input, "AF") %>%
  rm_missing_rows() %>%
  rename(`GERP.._RS` = `GERP++_RS`,
         `M.CAP_score` = `M-CAP_score`,
         `fathmm.MKL_coding_score` = `fathmm-MKL_coding_score`)


# Predict -----------------------------------------------------------------

preds <- add_predictions(input_feat,
                         mod,
                         type = "prob") %>%
  unnest(pred) %>%
  select(ID, `.pred_P_LP`)

final_preds <- input %>%
  merge(preds,
        by = "ID") %>%
  arrange(desc(`.pred_P_LP`)) %>%
  select(-ID)

write_tsv(final_preds,
  file = args[3])