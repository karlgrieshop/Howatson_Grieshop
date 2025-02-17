#__________________________----
#PACKAGES ----
library(tidyverse) # Tidy data packages
library(janitor) # Clean data names
library(dplyr) #Data manipulation.
library(stringr) # Simplifies string operations
library(tidyr) # Tidying Data
library(biomaRt)
#__________________________----
#IMPORT DATA ----
setwd("/Users/elliehowatson/Desktop/Uni/Year 3/Final Year Project /R/data")

flybase_orthologs_PCOS <- read_csv("/Users/elliehowatson/Desktop/Uni/Year 3/Final Year Project /R/data/diopt_results_PCOS_genes.csv")
flybase_human_SA_orthologs <- read_csv("/Users/elliehowatson/Desktop/Uni/Year 3/Final Year Project /R/data/diopt_results.csv")

disease_files <- list.files(path = "/Users/elliehowatson/Desktop/Uni/Year 3/Final Year Project /R/data", pattern = "*.txt", full.names = TRUE)

disease_data <- lapply(disease_files, function(file) {
  disease_genes <- read_lines(file)
  disease_name <- tools::file_path_sans_ext(basename(file))  # Extract disease name from the file name
  data.frame(FlyBase_gene_ID = disease_genes, Disease = disease_name)
})

# Combine all disease data into one dataframe
all_diseases_data <- bind_rows(disease_data)

#__________________________----
#CLEAN_DATA ----


df <- flybase_human_SA_orthologs %>% select(-input_order)
