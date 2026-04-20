#PACKAGES ----
library(tidyverse) # Tidy data packages
library(janitor) # Clean data names
library(dplyr) #Data manipulation.
library(stringr) # Simplifies string operations
library(tidyr) # Tidying Data
library(biomaRt) # BioMart databases like Ensembl to query genomic data
library(ggplot2) # Data visualisation
library(clipr) # Copying data between R and other applications
library(broom) # Statistial analyses
library(patchwork)
#__________________________----
#IMPORTING DATA (FlyBase GeneIDs) ----

setwd("/Users/elliehowatson/Final_Year_Project/R/data/Sex specific diseases")

# Importing all disease genes from Flybase----

all_disease_genes <-read.csv("/Users/elliehowatson/Final_Year_Project/R/data/all_disease_genes.csv")
colnames(all_disease_genes)[1] <- "FlyBaseID"

# Importing all human ortholog genes----

all_orthologs_genes <-read.csv("/Users/elliehowatson/Final_Year_Project/R/data/all_orthologs.csv")

# Imported dataset for Mishra et al. 2024 ----

# SBGE contrasts

SBG_whole_bodies <-read.csv("/Users/elliehowatson/Final_Year_Project/R/data/Sex.Biased.Genes.Whole.csv")
colnames(SBG_whole_bodies)[1] <- "FlyBaseID"

# Imported dataset from Fraisse et al. 2018 File S7 -----

Fra <- read.delim(file="/Users/elliehowatson/Final_Year_Project/R/data/FileS7.txt", header=TRUE)
colnames(Fra)[2] <- "FlyBaseID"
Fra[sapply(Fra, is.infinite)] <- NA

# Imported dataset from Singh and Agrawal 2023 ----

Singh <- read.csv(file="/Users/elliehowatson/Final_Year_Project/R/data/SBGEandSSSdataForMBE.csv", header=TRUE)
colnames(Singh)[2] <- "FlyBaseID"

# Imported dataset my disease dataset that merged Mishra et al. 2024

SBGE_disease_data <- read.csv("/Users/elliehowatson/Final_Year_Project/R/data/SBGE_disease_data.csv")
SBGE_disease_data$Tissue <- NULL
SBGE_disease_data$GeneType<- "Candidate"
# Imported dataset my disease dataset that merged Mishra et al. 2024 and Fraisse et al. 2018

mish_fra_SBGE<- read.csv("/Users/elliehowatson/Final_Year_Project/R/data/mish_fra_SBGE.csv")
mish_fra_SBGE$Tissue <- NULL
mish_fra_SBGE$Fly.Symbol <- NULL
mish_fra_SBGE$GeneType<- "Candidate"


#__________________________----
#CLEAN DATA

# Cleaning all orthologs df ----
# Keep only the FlyBase gene ID column
all_orthologs_list<- all_orthologs_genes[ , "X..Dmel_gene_ID", drop = FALSE]

# Remove duplicates
all_orthologs_list<- unique(all_orthologs_list)

# Rename the column
names(all_orthologs_list) <- "FlyBaseID"

# Cleaning all disease df ----
# Keep only the FlyBase gene ID column
all_disease_genes_list<- all_disease_genes[ , "FlyBaseID", drop = FALSE]

# Remove duplicates
all_disease_genes_list<- unique(all_disease_genes_list)

# Rename the column
names(all_disease_genes_list) <- "FlyBaseID"
#__________________________----
# DATA ANALYSIS ----

orthologs_expr<- all_orthologs_list %>%
  inner_join(SBG_whole_bodies, by = "FlyBaseID")

names(orthologs_expr)[names(orthologs_expr) == "log2FoldChange"] <- "exp_SBGE_ase"

orthologs_Fra<- orthologs_expr %>%
  inner_join(Fra, by = "FlyBaseID")

orthologs_Singh<- orthologs_expr %>%
  inner_join(Singh, by = "FlyBaseID")

# Orthologs ----
# Fraisse -----
# 1) Add the type labels

orthologs_Fra$GeneType<- "Ortholog"

orthologs_Fra<- orthologs_Fra[, c(1, 3, 2, 4:ncol(orthologs_Fra))]

# 2) Remove orthologues that are also in the candidate disease set
orthologs_Fra <- orthologs_Fra[ !(orthologs_Fra$FlyBaseID %in% mish_fra_SBGE$FlyBaseID), ]

# 3) Bind together
candidate_ortholog_Fra <- rbind(mish_fra_SBGE, orthologs_Fra)
candidate_ortholog_Fra<- candidate_ortholog_Fra[!duplicated(candidate_ortholog_Fra), ]

# Singh -----

# 1) Add the gene type labels

orthologs_Singh$GeneType<- "Ortholog"

orthologs_Singh<- orthologs_Singh[, c(1, 3, 2, 4:ncol(orthologs_Singh))]

# 2) Remove orthologues that are also in the disease set

mish_singh_disease<- SBGE_disease_data %>%
  inner_join(Singh, by = "FlyBaseID")

mish_singh_disease<- mish_singh_disease[, c(3, 1, 2, 5:ncol(mish_singh_disease), 4)]

orthologs_Singh <- orthologs_Singh[ !(orthologs_Singh$FlyBaseID %in% mish_singh_disease$FlyBaseID), ]

# 3) Bind together 

candidate_ortholog_Singh <- rbind(mish_singh_disease, orthologs_Singh)


#All disease genes----

all_disease_genes_expr<- all_disease_genes_list %>%
  inner_join(SBG_whole_bodies, by = "FlyBaseID")

names(all_disease_genes_expr)[names(all_disease_genes_expr) == "log2FoldChange"] <- "exp_SBGE_ase"

all_disease_genes_Fra<- all_disease_genes_expr %>%
  inner_join(Fra, by = "FlyBaseID")

all_disease_genes_Singh<- all_disease_genes_expr %>%
  inner_join(Singh, by = "FlyBaseID")


# Fraisse -----
# 1) Add the gene type labels

all_disease_genes_Fra$GeneType<- "Disease"

all_disease_genes_Fra<- all_disease_genes_Fra[, c(1, 3, 2, 4:ncol(all_disease_genes_Fra))]

# 2) Remove orthologues that are also in the disease set
all_disease_genes_Fra <- all_disease_genes_Fra[ !(all_disease_genes_Fra$FlyBaseID %in% mish_fra_SBGE$FlyBaseID), ]

# 3) Bind together
candidate_all_disease_Fra <- rbind(mish_fra_SBGE, all_disease_genes_Fra )
candidate_all_disease_Fra <- candidate_all_disease_Fra [!duplicated(candidate_all_disease_Fra ), ]

# Singh -----

# 1) Add the gene type labels

all_disease_genes_Singh$GeneType<- "Disease"

all_disease_genes_Singh<- all_disease_genes_Singh[, c(1, 3, 2, 4:ncol(all_disease_genes_Singh))]

# 2) Remove orthologues that are also in the disease set

mish_singh_disease<- SBGE_disease_data %>%
  inner_join(Singh, by = "FlyBaseID")

mish_singh_disease<- mish_singh_disease[, c(3, 1, 2, 5:ncol(mish_singh_disease), 4)]

all_disease_genes_Singh <- all_disease_genes_Singh[ !(all_disease_genes_Singh$FlyBaseID %in% mish_singh_disease$FlyBaseID), ]

# 3) Bind together

candidate_all_disease_Singh <- rbind(mish_singh_disease, all_disease_genes_Singh)


#ORTHOLOGS----

# Tajimas D Summary table ----

SBGE_tajDS_compare_orth<- candidate_ortholog_Singh%>%
  select(FlyBaseID,Sex.Bias.Category,exp_SBGE_ase,tajD.S,GeneType)%>% 
  distinct() # Remove disease variable to avoid repeats

SBGE_tajDS_compare_orth<- SBGE_tajDS_compare_orth %>%
  mutate(exp_SBGE_ase = as.numeric(exp_SBGE_ase))

# Creating Bins

SBGE_tajDS_bins_compare_orth<-SBGE_tajDS_compare_orth %>%
  mutate(Sex.Bias.Category = case_when(
    exp_SBGE_ase >= 2 ~ "Highly male-biased",
    exp_SBGE_ase >= 0.5 & exp_SBGE_ase < 2 ~ "Moderately male-biased",
    exp_SBGE_ase >= -0.5 & exp_SBGE_ase < 0.5 ~ "Unbiased",
    exp_SBGE_ase >= -2 & exp_SBGE_ase < -0.5 ~ "Moderately female-biased",
    exp_SBGE_ase < -2 ~ "Highly female-biased"
  ))

SBGE_tajDS_bins_compare_orth <- SBGE_tajDS_bins_compare_orth[!is.na(SBGE_tajDS_bins_compare_orth$tajD.S), ]

# Compute Mean and Standard Error (SE) 

# Calculate mean
mean_tajD_orth <- aggregate(tajD.S ~ Sex.Bias.Category + GeneType, 
                       data = SBGE_tajDS_bins_compare_orth, 
                       FUN = function(x) mean(x, na.rm = TRUE))

# Calculate standard error
se_tajD_orth <- aggregate(tajD.S ~ Sex.Bias.Category + GeneType, 
                     data = SBGE_tajDS_bins_compare_orth, 
                     FUN = function(x) sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x))))

# Combine into one data frame
summary_tajDS_data_comparative_orth <- merge(mean_tajD_orth, se_tajD_orth, 
                                        by = c("Sex.Bias.Category", "GeneType"))

# Rename columns
colnames(summary_tajDS_data_comparative_orth) <- c("Sex.Bias.Category", "GeneType", 
                                              "Mean_tajD.S", "SE_tajD.S")

sample_sizes_tadj_orth <- as.data.frame(table(SBGE_tajDS_bins_compare_orth$Sex.Bias.Category,
                                         SBGE_tajDS_bins_compare_orth$GeneType))

# Rename columns
colnames(sample_sizes_tadj_orth) <- c("Sex.Bias.Category", "GeneType", "SampleSize")

summary_tajDS_data_comparative_orth <- left_join(summary_tajDS_data_comparative_orth, sample_sizes_tadj_orth, by = c("Sex.Bias.Category", "GeneType"))

# Then compute the "All" category stats

split_data<- split(summary_tajDS_data_comparative_orth, summary_tajDS_data_comparative_orth$GeneType)

all_category <- do.call(rbind, lapply(split_data, function(df) {
  data.frame(
    GeneType = unique(df$GeneType),
    Mean_tajD.S = mean(df$Mean_tajD.S, na.rm = TRUE),
    SE_tajD.S = sqrt(sum(df$SE_tajD.S^2, na.rm = TRUE) / nrow(df)),
    SampleSize = sum(df$SampleSize, na.rm = TRUE),
    Sex.Bias.Category = "All"
  )
}))

summary_tajDS_data_comparative_orth <- bind_rows(summary_tajDS_data_comparative_orth, all_category)

# rMF Summary table ----

# Change into numerical 

SBGE_rMF_compare_orth<- candidate_ortholog_Singh%>%
  select(FlyBaseID,Sex.Bias.Category,exp_SBGE_ase,rmf,GeneType)%>% 
  distinct() # Remove disease variable to avoid repeats

SBGE_rMF_compare_orth<- SBGE_rMF_compare_orth%>%
  mutate(exp_SBGE_ase = as.numeric(exp_SBGE_ase))

# Creating Bins

SBGE_rMF_bins_compare_orth<-SBGE_rMF_compare_orth%>%
  mutate(Sex.Bias.Category = case_when(
    exp_SBGE_ase >= 2 ~ "Highly male-biased",
    exp_SBGE_ase >= 0.5 & exp_SBGE_ase < 2 ~ "Moderately male-biased",
    exp_SBGE_ase >= -0.5 & exp_SBGE_ase < 0.5 ~ "Unbiased",
    exp_SBGE_ase >= -2 & exp_SBGE_ase < -0.5 ~ "Moderately female-biased",
    exp_SBGE_ase < -2 ~ "Highly female-biased"
  ))

SBGE_rMF_bins_compare_orth <- SBGE_rMF_bins_compare_orth[!is.na(SBGE_rMF_bins_compare_orth$rmf), ]

# Compute Mean and Standard Error (SE) 

# Compute mean
mean_rmf_orth <- aggregate(rmf ~ Sex.Bias.Category + GeneType, 
                      data = SBGE_rMF_bins_compare_orth, 
                      FUN = function(x) mean(x, na.rm = TRUE))

# Compute standard error
se_rmf_orth <- aggregate(rmf ~ Sex.Bias.Category + GeneType, 
                    data = SBGE_rMF_bins_compare_orth, 
                    FUN = function(x) sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x))))

# Merge the two summaries
summary_rMF_data_comparative_orth <- merge(mean_rmf_orth, se_rmf_orth, 
                                      by = c("Sex.Bias.Category", "GeneType"))

# Rename columns
colnames(summary_rMF_data_comparative_orth) <- c("Sex.Bias.Category", "GeneType", 
                                            "Mean_rmf", "SE_rmf")

rMF_sample_sizes_orth <- as.data.frame(table(SBGE_rMF_bins_compare_orth$Sex.Bias.Category,
                                        SBGE_rMF_bins_compare_orth$GeneType))

# Rename columns 
colnames(rMF_sample_sizes_orth) <- c("Sex.Bias.Category", "GeneType", "SampleSize")
summary_rMF_data_comparative_orth <- left_join(summary_rMF_data_comparative_orth, rMF_sample_sizes_orth, by = c("Sex.Bias.Category", "GeneType"))

# Split by GeneType
split_data_rmf_orth <- split(summary_rMF_data_comparative_orth, summary_rMF_data_comparative_orth$GeneType)

# For each GeneType, compute "All" summary row
all_category_rmf_orth <- do.call(rbind, lapply(split_data_rmf_orth , function(df) {
  data.frame(
    GeneType = unique(df$GeneType),
    Mean_rmf = mean(df$Mean_rmf, na.rm = TRUE),
    SE_rmf = sqrt(sum(df$SE_rmf^2, na.rm = TRUE) / nrow(df)),
    SampleSize = sum(df$SampleSize, na.rm = TRUE),
    Sex.Bias.Category = "All"
  )
}))

summary_rMF_data_comparative_orth <- bind_rows(summary_rMF_data_comparative_orth,all_category_rmf_orth)


# DoS Summary table ----
# Change into numerical 

SBGE_DoS_compare_orth<-candidate_ortholog_Singh %>%
  select(FlyBaseID,Sex.Bias.Category,exp_SBGE_ase,DoS,GeneType)%>% 
  distinct() # Remove disease variable to avoid repeats


SBGE_DoS_compare_orth<- SBGE_DoS_compare_orth%>%
  mutate(exp_SBGE_ase = as.numeric(exp_SBGE_ase))

SBGE_DoS_compare_orth <- SBGE_DoS_compare_orth[!is.na(SBGE_DoS_compare_orth$DoS), ]

# Creating Bins

SBGE_DoS_bins_compare_orth<-SBGE_DoS_compare_orth%>%
  mutate(Sex.Bias.Category = case_when(
    exp_SBGE_ase >= 2 ~ "Highly male-biased",
    exp_SBGE_ase >= 0.5 & exp_SBGE_ase < 2 ~ "Moderately male-biased",
    exp_SBGE_ase >= -0.5 & exp_SBGE_ase < 0.5 ~ "Unbiased",
    exp_SBGE_ase >= -2 & exp_SBGE_ase < -0.5 ~ "Moderately female-biased",
    exp_SBGE_ase < -2 ~ "Highly female-biased"
  ))

# Compute Mean and Standard Error (SE) 

# Calculate mean
mean_DoS_orth <- aggregate(DoS ~ Sex.Bias.Category + GeneType, 
                      data = SBGE_DoS_bins_compare_orth, 
                      FUN = function(x) mean(x, na.rm = TRUE))

# Calculate standard error
se_DoS_orth <- aggregate(DoS ~ Sex.Bias.Category + GeneType, 
                    data = SBGE_DoS_bins_compare_orth, 
                    FUN = function(x) sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x))))

# Merge both summaries
summary_DoS_data_comparative_orth <- merge(mean_DoS_orth, se_DoS_orth,
                                      by = c("Sex.Bias.Category", "GeneType"))

# Rename columns
colnames(summary_DoS_data_comparative_orth) <- c("Sex.Bias.Category", "GeneType", "Mean_DoS", "SE_DoS")


DoS_sample_sizes_orth <- as.data.frame(table(SBGE_DoS_bins_compare_orth$Sex.Bias.Category,
                                        SBGE_DoS_bins_compare_orth$GeneType))

# Rename columns
colnames(DoS_sample_sizes_orth) <- c("Sex.Bias.Category", "GeneType", "SampleSize")

summary_DoS_data_comparative_orth <- left_join(summary_DoS_data_comparative_orth, DoS_sample_sizes_orth, by = c("Sex.Bias.Category", "GeneType"))


# Split by GeneType
split_data_DoS_orth <- split(summary_DoS_data_comparative_orth, summary_DoS_data_comparative_orth$GeneType)

# For each GeneType, compute "All" summary row
all_category_DoS_orth <- do.call(rbind, lapply(split_data_DoS_orth , function(df) {
  data.frame(
    GeneType = unique(df$GeneType),
    Mean_DoS = mean(df$Mean_DoS, na.rm = TRUE),
    SE_DoS = sqrt(sum(df$SE_DoS^2, na.rm = TRUE) / nrow(df)),
    SampleSize = sum(df$SampleSize, na.rm = TRUE),
    Sex.Bias.Category = "All"
  )
}))

# Reorder columns to match the expected format
all_category_DoS_orth<- all_category_DoS_orth [, c("Sex.Bias.Category", "GeneType", "Mean_DoS", "SE_DoS", "SampleSize")]


summary_DoS_data_comparative_orth <- bind_rows(summary_DoS_data_comparative_orth,all_category_DoS_orth)


# Tissue Specificity Summary table  ----

# Change into numerical 

SBGE_tissuespec_compare_orth<- candidate_ortholog_Fra%>%
  select(FlyBaseID,Sex.Bias.Category,exp_SBGE_ase,tissueSpecificityAdults,GeneType)%>% 
  distinct() # Remove disease variable to avoid repeats

SBGE_tissuespec_compare_orth<- SBGE_tissuespec_compare_orth%>%
  mutate(exp_SBGE_ase = as.numeric(exp_SBGE_ase))

SBGE_tissuespec_compare_orth <- SBGE_tissuespec_compare_orth[!is.na(SBGE_tissuespec_compare_orth$tissueSpecificityAdults), ]


# Creating Bins

SBGE_tissuespec_bins_compare_orth<-SBGE_tissuespec_compare_orth%>%
  mutate(Sex.Bias.Category = case_when(
    exp_SBGE_ase >= 2 ~ "Highly male-biased",
    exp_SBGE_ase >= 0.5 & exp_SBGE_ase < 2 ~ "Moderately male-biased",
    exp_SBGE_ase >= -0.5 & exp_SBGE_ase < 0.5 ~ "Unbiased",
    exp_SBGE_ase >= -2 & exp_SBGE_ase < -0.5 ~ "Moderately female-biased",
    exp_SBGE_ase < -2 ~ "Highly female-biased"
  ))

SBGE_tissuespec_bins_compare_orth <- SBGE_tissuespec_bins_compare_orth[!is.na(SBGE_tissuespec_bins_compare_orth$tissueSpecificityAdults), ]

# Compute Mean and Standard Error (SE) 

# Calculate mean
mean_tissueSpecificityAdults_orth<- aggregate(tissueSpecificityAdults~ Sex.Bias.Category + GeneType, 
                                         data = SBGE_tissuespec_bins_compare_orth, 
                                         FUN = function(x) mean(x, na.rm = TRUE))

# Calculate standard error
se_tissueSpecificityAdults_orth<- aggregate(tissueSpecificityAdults ~ Sex.Bias.Category + GeneType, 
                                       data = SBGE_tissuespec_bins_compare_orth, 
                                       FUN = function(x) sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x))))

# Merge both summaries
summary_tissuespec_data_comparative_orth <- merge(mean_tissueSpecificityAdults_orth,se_tissueSpecificityAdults_orth,
                                             by = c("Sex.Bias.Category", "GeneType"))

# Rename columns
colnames(summary_tissuespec_data_comparative_orth) <- c("Sex.Bias.Category", "GeneType", "Mean_tissueSpecificityAdults", "SE_tissueSpecificityAdults")


tissuespec_sample_sizes_orth <- as.data.frame(table(SBGE_tissuespec_bins_compare_orth$Sex.Bias.Category,
                                               SBGE_tissuespec_bins_compare_orth$GeneType))

# Rename columns
colnames(tissuespec_sample_sizes_orth) <- c("Sex.Bias.Category", "GeneType", "SampleSize")

summary_tissuespec_data_comparative_orth <- left_join(summary_tissuespec_data_comparative_orth, tissuespec_sample_sizes_orth, by = c("Sex.Bias.Category", "GeneType"))


# Split by GeneType
split_data_tissuespec_orth <- split(summary_tissuespec_data_comparative_orth, summary_tissuespec_data_comparative_orth$GeneType)

# For each GeneType, compute "All" summary row
all_category_tissuespec_orth <- do.call(rbind, lapply(split_data_tissuespec_orth , function(df) {
  data.frame(
    GeneType = unique(df$GeneType),
    Mean_tissueSpecificityAdults = mean(df$Mean_tissueSpecificityAdults, na.rm = TRUE),
    SE_tissueSpecificityAdults  = sqrt(sum(df$SE_tissueSpecificityAdults^2, na.rm = TRUE) / nrow(df)),
    SampleSize = sum(df$SampleSize, na.rm = TRUE),
    Sex.Bias.Category = "All"
  )
}))

# Reorder columns to match the expected format
all_category_tissuespec_orth <- all_category_tissuespec_orth[, c("Sex.Bias.Category", "GeneType", "Mean_tissueSpecificityAdults", "SE_tissueSpecificityAdults", "SampleSize")]


summary_tissuespec_data_comparative_orth <- bind_rows(summary_tissuespec_data_comparative_orth,all_category_tissuespec_orth)

# Relative Level of Nonsynonymous Diversity Summary table  ----

# Change into numerical 

SBGE_rlnd_compare_orth<- candidate_ortholog_Singh%>%
  select(FlyBaseID,Sex.Bias.Category,exp_SBGE_ase,piNfract,GeneType)%>% 
  distinct() # Remove disease variable to avoid repeats

SBGE_rlnd_compare_orth<- SBGE_rlnd_compare_orth%>%
  mutate(exp_SBGE_ase = as.numeric(exp_SBGE_ase))

# Creating Bins

SBGE_rlnd_bins_compare_orth<-SBGE_rlnd_compare_orth%>%
  mutate(Sex.Bias.Category = case_when(
    exp_SBGE_ase >= 2 ~ "Highly male-biased",
    exp_SBGE_ase >= 0.5 & exp_SBGE_ase < 2 ~ "Moderately male-biased",
    exp_SBGE_ase >= -0.5 & exp_SBGE_ase < 0.5 ~ "Unbiased",
    exp_SBGE_ase >= -2 & exp_SBGE_ase < -0.5 ~ "Moderately female-biased",
    exp_SBGE_ase < -2 ~ "Highly female-biased"
  ))

SBGE_rlnd_bins_compare_orth <- SBGE_rlnd_bins_compare_orth[!is.na(SBGE_rlnd_bins_compare_orth$piNfract), ]

# Compute Mean and Standard Error (SE) 

# Calculate mean
mean_rlnd_orth<- aggregate(piNfract~ Sex.Bias.Category + GeneType, 
                      data = SBGE_rlnd_bins_compare_orth, 
                      FUN = function(x) mean(x, na.rm = TRUE))

# Calculate standard error
se_rlnd_orth<- aggregate(piNfract ~ Sex.Bias.Category + GeneType, 
                    data = SBGE_rlnd_bins_compare_orth, 
                    FUN = function(x) sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x))))

# Merge both summaries
summary_rlnd_data_comparative_orth <- merge(mean_rlnd_orth,se_rlnd_orth,
                                       by = c("Sex.Bias.Category", "GeneType"))

# Rename columns
colnames(summary_rlnd_data_comparative_orth) <- c("Sex.Bias.Category", "GeneType", "Mean_rlnd", "SE_rlnd")


rlnd_sample_sizes_orth <- as.data.frame(table(SBGE_rlnd_bins_compare_orth$Sex.Bias.Category,
                                         SBGE_rlnd_bins_compare_orth$GeneType))

# Rename columns
colnames(rlnd_sample_sizes_orth) <- c("Sex.Bias.Category", "GeneType", "SampleSize")

summary_rlnd_data_comparative_orth <- left_join(summary_rlnd_data_comparative_orth, rlnd_sample_sizes_orth, by = c("Sex.Bias.Category", "GeneType"))


# Split by GeneType
split_data_rlnd_orth <- split(summary_rlnd_data_comparative_orth, summary_rlnd_data_comparative_orth$GeneType)

# For each GeneType, compute "All" summary row
all_category_rlnd_orth <- do.call(rbind, lapply(split_data_rlnd_orth , function(df) {
  data.frame(
    GeneType = unique(df$GeneType),
    Mean_rlnd = mean(df$Mean_rlnd, na.rm = TRUE),
    SE_rlnd = sqrt(sum(df$SE_rlnd^2, na.rm = TRUE) / nrow(df)),
    SampleSize = sum(df$SampleSize, na.rm = TRUE),
    Sex.Bias.Category = "All"
  )
}))

# Reorder columns to match the expected format

all_category_rlnd_orth<- all_category_rlnd_orth[, c("Sex.Bias.Category", "GeneType", "Mean_rlnd", "SE_rlnd", "SampleSize")]

summary_rlnd_data_comparative_orth <- bind_rows(summary_rlnd_data_comparative_orth,all_category_rlnd_orth)

# ALL DISEASE GENES----
# Tajimas D Summary table ----

SBGE_tajDS_compare_AD<- candidate_all_disease_Singh%>%
  select(FlyBaseID,Sex.Bias.Category,exp_SBGE_ase,tajD.S,GeneType)%>% 
  distinct() # Remove disease variable to avoid repeats

SBGE_tajDS_compare_AD<- SBGE_tajDS_compare_AD %>%
  mutate(exp_SBGE_ase = as.numeric(exp_SBGE_ase))

# Creating Bins

SBGE_tajDS_bins_compare_AD<-SBGE_tajDS_compare_AD %>%
  mutate(Sex.Bias.Category = case_when(
    exp_SBGE_ase >= 2 ~ "Highly male-biased",
    exp_SBGE_ase >= 0.5 & exp_SBGE_ase < 2 ~ "Moderately male-biased",
    exp_SBGE_ase >= -0.5 & exp_SBGE_ase < 0.5 ~ "Unbiased",
    exp_SBGE_ase >= -2 & exp_SBGE_ase < -0.5 ~ "Moderately female-biased",
    exp_SBGE_ase < -2 ~ "Highly female-biased"
  ))

SBGE_tajDS_bins_compare_AD <- SBGE_tajDS_bins_compare_AD[!is.na(SBGE_tajDS_bins_compare_AD$tajD.S), ]

# Compute Mean and Standard Error (SE) 

# Calculate mean
mean_tajD_AD <- aggregate(tajD.S ~ Sex.Bias.Category + GeneType, 
                            data = SBGE_tajDS_bins_compare_AD, 
                            FUN = function(x) mean(x, na.rm = TRUE))

# Calculate standard error
se_tajD_AD <- aggregate(tajD.S ~ Sex.Bias.Category + GeneType, 
                          data = SBGE_tajDS_bins_compare_AD, 
                          FUN = function(x) sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x))))

# Combine into one data frame
summary_tajDS_data_comparative_AD <- merge(mean_tajD_AD, se_tajD_AD, 
                                             by = c("Sex.Bias.Category", "GeneType"))

# Rename columns
colnames(summary_tajDS_data_comparative_AD) <- c("Sex.Bias.Category", "GeneType", 
                                                   "Mean_tajD.S", "SE_tajD.S")

sample_sizes_tadj_AD <- as.data.frame(table(SBGE_tajDS_bins_compare_AD$Sex.Bias.Category,
                                              SBGE_tajDS_bins_compare_AD$GeneType))

# Rename columns
colnames(sample_sizes_tadj_AD) <- c("Sex.Bias.Category", "GeneType", "SampleSize")

summary_tajDS_data_comparative_AD <- left_join(summary_tajDS_data_comparative_AD, sample_sizes_tadj_AD, by = c("Sex.Bias.Category", "GeneType"))

# Then compute the "All" category stats

split_data_AD<- split(summary_tajDS_data_comparative_AD, summary_tajDS_data_comparative_AD$GeneType)

all_category_AD <- do.call(rbind, lapply(split_data_AD, function(df) {
  data.frame(
    GeneType = unique(df$GeneType),
    Mean_tajD.S = mean(df$Mean_tajD.S, na.rm = TRUE),
    SE_tajD.S = sqrt(sum(df$SE_tajD.S^2, na.rm = TRUE) / nrow(df)),
    SampleSize = sum(df$SampleSize, na.rm = TRUE),
    Sex.Bias.Category = "All"
  )
}))

summary_tajDS_data_comparative_AD <- bind_rows(summary_tajDS_data_comparative_AD, all_category_AD)

# rMF Summary table ----

# Change into numerical 

SBGE_rMF_compare_AD<- candidate_all_disease_Singh%>%
  select(FlyBaseID,Sex.Bias.Category,exp_SBGE_ase,rmf,GeneType)%>% 
  distinct() # Remove disease variable to avoid repeats

SBGE_rMF_compare_AD<- SBGE_rMF_compare_AD%>%
  mutate(exp_SBGE_ase = as.numeric(exp_SBGE_ase))

# Creating Bins

SBGE_rMF_bins_compare_AD<-SBGE_rMF_compare_AD%>%
  mutate(Sex.Bias.Category = case_when(
    exp_SBGE_ase >= 2 ~ "Highly male-biased",
    exp_SBGE_ase >= 0.5 & exp_SBGE_ase < 2 ~ "Moderately male-biased",
    exp_SBGE_ase >= -0.5 & exp_SBGE_ase < 0.5 ~ "Unbiased",
    exp_SBGE_ase >= -2 & exp_SBGE_ase < -0.5 ~ "Moderately female-biased",
    exp_SBGE_ase < -2 ~ "Highly female-biased"
  ))

SBGE_rMF_bins_compare_AD <- SBGE_rMF_bins_compare_AD[!is.na(SBGE_rMF_bins_compare_AD$rmf), ]

# Compute Mean and Standard Error (SE) 

# Compute mean
mean_rmf_AD <- aggregate(rmf ~ Sex.Bias.Category + GeneType, 
                           data = SBGE_rMF_bins_compare_AD, 
                           FUN = function(x) mean(x, na.rm = TRUE))

# Compute standard error
se_rmf_AD <- aggregate(rmf ~ Sex.Bias.Category + GeneType, 
                         data = SBGE_rMF_bins_compare_AD, 
                         FUN = function(x) sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x))))

# Merge the two summaries
summary_rMF_data_comparative_AD <- merge(mean_rmf_AD, se_rmf_AD, 
                                           by = c("Sex.Bias.Category", "GeneType"))

# Rename columns
colnames(summary_rMF_data_comparative_AD) <- c("Sex.Bias.Category", "GeneType", 
                                                 "Mean_rmf", "SE_rmf")

rMF_sample_sizes_AD <- as.data.frame(table(SBGE_rMF_bins_compare_AD$Sex.Bias.Category,
                                             SBGE_rMF_bins_compare_AD$GeneType))

# Rename columns 
colnames(rMF_sample_sizes_AD) <- c("Sex.Bias.Category", "GeneType", "SampleSize")
summary_rMF_data_comparative_AD <- left_join(summary_rMF_data_comparative_AD, rMF_sample_sizes_AD, by = c("Sex.Bias.Category", "GeneType"))

# Split by GeneType
split_data_rmf_AD <- split(summary_rMF_data_comparative_AD, summary_rMF_data_comparative_AD$GeneType)

# For each GeneType, compute "All" summary row
all_category_rmf_AD <- do.call(rbind, lapply(split_data_rmf_AD , function(df) {
  data.frame(
    GeneType = unique(df$GeneType),
    Mean_rmf = mean(df$Mean_rmf, na.rm = TRUE),
    SE_rmf = sqrt(sum(df$SE_rmf^2, na.rm = TRUE) / nrow(df)),
    SampleSize = sum(df$SampleSize, na.rm = TRUE),
    Sex.Bias.Category = "All"
  )
}))

summary_rMF_data_comparative_AD <- bind_rows(summary_rMF_data_comparative_AD,all_category_rmf_AD)


# DoS Summary table ----
# Change into numerical 

SBGE_DoS_compare_AD<-candidate_all_disease_Singh %>%
  select(FlyBaseID,Sex.Bias.Category,exp_SBGE_ase,DoS,GeneType)%>% 
  distinct() # Remove disease variable to avoid repeats


SBGE_DoS_compare_AD<- SBGE_DoS_compare_AD%>%
  mutate(exp_SBGE_ase = as.numeric(exp_SBGE_ase))

SBGE_DoS_compare_AD <- SBGE_DoS_compare_AD[!is.na(SBGE_DoS_compare_AD$DoS), ]

# Creating Bins

SBGE_DoS_bins_compare_AD<-SBGE_DoS_compare_AD%>%
  mutate(Sex.Bias.Category = case_when(
    exp_SBGE_ase >= 2 ~ "Highly male-biased",
    exp_SBGE_ase >= 0.5 & exp_SBGE_ase < 2 ~ "Moderately male-biased",
    exp_SBGE_ase >= -0.5 & exp_SBGE_ase < 0.5 ~ "Unbiased",
    exp_SBGE_ase >= -2 & exp_SBGE_ase < -0.5 ~ "Moderately female-biased",
    exp_SBGE_ase < -2 ~ "Highly female-biased"
  ))

# Compute Mean and Standard Error (SE) 

# Calculate mean
mean_DoS_AD <- aggregate(DoS ~ Sex.Bias.Category + GeneType, 
                           data = SBGE_DoS_bins_compare_AD, 
                           FUN = function(x) mean(x, na.rm = TRUE))

# Calculate standard error
se_DoS_AD <- aggregate(DoS ~ Sex.Bias.Category + GeneType, 
                         data = SBGE_DoS_bins_compare_AD, 
                         FUN = function(x) sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x))))

# Merge both summaries
summary_DoS_data_comparative_AD <- merge(mean_DoS_AD, se_DoS_AD,
                                           by = c("Sex.Bias.Category", "GeneType"))

# Rename columns
colnames(summary_DoS_data_comparative_AD) <- c("Sex.Bias.Category", "GeneType", "Mean_DoS", "SE_DoS")


DoS_sample_sizes_AD <- as.data.frame(table(SBGE_DoS_bins_compare_AD$Sex.Bias.Category,
                                             SBGE_DoS_bins_compare_AD$GeneType))

# Rename columns
colnames(DoS_sample_sizes_AD) <- c("Sex.Bias.Category", "GeneType", "SampleSize")

summary_DoS_data_comparative_AD <- left_join(summary_DoS_data_comparative_AD, DoS_sample_sizes_AD, by = c("Sex.Bias.Category", "GeneType"))


# Split by GeneType
split_data_DoS_AD <- split(summary_DoS_data_comparative_AD, summary_DoS_data_comparative_AD$GeneType)

# For each GeneType, compute "All" summary row
all_category_DoS_AD <- do.call(rbind, lapply(split_data_DoS_AD , function(df) {
  data.frame(
    GeneType = unique(df$GeneType),
    Mean_DoS = mean(df$Mean_DoS, na.rm = TRUE),
    SE_DoS = sqrt(sum(df$SE_DoS^2, na.rm = TRUE) / nrow(df)),
    SampleSize = sum(df$SampleSize, na.rm = TRUE),
    Sex.Bias.Category = "All"
  )
}))

# Reorder columns to match the expected format
all_category_DoS_AD<- all_category_DoS_AD [, c("Sex.Bias.Category", "GeneType", "Mean_DoS", "SE_DoS", "SampleSize")]


summary_DoS_data_comparative_AD <- bind_rows(summary_DoS_data_comparative_AD,all_category_DoS_AD)


# Tissue Specificity Summary table  ----

# Change into numerical 

SBGE_tissuespec_compare_AD<- candidate_all_disease_Fra%>%
  select(FlyBaseID,Sex.Bias.Category,exp_SBGE_ase,tissueSpecificityAdults,GeneType)%>% 
  distinct() # Remove disease variable to avoid repeats

SBGE_tissuespec_compare_AD<- SBGE_tissuespec_compare_AD%>%
  mutate(exp_SBGE_ase = as.numeric(exp_SBGE_ase))

SBGE_tissuespec_compare_AD <- SBGE_tissuespec_compare_AD[!is.na(SBGE_tissuespec_compare_AD$tissueSpecificityAdults), ]


# Creating Bins

SBGE_tissuespec_bins_compare_AD<-SBGE_tissuespec_compare_AD%>%
  mutate(Sex.Bias.Category = case_when(
    exp_SBGE_ase >= 2 ~ "Highly male-biased",
    exp_SBGE_ase >= 0.5 & exp_SBGE_ase < 2 ~ "Moderately male-biased",
    exp_SBGE_ase >= -0.5 & exp_SBGE_ase < 0.5 ~ "Unbiased",
    exp_SBGE_ase >= -2 & exp_SBGE_ase < -0.5 ~ "Moderately female-biased",
    exp_SBGE_ase < -2 ~ "Highly female-biased"
  ))

SBGE_tissuespec_bins_compare_AD <- SBGE_tissuespec_bins_compare_AD[!is.na(SBGE_tissuespec_bins_compare_AD$tissueSpecificityAdults), ]

# Compute Mean and Standard Error (SE) 

# Calculate mean
mean_tissueSpecificityAdults_AD<- aggregate(tissueSpecificityAdults~ Sex.Bias.Category + GeneType, 
                                              data = SBGE_tissuespec_bins_compare_AD, 
                                              FUN = function(x) mean(x, na.rm = TRUE))

# Calculate standard error
se_tissueSpecificityAdults_AD<- aggregate(tissueSpecificityAdults ~ Sex.Bias.Category + GeneType, 
                                            data = SBGE_tissuespec_bins_compare_AD, 
                                            FUN = function(x) sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x))))

# Merge both summaries
summary_tissuespec_data_comparative_AD <- merge(mean_tissueSpecificityAdults_AD,se_tissueSpecificityAdults_AD,
                                                  by = c("Sex.Bias.Category", "GeneType"))

# Rename columns
colnames(summary_tissuespec_data_comparative_AD) <- c("Sex.Bias.Category", "GeneType", "Mean_tissueSpecificityAdults", "SE_tissueSpecificityAdults")


tissuespec_sample_sizes_AD <- as.data.frame(table(SBGE_tissuespec_bins_compare_AD$Sex.Bias.Category,
                                                    SBGE_tissuespec_bins_compare_AD$GeneType))

# Rename columns
colnames(tissuespec_sample_sizes_AD) <- c("Sex.Bias.Category", "GeneType", "SampleSize")

summary_tissuespec_data_comparative_AD <- left_join(summary_tissuespec_data_comparative_AD, tissuespec_sample_sizes_AD, by = c("Sex.Bias.Category", "GeneType"))


# Split by GeneType
split_data_tissuespec_AD <- split(summary_tissuespec_data_comparative_AD, summary_tissuespec_data_comparative_AD$GeneType)

# For each GeneType, compute "All" summary row
all_category_tissuespec_AD <- do.call(rbind, lapply(split_data_tissuespec_AD , function(df) {
  data.frame(
    GeneType = unique(df$GeneType),
    Mean_tissueSpecificityAdults = mean(df$Mean_tissueSpecificityAdults, na.rm = TRUE),
    SE_tissueSpecificityAdults  = sqrt(sum(df$SE_tissueSpecificityAdults^2, na.rm = TRUE) / nrow(df)),
    SampleSize = sum(df$SampleSize, na.rm = TRUE),
    Sex.Bias.Category = "All"
  )
}))

# Reorder columns to match the expected format
all_category_tissuespec_AD <- all_category_tissuespec_AD[, c("Sex.Bias.Category", "GeneType", "Mean_tissueSpecificityAdults", "SE_tissueSpecificityAdults", "SampleSize")]


summary_tissuespec_data_comparative_AD <- bind_rows(summary_tissuespec_data_comparative_AD,all_category_tissuespec_AD)

# Relative Level of Nonsynonymous Diversity Summary table  ----

# Change into numerical 

SBGE_rlnd_compare_AD<- candidate_all_disease_Singh%>%
  select(FlyBaseID,Sex.Bias.Category,exp_SBGE_ase,piNfract,GeneType)%>% 
  distinct() # Remove disease variable to avoid repeats

SBGE_rlnd_compare_AD<- SBGE_rlnd_compare_AD%>%
  mutate(exp_SBGE_ase = as.numeric(exp_SBGE_ase))

# Creating Bins

SBGE_rlnd_bins_compare_AD<-SBGE_rlnd_compare_AD%>%
  mutate(Sex.Bias.Category = case_when(
    exp_SBGE_ase >= 2 ~ "Highly male-biased",
    exp_SBGE_ase >= 0.5 & exp_SBGE_ase < 2 ~ "Moderately male-biased",
    exp_SBGE_ase >= -0.5 & exp_SBGE_ase < 0.5 ~ "Unbiased",
    exp_SBGE_ase >= -2 & exp_SBGE_ase < -0.5 ~ "Moderately female-biased",
    exp_SBGE_ase < -2 ~ "Highly female-biased"
  ))

SBGE_rlnd_bins_compare_AD <- SBGE_rlnd_bins_compare_AD[!is.na(SBGE_rlnd_bins_compare_AD$piNfract), ]

# Compute Mean and Standard Error (SE) 

# Calculate mean
mean_rlnd_AD<- aggregate(piNfract~ Sex.Bias.Category + GeneType, 
                           data = SBGE_rlnd_bins_compare_AD, 
                           FUN = function(x) mean(x, na.rm = TRUE))

# Calculate standard error
se_rlnd_AD<- aggregate(piNfract ~ Sex.Bias.Category + GeneType, 
                         data = SBGE_rlnd_bins_compare_AD, 
                         FUN = function(x) sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x))))

# Merge both summaries
summary_rlnd_data_comparative_AD <- merge(mean_rlnd_AD,se_rlnd_AD,
                                            by = c("Sex.Bias.Category", "GeneType"))

# Rename columns
colnames(summary_rlnd_data_comparative_AD) <- c("Sex.Bias.Category", "GeneType", "Mean_rlnd", "SE_rlnd")


rlnd_sample_sizes_AD <- as.data.frame(table(SBGE_rlnd_bins_compare_AD$Sex.Bias.Category,
                                              SBGE_rlnd_bins_compare_AD$GeneType))

# Rename columns
colnames(rlnd_sample_sizes_AD) <- c("Sex.Bias.Category", "GeneType", "SampleSize")

summary_rlnd_data_comparative_AD <- left_join(summary_rlnd_data_comparative_AD, rlnd_sample_sizes_AD, by = c("Sex.Bias.Category", "GeneType"))


# Split by GeneType
split_data_rlnd_AD <- split(summary_rlnd_data_comparative_AD, summary_rlnd_data_comparative_AD$GeneType)

# For each GeneType, compute "All" summary row
all_category_rlnd_AD <- do.call(rbind, lapply(split_data_rlnd_AD , function(df) {
  data.frame(
    GeneType = unique(df$GeneType),
    Mean_rlnd = mean(df$Mean_rlnd, na.rm = TRUE),
    SE_rlnd = sqrt(sum(df$SE_rlnd^2, na.rm = TRUE) / nrow(df)),
    SampleSize = sum(df$SampleSize, na.rm = TRUE),
    Sex.Bias.Category = "All"
  )
}))

# Reorder columns to match the expected format

all_category_rlnd_AD<- all_category_rlnd_AD[, c("Sex.Bias.Category", "GeneType", "Mean_rlnd", "SE_rlnd", "SampleSize")]

summary_rlnd_data_comparative_AD <- bind_rows(summary_rlnd_data_comparative_AD,all_category_rlnd_AD)

#___________________________-----
# STATISTICAL ANALYSES----

# PERMUTATION TESTS ----

# Adding logical variable enables function to work

summary_tajDS_data_comparative_orth$GeneType_logical <- summary_tajDS_data_comparative_orth$GeneType == "Ortholog"
str(summary_tajDS_data_comparative_orth)

SBGE_tajDS_bins_compare_orth$GeneType_logical <- SBGE_tajDS_bins_compare_orth$GeneType == "Ortholog"
str(SBGE_tajDS_bins_compare_orth)

summary_rMF_data_comparative_orth$GeneType_logical <- summary_rMF_data_comparative_orth$GeneType == "Ortholog"
str(summary_rMF_data_comparative_orth)

SBGE_rMF_bins_compare_orth$GeneType_logical <- SBGE_rMF_bins_compare_orth$GeneType == "Ortholog"
str(SBGE_rMF_bins_compare_orth)

summary_DoS_data_comparative_orth$GeneType_logical <- summary_DoS_data_comparative_orth$GeneType == "Ortholog"
str(summary_DoS_data_comparative_orth)

SBGE_DoS_bins_compare_orth$GeneType_logical <- SBGE_DoS_bins_compare_orth$GeneType == "Ortholog"
str(SBGE_DoS_bins_compare_orth)

summary_tissuespec_data_comparative_orth$GeneType_logical <- summary_tissuespec_data_comparative_orth$GeneType == "Ortholog"
str(summary_tissuespec_data_comparative_orth)

SBGE_tissuespec_bins_compare_orth$GeneType_logical <- SBGE_tissuespec_bins_compare_orth$GeneType == "Ortholog"
str(SBGE_tissuespec_bins_compare_orth)

summary_rlnd_data_comparative_orth$GeneType_logical <- summary_rlnd_data_comparative_orth$GeneType == "Ortholog"
str(summary_rlnd_data_comparative_orth)

SBGE_rlnd_bins_compare_orth$GeneType_logical <- SBGE_rlnd_bins_compare_orth$GeneType == "Ortholog"
str(SBGE_rlnd_bins_compare_orth)

summary_tajDS_data_comparative_AD$GeneType_logical <- summary_tajDS_data_comparative_AD$GeneType == "Disease"
str(summary_tajDS_data_comparative_AD)

SBGE_tajDS_bins_compare_AD$GeneType_logical <- SBGE_tajDS_bins_compare_AD$GeneType == "Disease"
str(SBGE_tajDS_bins_compare_AD)

summary_rMF_data_comparative_AD$GeneType_logical <- summary_rMF_data_comparative_AD$GeneType == "Disease"
str(summary_rMF_data_comparative_AD)

SBGE_rMF_bins_compare_AD$GeneType_logical <- SBGE_rMF_bins_compare_AD$GeneType == "Disease"
str(SBGE_rMF_bins_compare_AD)

summary_DoS_data_comparative_AD$GeneType_logical <- summary_DoS_data_comparative_AD$GeneType == "Disease"
str(summary_DoS_data_comparative_AD)

SBGE_DoS_bins_compare_AD$GeneType_logical <- SBGE_DoS_bins_compare_AD$GeneType == "Disease"
str(SBGE_DoS_bins_compare_AD)

summary_tissuespec_data_comparative_AD$GeneType_logical <- summary_tissuespec_data_comparative_AD$GeneType == "Disease"
str(summary_tissuespec_data_comparative_AD)

SBGE_tissuespec_bins_compare_AD$GeneType_logical <- SBGE_tissuespec_bins_compare_AD$GeneType == "Disease"
str(SBGE_tissuespec_bins_compare_AD)

summary_rlnd_data_comparative_AD$GeneType_logical <- summary_rlnd_data_comparative_AD$GeneType == "Disease"
str(summary_rlnd_data_comparative_AD)

SBGE_rlnd_bins_compare_AD$GeneType_logical <- SBGE_rlnd_bins_compare_AD$GeneType == "Disease"
str(SBGE_rlnd_bins_compare_AD)

# TwoPerm
# Two sample permutation test for difference in means
# args: perm_dat = the dataset containing a numeric column to permute
#       x_col = string denoting the column containing the data to permute
#       groupBy = logical column in the dataset which separates group 1 and group 2
#       n_perm = number of permutations. default = 10,000
#       alternative = alt. for hypothesis testing
TwoPerm <- function(perm_dat, x_col, 
                    groupBy, n_perm = 10000,
                    alternative = c("two.sided","less","greater")){
  # separate data for groups 1 and 2
  trt1 <- perm_dat[perm_dat[[groupBy]],][[x_col]]
  trt2 <- perm_dat[!perm_dat[[groupBy]],][[x_col]]
  # number of observations
  n1 <- length(trt1)
  n2 <- length(trt2)
  
  # observed difference in mean
  obs.diff <- mean(trt1) - mean(trt2)
  # shuffle data from trt1 and trt2
  rand.samples <- sample(c(trt1, trt2), n1 + n2, replace = F)
  
  # initialize dataframe for permuted groups
  permuted.trt1 <- c()
  permuted.trt2 <- c()
  diff.permuted <- c()
  for (i in 1:n_perm){
    # shuffle and build permuted test stat for each group
    permuted.trt1 <- sample(rand.samples, n1, replace = F)
    permuted.trt2 <- sample(rand.samples, n2, replace = F)
    
    # calculated difference in means for each permutation
    diff.permuted <- append(diff.permuted, mean(permuted.trt1) - mean(permuted.trt2))
  }
  
  # assign p-value
  if(alternative[1]=="less"){
    pval = (sum(diff.permuted <= obs.diff) + 1) / (n_perm + 1)
  } else if(alternative[1]=="greater"){
    pval = (sum(diff.permuted >= obs.diff) + 1) / (n_perm + 1)
  } else{
    pval = (sum(abs(diff.permuted) > abs(obs.diff)) + 1) / (n_perm + 1)
  }
  
  
  summary.test <- c(obs.diff, n1, n2, pval, (pval < 0.05 & n1 > 30 & n2 > 30))
  return(summary.test)
}


PairedTwoPerm <- function(perm_dat, x1, x2, n_perm = 10000,
                          alternative = c("two.sided","less","greater")){
  
  summary.test <- NULL # return object
  perm_dat <- perm_dat[!is.na(perm_dat[[x1]]) & 
                         !is.na(perm_dat[[x2]]),] # only keep paired observations
  
  n1 <- length(na.omit(perm_dat[[x1]]))
  n2 <- length(na.omit(perm_dat[[x2]]))
  
  if(n1 == n2) {
    # find difference for each paired observation
    perm_dat <- perm_dat %>%
      dplyr::mutate(Diff = .[[x1]] - .[[x2]])
    
    # observed mean difference
    obs.diff <- mean(perm_dat$Diff)
    
    diff.permuted <- rep(12345, n_perm)
    
    for (i in 1:n_perm){
      signs <- sample(c(-1,1), size = n1, replace = TRUE)
      perm_diffs = signs*perm_dat$Diff
      diff.permuted[i] = mean(perm_diffs)
    }
    
    # assign p-value
    if(alternative[1]=="less"){
      pval = (sum(diff.permuted <= obs.diff) + 1) / (n_perm + 1)
    } else if(alternative[1]=="greater"){
      pval = (sum(diff.permuted >= obs.diff) + 1) / (n_perm + 1)
    } else{
      pval = (sum(abs(diff.permuted) > abs(obs.diff)) + 1) / (n_perm + 1)
    }
    
    
    summary.test <- data.frame(obs.stat = obs.diff, 
                               N = n1, 
                               pval = pval, 
                               Sig = (pval < 0.05 & n1 > 30 & n2 > 30))
  }
  
  return(summary.test)
}


# Function to permute based on SBGE category
# args: perm_dat = the dataset containing a numeric column to permute
#       x_col = string denoting denoting the column containing the data to permute
#       groupBy = logical column which separates group 1 and group 2
#       SBGE_cat = column specifying SBGE bins
TwoPerm_SBGE <- function(perm_dat, x_col, groupBy, SBGE_cat){
  dat <- data.frame() # initialize data.frame object to store results
  # make sure the column of categories is set to factor
  perm_dat[[SBGE_cat]] <- factor(perm_dat[[SBGE_cat]])
  
  # loop through each sex-bias category
  for(i in levels(perm_dat[[SBGE_cat]])){
    # for each SBGE category, call the permutation function
    comp <- c(i, TwoPerm(perm_dat[perm_dat[[SBGE_cat]] == i,], 
                         x_col = x_col, 
                         groupBy = groupBy, 
                         alternative = "two.sided"))
    
    # concatenate results in the data.frame to return
    dat <- rbind(dat, comp)
  }
  colnames(dat) <- c(SBGE_cat, "obs.diff", "n_TRUE", "n_FALSE", "pval", "Sig")
  dat$obs.diff <- as.numeric(dat$obs.diff)
  dat$n_TRUE <- as.numeric(dat$n_TRUE)
  dat$n_FALSE <- as.numeric(dat$n_FALSE)
  return(dat)
}

#ORTHOLOGS ----
#TajD perm----

# Now run TwoPerm
result_tajD_orth <- TwoPerm(
  perm_dat = SBGE_tajDS_bins_compare_orth,
  x_col = "tajD.S",
  groupBy = "GeneType_logical"
)

# Now run TwoPerm_SBGE() — it will include "All"
result_tajD_SBGE_cat_orth <- TwoPerm_SBGE(
  perm_dat = SBGE_tajDS_bins_compare_orth,
  x_col = "tajD.S",
  groupBy = "GeneType_logical",
  SBGE_cat = "Sex.Bias.Category"
)

result_tajD_all_orth <- data.frame(
  Sex.Bias.Category = "All",
  obs.diff = as.numeric(result_tajD_orth[1]),
  n_TRUE = as.numeric(result_tajD_orth[2]),
  n_FALSE = as.numeric(result_tajD_orth[3]),
  pval = as.character(result_tajD_orth[4]),
  Sig = as.character(result_tajD_orth[5])
)

tajD_perm_results_orth<- bind_rows(result_tajD_all_orth, result_tajD_SBGE_cat_orth )

#rMF perm----
# Now run TwoPerm
result_rmf_orth<- TwoPerm(
  perm_dat = SBGE_rMF_bins_compare_orth,
  x_col = "rmf",
  groupBy = "GeneType_logical"
)

print(result_rmf_orth)

# Now run TwoPerm_SBGE
result_rmf_SBGE_cat_orth <- TwoPerm_SBGE(
  perm_dat = SBGE_rMF_bins_compare_orth,
  x_col = "rmf",
  groupBy = "GeneType_logical",
  SBGE_cat = "Sex.Bias.Category"
)

result_rmf_all_orth <- data.frame(
  Sex.Bias.Category = "All",
  obs.diff = as.numeric(result_rmf_orth[1]),
  n_TRUE = as.numeric(result_rmf_orth[2]),
  n_FALSE = as.numeric(result_rmf_orth[3]),
  pval = as.character(result_rmf_orth[4]),
  Sig = as.character(result_rmf_orth[5])
)

rmf_perm_results_orth<- bind_rows(result_rmf_all_orth, result_rmf_SBGE_cat_orth)



#DoS perm----

# Now run TwoPerm
result_DoS_orth<- TwoPerm(
  perm_dat = SBGE_DoS_bins_compare_orth,
  x_col = "DoS",
  groupBy = "GeneType_logical"
)

# Now run TwoPerm_SBGE() — it will include "All"
result_DoS_SBGE_cat_orth <- TwoPerm_SBGE(
  perm_dat = SBGE_DoS_bins_compare_orth,
  x_col = "DoS",
  groupBy = "GeneType_logical",
  SBGE_cat = "Sex.Bias.Category"
)

result_DoS_all_orth <- data.frame(
  Sex.Bias.Category = "All",
  obs.diff = as.numeric(result_DoS_orth[1]),
  n_TRUE = as.numeric(result_DoS_orth[2]),
  n_FALSE = as.numeric(result_DoS_orth[3]),
  pval = as.character(result_DoS_orth[4]),
  Sig = as.character(result_DoS_orth[5])
)

DoS_perm_results_orth<- bind_rows(result_DoS_all_orth, result_DoS_SBGE_cat_orth )
print(result_DoS_orth)

# Tissue Specificity perm----
# Now run TwoPerm
result_TissueSpecificity_orth<- TwoPerm(
  perm_dat = SBGE_tissuespec_bins_compare_orth,
  x_col = "tissueSpecificityAdults",
  groupBy = "GeneType_logical"
)

# Now run TwoPerm_SBGE() — it will include "All"
result_TissueSpecificity_SBGE_cat_orth <- TwoPerm_SBGE(
  perm_dat = SBGE_tissuespec_bins_compare_orth,
  x_col = "tissueSpecificityAdults",
  groupBy = "GeneType_logical",
  SBGE_cat = "Sex.Bias.Category"
)

result_TissueSpecificity_all_orth <- data.frame(
  Sex.Bias.Category = "All",
  obs.diff = as.numeric(result_TissueSpecificity_orth[1]),
  n_TRUE = as.numeric(result_TissueSpecificity_orth[2]),
  n_FALSE = as.numeric(result_TissueSpecificity_orth[3]),
  pval = as.character(result_TissueSpecificity_orth[4]),
  Sig = as.character(result_TissueSpecificity_orth[5])
)

TissueSpecificity_perm_results_orth<- bind_rows(result_TissueSpecificity_all_orth, result_TissueSpecificity_SBGE_cat_orth )

# Relative level of Nonsynonymous Diversity----

# Now run TwoPerm
result_rlnd_orth<- TwoPerm(
  perm_dat = SBGE_rlnd_bins_compare_orth,
  x_col = "piNfract",
  groupBy = "GeneType_logical"
)

# Now run TwoPerm_SBGE() — it will include "All"
result_rlnd_SBGE_cat_orth <- TwoPerm_SBGE(
  perm_dat = SBGE_rlnd_bins_compare_orth,
  x_col = "piNfract",
  groupBy = "GeneType_logical",
  SBGE_cat = "Sex.Bias.Category"
)

result_rlnd_all_orth <- data.frame(
  Sex.Bias.Category = "All",
  obs.diff = as.numeric(result_rlnd_orth[1]),
  n_TRUE = as.numeric(result_rlnd_orth[2]),
  n_FALSE = as.numeric(result_rlnd_orth[3]),
  pval = as.character(result_rlnd_orth[4]),
  Sig = as.character(result_rlnd_orth[5])
)

rlnd_perm_results_orth<- bind_rows(result_rlnd_all_orth , result_rlnd_SBGE_cat_orth )

# ALL DISEASE ----
# TajD perm----

# Now run TwoPerm
result_tajD_AD <- TwoPerm(
  perm_dat = SBGE_tajDS_bins_compare_AD,
  x_col = "tajD.S",
  groupBy = "GeneType_logical"
)

# Now run TwoPerm_SBGE() — it will include "All"
result_tajD_SBGE_cat_AD <- TwoPerm_SBGE(
  perm_dat = SBGE_tajDS_bins_compare_AD,
  x_col = "tajD.S",
  groupBy = "GeneType_logical",
  SBGE_cat = "Sex.Bias.Category"
)

result_tajD_all_AD <- data.frame(
  Sex.Bias.Category = "All",
  obs.diff = as.numeric(result_tajD_AD[1]),
  n_TRUE = as.numeric(result_tajD_AD[2]),
  n_FALSE = as.numeric(result_tajD_AD[3]),
  pval = as.character(result_tajD_AD[4]),
  Sig = as.character(result_tajD_AD[5])
)

tajD_perm_results_AD<- bind_rows(result_tajD_all_AD, result_tajD_SBGE_cat_AD )

#rMF perm----
# Now run TwoPerm
result_rmf_AD<- TwoPerm(
  perm_dat = SBGE_rMF_bins_compare_AD,
  x_col = "rmf",
  groupBy = "GeneType_logical"
)

print(result_rmf_AD)

# Now run TwoPerm_SBGE
result_rmf_SBGE_cat_AD <- TwoPerm_SBGE(
  perm_dat = SBGE_rMF_bins_compare_AD,
  x_col = "rmf",
  groupBy = "GeneType_logical",
  SBGE_cat = "Sex.Bias.Category"
)

result_rmf_all_AD <- data.frame(
  Sex.Bias.Category = "All",
  obs.diff = as.numeric(result_rmf_AD[1]),
  n_TRUE = as.numeric(result_rmf_AD[2]),
  n_FALSE = as.numeric(result_rmf_AD[3]),
  pval = as.character(result_rmf_AD[4]),
  Sig = as.character(result_rmf_AD[5])
)

rmf_perm_results_AD<- bind_rows(result_rmf_all_AD, result_rmf_SBGE_cat_AD)



#DoS perm----

# Now run TwoPerm
result_DoS_AD<- TwoPerm(
  perm_dat = SBGE_DoS_bins_compare_AD,
  x_col = "DoS",
  groupBy = "GeneType_logical"
)

# Now run TwoPerm_SBGE() — it will include "All"
result_DoS_SBGE_cat_AD <- TwoPerm_SBGE(
  perm_dat = SBGE_DoS_bins_compare_AD,
  x_col = "DoS",
  groupBy = "GeneType_logical",
  SBGE_cat = "Sex.Bias.Category"
)

result_DoS_all_AD <- data.frame(
  Sex.Bias.Category = "All",
  obs.diff = as.numeric(result_DoS_AD[1]),
  n_TRUE = as.numeric(result_DoS_AD[2]),
  n_FALSE = as.numeric(result_DoS_AD[3]),
  pval = as.character(result_DoS_AD[4]),
  Sig = as.character(result_DoS_AD[5])
)

DoS_perm_results_AD<- bind_rows(result_DoS_all_AD, result_DoS_SBGE_cat_AD )
print(result_DoS_AD)

# Tissue Specificity perm----
# Now run TwoPerm
result_TissueSpecificity_AD<- TwoPerm(
  perm_dat = SBGE_tissuespec_bins_compare_AD,
  x_col = "tissueSpecificityAdults",
  groupBy = "GeneType_logical"
)

# Now run TwoPerm_SBGE() — it will include "All"
result_TissueSpecificity_SBGE_cat_AD <- TwoPerm_SBGE(
  perm_dat = SBGE_tissuespec_bins_compare_AD,
  x_col = "tissueSpecificityAdults",
  groupBy = "GeneType_logical",
  SBGE_cat = "Sex.Bias.Category"
)

result_TissueSpecificity_all_AD <- data.frame(
  Sex.Bias.Category = "All",
  obs.diff = as.numeric(result_TissueSpecificity_AD[1]),
  n_TRUE = as.numeric(result_TissueSpecificity_AD[2]),
  n_FALSE = as.numeric(result_TissueSpecificity_AD[3]),
  pval = as.character(result_TissueSpecificity_AD[4]),
  Sig = as.character(result_TissueSpecificity_AD[5])
)

TissueSpecificity_perm_results_AD<- bind_rows(result_TissueSpecificity_all_AD, result_TissueSpecificity_SBGE_cat_AD )

# Relative level of Nonsynonymous Diversity----

# Now run TwoPerm
result_rlnd_AD<- TwoPerm(
  perm_dat = SBGE_rlnd_bins_compare_AD,
  x_col = "piNfract",
  groupBy = "GeneType_logical"
)

# Now run TwoPerm_SBGE() — it will include "All"
result_rlnd_SBGE_cat_AD <- TwoPerm_SBGE(
  perm_dat = SBGE_rlnd_bins_compare_AD,
  x_col = "piNfract",
  groupBy = "GeneType_logical",
  SBGE_cat = "Sex.Bias.Category"
)

result_rlnd_all_AD <- data.frame(
  Sex.Bias.Category = "All",
  obs.diff = as.numeric(result_rlnd_AD[1]),
  n_TRUE = as.numeric(result_rlnd_AD[2]),
  n_FALSE = as.numeric(result_rlnd_AD[3]),
  pval = as.character(result_rlnd_AD[4]),
  Sig = as.character(result_rlnd_AD[5])
)

rlnd_perm_results_AD<- bind_rows(result_rlnd_all_AD , result_rlnd_SBGE_cat_AD )

# BOOTSTRAPPING FUNCTIONS ----

# Bootstrapped 95% confidence intervals for two groups
# args: boot_dat = dataset to bootstrap
#       x_col = string denoting the column containing the data to bootstrap
#       groupBy = logical column which separates group 1 and group 2
#       boot_n = number of bootstrap replicates. Default is 1,000
#       myfun = test stat function. Default is mean

TwoBoot <- function(boot_dat, x_col, groupBy, 
                    boot_n = 1000, 
                    myfun = mean){
  boot_tabs <- as_tibble(data_frame(bs = 1:boot_n) %>% # make boot_n bootstrap replicates
                           dplyr::group_by(bs) %>% # for each bootstrap,
                           # sample randomly x boot_n times from each group
                           dplyr::mutate(data = list(boot_dat %>% 
                                                       dplyr::group_by(.[[groupBy]]) %>% # separate group 1 and group 2
                                                       dplyr::sample_frac(size = 1, replace = T)))) %>% 
    unnest(c(bs, data)) %>% # create separate rows for each bootstrap replicate in the list
    dplyr::group_by(bs, .[[groupBy]]) %>% # group the data by bootstrap replicate and sig/non-sig
    # for each bootstrap replicate, calculate the test stat
    dplyr::do(tidy(myfun(.[[x_col]]))) %>%
    dplyr::rename({{groupBy}} := 2)
  
  # summarise bootstrap replicates
  boot_SE <- boot_tabs %>%
    dplyr::ungroup() %>%
    dplyr::group_by(.[[groupBy]]) %>%
    dplyr::summarise(q05 = quantile(x, 0.025, na.rm = TRUE),
                     q50 = quantile(x, 0.5, na.rm = TRUE),
                     q95 = quantile(x, 0.975, na.rm = TRUE)) %>%
    dplyr::rename({{groupBy}} := 1)
  return(boot_SE)
}


TwoBoot_SBGE <- function(boot_dat, x_col, groupBy, SBGE_cat){
  dat <- data.frame() # initialize data.frame object to store results
  # make sure the column of categories is set to factor
  boot_dat[[SBGE_cat]] <- factor(boot_dat[[SBGE_cat]])
  comp <- c()
  # loop through each sex-bias category
  for(i in levels(boot_dat[[SBGE_cat]])){
    # for each SBGE category, call the permutation function
    comp <- TwoBoot(boot_dat[boot_dat[[SBGE_cat]] == i,], 
                    x_col = x_col, 
                    groupBy = groupBy)
    comp[[SBGE_cat]] <- i
    # concatenate results in the data.frame to return
    dat <- rbind(dat, comp)
  }
  colnames(dat) <- c("Sig", "q05", "q50", "q95", SBGE_cat)
  return(dat)
}



#ORTHOLOGS----
#tajD boot----
tajD_boot_orth<-TwoBoot(SBGE_tajDS_bins_compare_orth, x_col = "tajD.S", groupBy = "GeneType", boot_n = 10)
tajD_boot_SBGE_orth <-TwoBoot_SBGE(SBGE_tajDS_bins_compare_orth, x_col = "tajD.S", groupBy = "GeneType",SBGE_cat = "Sex.Bias.Category")

# Add "All" category to first dataframe
tajD_boot_orth$Sex.Bias.Category <- "All"
names(tajD_boot_SBGE_orth)[names(tajD_boot_SBGE_orth) == "Sig"] <- "GeneType"

tajD_boot_results_orth<- bind_rows(tajD_boot_orth, tajD_boot_SBGE_orth )

#rMF boot----
rmf_boot_orth<-TwoBoot(SBGE_rMF_bins_compare_orth, x_col = "rmf", groupBy = "GeneType", boot_n = 10)
rmf_boot_SBGE_orth <-TwoBoot_SBGE(SBGE_rMF_bins_compare_orth, x_col = "rmf", groupBy = "GeneType",SBGE_cat = "Sex.Bias.Category")

# Add "All" category to first dataframe
rmf_boot_orth$Sex.Bias.Category <- "All"

names(rmf_boot_SBGE_orth)[names(rmf_boot_SBGE_orth) == "Sig"] <- "GeneType"


rmf_boot_results_orth<- bind_rows(rmf_boot_orth, rmf_boot_SBGE_orth )


#DoS boot----
DoS_boot_orth<-TwoBoot(SBGE_DoS_bins_compare_orth, x_col = "DoS", groupBy = "GeneType", boot_n = 10)
DoS_boot_SBGE_orth <-TwoBoot_SBGE(SBGE_DoS_bins_compare_orth, x_col = "DoS", groupBy = "GeneType",SBGE_cat = "Sex.Bias.Category")

# Add "All" category to first dataframe
DoS_boot_orth$Sex.Bias.Category <- "All"
names(DoS_boot_SBGE_orth)[names(DoS_boot_SBGE_orth) == "Sig"] <- "GeneType"

DoS_boot_results_orth<- bind_rows(DoS_boot_orth, DoS_boot_SBGE_orth )



#Tissue Specificity boot----
TissueSpecificity_boot_orth<-TwoBoot(SBGE_tissuespec_bins_compare_orth, x_col = "tissueSpecificityAdults", groupBy = "GeneType", boot_n = 10)
TissueSpecificity_boot_SBGE_orth <-TwoBoot_SBGE(SBGE_tissuespec_bins_compare_orth, x_col = "tissueSpecificityAdults", groupBy = "GeneType",SBGE_cat = "Sex.Bias.Category")

# Add "All" category to first dataframe
TissueSpecificity_boot_orth$Sex.Bias.Category <- "All"
names(TissueSpecificity_boot_SBGE_orth)[names(TissueSpecificity_boot_SBGE_orth) == "Sig"] <- "GeneType"

TissueSpecificity_boot_results_orth<- bind_rows(TissueSpecificity_boot_orth, TissueSpecificity_boot_SBGE_orth )

#Relative level of Nonsynonymous Diversity boot----
rlnd_boot_orth<-TwoBoot(SBGE_rlnd_bins_compare_orth, x_col = "piNfract", groupBy = "GeneType", boot_n = 10)
rlnd_boot_SBGE_orth <-TwoBoot_SBGE(SBGE_rlnd_bins_compare_orth, x_col = "piNfract", groupBy = "GeneType",SBGE_cat = "Sex.Bias.Category")

# Add "All" category to first dataframe
rlnd_boot_orth$Sex.Bias.Category <- "All"
names(rlnd_boot_SBGE_orth)[names(rlnd_boot_SBGE_orth) == "Sig"] <- "GeneType"

rlnd_boot_results_orth<- bind_rows(rlnd_boot_orth, rlnd_boot_SBGE_orth )

#ALL DISEASE----

#tajD boot----
tajD_boot_AD<-TwoBoot(SBGE_tajDS_bins_compare_AD, x_col = "tajD.S", groupBy = "GeneType", boot_n = 10)
tajD_boot_SBGE_AD <-TwoBoot_SBGE(SBGE_tajDS_bins_compare_AD, x_col = "tajD.S", groupBy = "GeneType",SBGE_cat = "Sex.Bias.Category")

# Add "All" category to first dataframe
tajD_boot_AD$Sex.Bias.Category <- "All"
names(tajD_boot_SBGE_AD)[names(tajD_boot_SBGE_AD) == "Sig"] <- "GeneType"

tajD_boot_results_AD<- bind_rows(tajD_boot_AD, tajD_boot_SBGE_AD )

#rMF boot----
rmf_boot_AD<-TwoBoot(SBGE_rMF_bins_compare_AD, x_col = "rmf", groupBy = "GeneType", boot_n = 10)
rmf_boot_SBGE_AD <-TwoBoot_SBGE(SBGE_rMF_bins_compare_AD, x_col = "rmf", groupBy = "GeneType",SBGE_cat = "Sex.Bias.Category")

# Add "All" category to first dataframe
rmf_boot_AD$Sex.Bias.Category <- "All"

names(rmf_boot_SBGE_AD)[names(rmf_boot_SBGE_AD) == "Sig"] <- "GeneType"


rmf_boot_results_AD<- bind_rows(rmf_boot_AD, rmf_boot_SBGE_AD )


#DoS boot----
DoS_boot_AD<-TwoBoot(SBGE_DoS_bins_compare_AD, x_col = "DoS", groupBy = "GeneType", boot_n = 10)
DoS_boot_SBGE_AD <-TwoBoot_SBGE(SBGE_DoS_bins_compare_AD, x_col = "DoS", groupBy = "GeneType",SBGE_cat = "Sex.Bias.Category")

# Add "All" category to first dataframe
DoS_boot_AD$Sex.Bias.Category <- "All"
names(DoS_boot_SBGE_AD)[names(DoS_boot_SBGE_AD) == "Sig"] <- "GeneType"

DoS_boot_results_AD<- bind_rows(DoS_boot_AD, DoS_boot_SBGE_AD )



#Tissue Specificity boot----
TissueSpecificity_boot_AD<-TwoBoot(SBGE_tissuespec_bins_compare_AD, x_col = "tissueSpecificityAdults", groupBy = "GeneType", boot_n = 10)
TissueSpecificity_boot_SBGE_AD <-TwoBoot_SBGE(SBGE_tissuespec_bins_compare_AD, x_col = "tissueSpecificityAdults", groupBy = "GeneType",SBGE_cat = "Sex.Bias.Category")

# Add "All" category to first dataframe
TissueSpecificity_boot_AD$Sex.Bias.Category <- "All"
names(TissueSpecificity_boot_SBGE_AD)[names(TissueSpecificity_boot_SBGE_AD) == "Sig"] <- "GeneType"

TissueSpecificity_boot_results_AD<- bind_rows(TissueSpecificity_boot_AD, TissueSpecificity_boot_SBGE_AD )

#Relative level of Nonsynonymous Diversity boot----
rlnd_boot_AD<-TwoBoot(SBGE_rlnd_bins_compare_AD, x_col = "piNfract", groupBy = "GeneType", boot_n = 10)
rlnd_boot_SBGE_AD <-TwoBoot_SBGE(SBGE_rlnd_bins_compare_AD, x_col = "piNfract", groupBy = "GeneType",SBGE_cat = "Sex.Bias.Category")

# Add "All" category to first dataframe
rlnd_boot_AD$Sex.Bias.Category <- "All"
names(rlnd_boot_SBGE_AD)[names(rlnd_boot_SBGE_AD) == "Sig"] <- "GeneType"

rlnd_boot_results_AD<- bind_rows(rlnd_boot_AD, rlnd_boot_SBGE_AD )

#_____________________----
#PLOTS----
#ALL 4 POINTS ON ONE PLOT----

# --- Tajimas D ----
tajD_cand_ref <- summary_tajDS_data_comparative_AD %>% filter(GeneType == "Candidate")
tajD_orth     <- summary_tajDS_data_comparative_orth %>% filter(GeneType == "Ortholog")
tajD_disease  <- summary_tajDS_data_comparative_AD %>% filter(GeneType == "Disease")
tajD_background <-summary_tajDS_data_comparative %>% filter(GeneType == "Background")

tajD_plot_4 <- bind_rows(tajD_cand_ref,tajD_orth, tajD_disease,tajD_background) %>%
  mutate(
    Sex.Bias.Category = factor(
      Sex.Bias.Category,
      levels = c("All","Highly female-biased","Moderately female-biased","Unbiased","Moderately male-biased","Highly male-biased")
    ),
    GeneType = factor(GeneType, levels = c("Candidate","Ortholog","Disease", "Background"))
  )

# --- Tajimas D Plot---
tajD_plot_merge <- ggplot(tajD_plot_4, 
                       aes(x = Sex.Bias.Category, y = q50, color = GeneType)) + 
  # Points and error bars
  geom_pointrange(aes(ymin = q05, ymax = q95), 
                  position = position_dodge(width = 0.5), 
                  size = 1.2)+
  geom_point(aes(color = GeneType),
             position = position_dodge(width = 0.5),
             size = 2) +  
  geom_errorbar(aes(ymin = q05, ymax = q95), 
                position = position_dodge(width = 0.5),
                width = 0.4, 
                size = 1) +  # Thicker error bars
  # Sample size labels
  geom_text(aes(y = q95+0.002, label = SampleSize, group = GeneType), 
            position = position_dodge(width = 0.5) ,
            color = "black", size = 4) +
  
  # Manual colors & shapes
  scale_color_manual(values = c("Disease" = "darkgreen", "Candidate" = "purple", "Ortholog"= "orange", "Background"= "blue")) +
  
  # Axis labels
  labs(x = "Sex Bias Category", y = expression("Tajima's " * italic(D)), color = "Gene Type") +
  
  # X-axis category labels
  scale_x_discrete(labels = c("All", "Highly FB", "Female", "Unbiased", "Male", "Highly MB")) +
  
  # Vertical dividers
  geom_vline(xintercept = 1.5, color = "black", size = 1.2) +
  geom_vline(xintercept = 2.5:5.5, color = "grey60", size = 0.5) +
  
  # Theme
  theme_classic(base_size = 16) +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 14),
    axis.text.y = element_text(size = 14),
    axis.title = element_text(size = 16, face = "bold"),
    legend.position = "top",
    legend.box = "horizontal",
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14, face = "bold"),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 1.2),
    panel.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(20, 20, 20, 20)
  )

print(tajD_plot_merge)

# --- rMF ----
rMF_cand_ref <- summary_rMF_data_comparative_AD %>% filter(GeneType == "Candidate")
rMF_orth     <- summary_rMF_data_comparative_orth %>% filter(GeneType == "Ortholog")
rMF_disease  <- summary_rMF_data_comparative_AD %>% filter(GeneType == "Disease")
rMF_background <-summary_rMF_data_comparative %>% filter(GeneType == "Background")

rMF_plot_4 <- bind_rows(rMF_cand_ref,rMF_orth, rMF_disease,rMF_background) %>%
  mutate(
    Sex.Bias.Category = factor(
      Sex.Bias.Category,
      levels = c("All","Highly female-biased","Moderately female-biased","Unbiased","Moderately male-biased","Highly male-biased")
    ),
    GeneType = factor(GeneType, levels = c("Candidate","Ortholog","Disease","Background"))
  )


# --- rMF plot ---
rMF_plot_merge <- ggplot(rMF_plot_4, 
                       aes(x = Sex.Bias.Category, y = q50, color = GeneType)) + 
  # Points and error bars
  geom_pointrange(aes(ymin = q05, ymax = q95), 
                  position = position_dodge(width = 0.5), 
                  size = 1.2)+
  geom_point(aes(color = GeneType),
             position = position_dodge(width = 0.5),
             size = 2) +  
  geom_errorbar(aes(ymin = q05, ymax = q95), 
                position = position_dodge(width = 0.5),
                width = 0.4, 
                size = 1) +  # Thicker error bars
  # Sample size labels
  geom_text(aes(y = q95+0.05, label = SampleSize, group = GeneType), 
            position = position_dodge(width = 0.5) ,
            color = "black", size = 4) +
  
  # Manual colors & shapes
  scale_color_manual(values = c("Disease" = "darkgreen", "Candidate" = "purple", "Ortholog"= "orange", "Background"= "blue")) +
  
  # Axis labels
  labs(x = "Sex Bias Category", y = expression(italic(r)[MF]), color = "Gene Type") +
  
  # X-axis category labels
  scale_x_discrete(labels = c("All", "Highly FB", "Female", "Unbiased", "Male", "Highly MB")) +
  
  # Vertical dividers
  geom_vline(xintercept = 1.5, color = "black", size = 1.2) +
  geom_vline(xintercept = 2.5:5.5, color = "grey60", size = 0.5) +
  
  # Theme
  theme_classic(base_size = 16) +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 14),
    axis.text.y = element_text(size = 14),
    axis.title = element_text(size = 16, face = "bold"),
    legend.position = "top",
    legend.box = "horizontal",
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14, face = "bold"),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 1.2),
    panel.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(20, 20, 20, 20)
  )

print(rMF_plot_merge)

# DoS ----
DoS_cand_ref <- summary_DoS_data_comparative_AD %>% filter(GeneType == "Candidate")
DoS_orth     <- summary_DoS_data_comparative_orth %>% filter(GeneType == "Ortholog")
DoS_disease  <- summary_DoS_data_comparative_AD %>% filter(GeneType == "Disease")
DoS_background <-summary_DoS_data_comparative %>% filter(GeneType == "Background")

DoS_plot_4 <- bind_rows(DoS_cand_ref,DoS_orth, DoS_disease,DoS_background) %>%
  mutate(
    Sex.Bias.Category = factor(
      Sex.Bias.Category,
      levels = c("All","Highly female-biased","Moderately female-biased","Unbiased","Moderately male-biased","Highly male-biased")
    ),
    GeneType = factor(GeneType, levels = c("Candidate","Ortholog","Disease","Background"))
  )

DoS_plot_merge <- ggplot(DoS_plot_4 , 
                       aes(x = Sex.Bias.Category, y = q50, color = GeneType)) + 
  # Points and error bars
  geom_pointrange(aes(ymin = q05, ymax = q95), 
                  position = position_dodge(width = 0.5), 
                  size = 1.2)+
  geom_point(aes(color = GeneType),
             position = position_dodge(width = 0.5),
             size = 2) +  
  geom_errorbar(aes(ymin = q05, ymax = q95), 
                position = position_dodge(width = 0.5),
                width = 0.4, 
                size = 1) +  # Thicker error bars
  # Sample size labels
  geom_text(aes(y = q95+0.01, label = SampleSize, group = GeneType), 
            position = position_dodge(width = 0.5) ,
            color = "black", size = 4) +
  
  # Manual colors & shapes
  scale_color_manual(values = c("Disease" = "darkgreen", "Candidate" = "purple","Ortholog"= "orange","Background"="blue")) +
  
  # Axis labels
  labs(x = "Sex Bias Category", y = expression(italic(DoS)), color = "Gene Type") +
  
  # X-axis category labels
  scale_x_discrete(labels = c("All", "Highly FB", "Female", "Unbiased", "Male", "Highly MB")) +
  
  # Vertical dividers
  geom_vline(xintercept = 1.5, color = "black", size = 1.2) +
  geom_vline(xintercept = 2.5:5.5, color = "grey60", size = 0.5) +
  
  # Theme
  theme_classic(base_size = 16) +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 14),
    axis.text.y = element_text(size = 14),
    axis.title = element_text(size = 16, face = "bold"),
    legend.position = "top",
    legend.box = "horizontal",
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14, face = "bold"),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 1.2),
    panel.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(20, 20, 20, 20)
  )

print(DoS_plot_merge)

# TissueSpecificty ----

tissuespec_cand_ref <- summary_tissuespec_data_comparative_AD %>% filter(GeneType == "Candidate")
tissuespec_orth     <- summary_tissuespec_data_comparative_orth %>% filter(GeneType == "Ortholog")
tissuespec_disease  <- summary_tissuespec_data_comparative_AD %>% filter(GeneType == "Disease")
tissuespec_background <-summary_tissuespec_data_comparative %>% filter(GeneType == "Background")

tissuespec_plot_4 <- bind_rows(tissuespec_cand_ref,tissuespec_orth, tissuespec_disease, tissuespec_background) %>%
  mutate(
    Sex.Bias.Category = factor(
      Sex.Bias.Category,
      levels = c("All","Highly female-biased","Moderately female-biased","Unbiased","Moderately male-biased","Highly male-biased")
    ),
    GeneType = factor(GeneType, levels = c("Candidate","Ortholog","Disease", "Background"))
  )

# --- tissuespec Plot---
tissuespec_plot_merge <- ggplot(tissuespec_plot_4, 
                          aes(x = Sex.Bias.Category, y = q50, color = GeneType)) + 
  # Points and error bars
  geom_pointrange(aes(ymin = q05, ymax = q95), 
                  position = position_dodge(width = 0.5), 
                  size = 1.2)+
  geom_point(aes(color = GeneType),
             position = position_dodge(width = 0.5),
             size = 2) +  
  geom_errorbar(aes(ymin = q05, ymax = q95), 
                position = position_dodge(width = 0.5),
                width = 0.4, 
                size = 1) +  # Thicker error bars
  # Sample size labels
  geom_text(aes(y = q95+0.002, label = SampleSize, group = GeneType), 
            position = position_dodge(width = 0.5) ,
            color = "black", size = 4) +
  
  # Manual colors & shapes
  scale_color_manual(values = c("Disease" = "darkgreen", "Candidate" = "purple", "Ortholog"= "orange", "Background"= "blue")) +
  
  # Axis labels
  labs(x = "Sex Bias Category", y = "Tissue Specificity", color = "Gene Type") +
  
  # X-axis category labels
  scale_x_discrete(labels = c("All", "Highly FB", "Female", "Unbiased", "Male", "Highly MB")) +
  
  # Vertical dividers
  geom_vline(xintercept = 1.5, color = "black", size = 1.2) +
  geom_vline(xintercept = 2.5:5.5, color = "grey60", size = 0.5) +
  
  # Theme
  theme_classic(base_size = 16) +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 14),
    axis.text.y = element_text(size = 14),
    axis.title = element_text(size = 16, face = "bold"),
    legend.position = "top",
    legend.box = "horizontal",
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14, face = "bold"),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 1.2),
    panel.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(20, 20, 20, 20)
  )

print(tissuespec_plot_merge)






tissuespec_cand_ref <- summary_tissuespec_data_comparative_AD %>% filter(GeneType == "Candidate")
tissuespec_orth     <- summary_tissuespec_data_comparative_orth %>% filter(GeneType == "Ortholog")
tissuespec_disease  <- summary_tissuespec_data_comparative_AD %>% filter(GeneType == "Disease")

tissuespec_plot <- bind_rows(tissuespec_cand_ref,tissuespec_orth, tissuespec_disease) %>%
  mutate(
    Sex.Bias.Category = factor(
      Sex.Bias.Category,
      levels = c("All","Highly female-biased","Moderately female-biased","Unbiased","Moderately male-biased","Highly male-biased")
    ),
    GeneType = factor(GeneType, levels = c("Candidate","Ortholog","Disease"))
  )

tissuespec_plot_merge <- ggplot(tissuespec_plot, 
                        aes(x = Sex.Bias.Category, y = q50, color = GeneType)) + 
  # Points and error bars
  geom_pointrange(aes(ymin = q05, ymax = q95), 
                  position = position_dodge(width = 0.5), 
                  size = 1.2)+
  geom_point(aes(color = GeneType),
             position = position_dodge(width = 0.5),
             size = 2) +  
  geom_errorbar(aes(ymin = q05, ymax = q95), 
                position = position_dodge(width = 0.5),
                width = 0.4, 
                size = 1) +  # Thicker error bars
  # Sample size labels
  geom_text(aes(y = q95+0.02, label = SampleSize, group = GeneType), 
            position = position_dodge(width = 0.5) ,
            color = "black", size = 4) +
  
  # Manual colors & shapes
  scale_color_manual(values = c("Disease" = "darkgreen", "Candidate" = "purple","Ortholog"= "orange")) +
  
  # Axis labels
  labs(x = "Sex Bias Category", y = "Tissue Specificity", color = "Gene Type") +
  
  # X-axis category labels
  scale_x_discrete(labels = c("All", "Highly FB", "Female", "Unbiased", "Male", "Highly MB")) +
  
  # Vertical dividers
  geom_vline(xintercept = 1.5, color = "black", size = 1.2) +
  geom_vline(xintercept = 2.5:5.5, color = "grey60", size = 0.5) +
  
  # Theme
  theme_classic(base_size = 16) +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 14),
    axis.text.y = element_text(size = 14),
    axis.title = element_text(size = 16, face = "bold"),
    legend.position = "top",
    legend.box = "horizontal",
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14, face = "bold"),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 1.2),
    panel.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(20, 20, 20, 20)
  )

print(tissuespec_plot_merge)


# rlnd ----
rlnd_cand_ref <- summary_rlnd_data_comparative_AD %>% filter(GeneType == "Candidate")
rlnd_orth     <- summary_rlnd_data_comparative_orth %>% filter(GeneType == "Ortholog")
rlnd_disease  <- summary_rlnd_data_comparative_AD %>% filter(GeneType == "Disease")
rlnd_background <-summary_rlnd_data_comparative %>% filter(GeneType == "Background")

rlnd_plot_4 <- bind_rows(rlnd_cand_ref,rlnd_orth, rlnd_disease,rlnd_background ) %>%
  mutate(
    Sex.Bias.Category = factor(
      Sex.Bias.Category,
      levels = c("All","Highly female-biased","Moderately female-biased","Unbiased","Moderately male-biased","Highly male-biased")
    ),
    GeneType = factor(GeneType, levels = c("Candidate","Ortholog","Disease","Background"))
  )


rlnd_plot_merge <- ggplot(rlnd_plot_4  , 
                       aes(x = Sex.Bias.Category, y = q50, color = GeneType)) + 
  # Points and error bars
  geom_pointrange(aes(ymin = q05, ymax = q95), 
                  position = position_dodge(width = 0.5), 
                  size = 1.2)+
  geom_point(aes(color = GeneType),
             position = position_dodge(width = 0.5),
             size = 2) +  
  geom_errorbar(aes(ymin = q05, ymax = q95), 
                position = position_dodge(width = 0.5),
                width = 0.4, 
                size = 1) +  # Thicker error bars
  # Sample size labels
  geom_text(aes(y = q95+0.006, label = SampleSize, group = GeneType), 
            position = position_dodge(width = 0.5) ,
            color = "black", size = 4) +
  
  # Manual colors & shapes
  scale_color_manual(values = c("Disease" = "darkgreen", "Candidate" = "purple", "Ortholog"="orange","Background"="blue")) +
  
  # Axis labels
  labs(x = "Sex Bias Category", y = expression(tilde(pi)[italic(r)*N]), color = "Gene Type") +
  
  # X-axis category labels
  scale_x_discrete(labels = c("All", "Highly FB", "Female", "Unbiased", "Male", "Highly MB")) +
  
  # Vertical dividers
  geom_vline(xintercept = 1.5, color = "black", size = 1.2) +
  geom_vline(xintercept = 2.5:5.5, color = "grey60", size = 0.5) +
  
  # Theme
  theme_classic(base_size = 16) +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 14),
    axis.text.y = element_text(size = 14),
    axis.title = element_text(size = 16, face = "bold"),
    legend.position = "top",
    legend.box = "horizontal",
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14, face = "bold"),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 1.2),
    panel.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(20, 20, 20, 20)
  )

print(rlnd_plot_merge)


#Combining 4 genomic signatures plots ----

(rMF_plot_merge / DoS_plot_merge /rlnd_plot_merge /tajD_plot_merge ) + plot_layout(guides = "collect") + plot_annotation(title = "My Combined Plot")

combined_plot_4 <- (rMF_plot_merge  / DoS_plot_merge / rlnd_plot_merge / tajD_plot_merge)+ 
  plot_layout(guides = "collect") + 
  plot_annotation(tag_levels = "A") & 
  theme(
    legend.position = "top",
    plot.tag = element_text(face = "bold", size = 14),
    plot.tag.position = c(0, 1)
  )
print(combined_plot_4)


ggsave("combined_plot_4 .png", combined_plot_4 , width = 17, height = 24, dpi = 300)


# Plots with Background as dashed line ----

# Base plot with everything except Background
rMF_plot_merge_dash <- ggplot(
  subset(rMF_plot_4, GeneType != "Background"), 
  aes(x = Sex.Bias.Category, y = q50, color = GeneType)
) +
  # Points and error bars
  geom_pointrange(aes(ymin = q05, ymax = q95),
                  position = position_dodge(width = 0.5),
                  size = 1.2) +
  geom_point(position = position_dodge(width = 0.5), size = 2) +
  geom_errorbar(aes(ymin = q05, ymax = q95),
                position = position_dodge(width = 0.5),
                width = 0.4,
                size = 1) +
  # Sample sizes
  geom_text(aes(y = q95 + 0.006, label = SampleSize, group = GeneType),
            position = position_dodge(width = 0.5),
            color = "black", size = 4) +
  
  # Background gene as black dashed line across bins
  geom_segment(
    data = subset(rMF_plot_4, GeneType == "Background"),
    aes(x = as.numeric(Sex.Bias.Category) - 0.4,
        xend = as.numeric(Sex.Bias.Category) + 0.4,
        y = q50, yend = q50),
    color = "black", linetype = "longdash", size = 1.2,
    inherit.aes = FALSE
  ) +
  
  # Colours
  scale_color_manual(values = c("Disease" = "darkgreen",
                                "Candidate" = "purple",
                                "Ortholog" = "orange")) +
  
  # Labels and Theme
  labs(x = "Sex Bias Category",
       y = expression(italic(r)[MF]),
       color = "Gene Type") +
  scale_x_discrete(labels = c("All", "Highly FB", "Female", "Unbiased", "Male", "Highly MB")) +
  geom_vline(xintercept = 1.5, color = "black", size = 1.2) +
  geom_vline(xintercept = 2.5:5.5, color = "grey60", size = 0.5) +
  theme_classic(base_size = 16) +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 14),
    axis.text.y = element_text(size = 14),
    axis.title = element_text(size = 16, face = "bold"),
    legend.position = "top",
    legend.box = "horizontal",
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14, face = "bold"),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 1.2),
    panel.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(20, 20, 20, 20)
  )

print(rMF_plot_merge_dash)

# Base plot with everything except Background

DoS_plot_merge_dash <- ggplot(
  subset(DoS_plot_4, GeneType != "Background"), 
  aes(x = Sex.Bias.Category, y = q50, color = GeneType)
) +
  # Points and error bars
  geom_pointrange(aes(ymin = q05, ymax = q95),
                  position = position_dodge(width = 0.5),
                  size = 1.2) +
  geom_point(position = position_dodge(width = 0.5), size = 2) +
  geom_errorbar(aes(ymin = q05, ymax = q95),
                position = position_dodge(width = 0.5),
                width = 0.4,
                size = 1) +
  # Sample sizes
  geom_text(aes(y = q95 + 0.006, label = SampleSize, group = GeneType),
            position = position_dodge(width = 0.5),
            color = "black", size = 4) +
  
  # Background gene as black dashed line across bins
  geom_segment(
    data = subset(DoS_plot_4, GeneType == "Background"),
    aes(x = as.numeric(Sex.Bias.Category) - 0.4,
        xend = as.numeric(Sex.Bias.Category) + 0.4,
        y = q50, yend = q50),
    color = "black", linetype = "longdash", size = 1.2,
    inherit.aes = FALSE
  ) +
  
  # Colours
  scale_color_manual(values = c("Disease" = "darkgreen",
                                "Candidate" = "purple",
                                "Ortholog" = "orange")) +
  
  # Labels and Theme
  labs(x = "Sex Bias Category",
       y = expression(italic(DoS)),
       color = "Gene Type") +
  scale_x_discrete(labels = c("All", "Highly FB", "Female", "Unbiased", "Male", "Highly MB")) +
  geom_vline(xintercept = 1.5, color = "black", size = 1.2) +
  geom_vline(xintercept = 2.5:5.5, color = "grey60", size = 0.5) +
  theme_classic(base_size = 16) +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 14),
    axis.text.y = element_text(size = 14),
    axis.title = element_text(size = 16, face = "bold"),
    legend.position = "top",
    legend.box = "horizontal",
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14, face = "bold"),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 1.2),
    panel.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(20, 20, 20, 20)
  )

print(DoS_plot_merge_dash)

# Base plot with everything except Background
tajD_plot_merge_dash <- ggplot(
  subset(tajD_plot_4, GeneType != "Background"), 
  aes(x = Sex.Bias.Category, y = q50, color = GeneType)
) +
  # Points and error bars
  geom_pointrange(aes(ymin = q05, ymax = q95),
                  position = position_dodge(width = 0.5),
                  size = 1.2) +
  geom_point(position = position_dodge(width = 0.5), size = 2) +
  geom_errorbar(aes(ymin = q05, ymax = q95),
                position = position_dodge(width = 0.5),
                width = 0.4,
                size = 1) +
  # Sample sizes
  geom_text(aes(y = q95 + 0.006, label = SampleSize, group = GeneType),
            position = position_dodge(width = 0.5),
            color = "black", size = 4) +
  
  # Background gene as black dashed line across bins
  geom_segment(
    data = subset(tajD_plot_4, GeneType == "Background"),
    aes(x = as.numeric(Sex.Bias.Category) - 0.4,
        xend = as.numeric(Sex.Bias.Category) + 0.4,
        y = q50, yend = q50),
    color = "black", linetype = "longdash", size = 1.2,
    inherit.aes = FALSE
  ) +
  
  # Manual colors for the others
  scale_color_manual(values = c("Disease" = "darkgreen",
                                "Candidate" = "purple",
                                "Ortholog" = "orange")) +
  
  # Labels, theme, etc (keep as in your code)
  labs(x = "Sex Bias Category",
       y =  expression("Tajima's " * italic(D)),
       color = "Gene Type") +
  scale_x_discrete(labels = c("All", "Highly FB", "Female", "Unbiased", "Male", "Highly MB")) +
  geom_vline(xintercept = 1.5, color = "black", size = 1.2) +
  geom_vline(xintercept = 2.5:5.5, color = "grey60", size = 0.5) +
  theme_classic(base_size = 16) +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 14),
    axis.text.y = element_text(size = 14),
    axis.title = element_text(size = 16, face = "bold"),
    legend.position = "top",
    legend.box = "horizontal",
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14, face = "bold"),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 1.2),
    panel.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(20, 20, 20, 20)
  )

print(tajD_plot_merge_dash)

# Base plot with everything except Background
rlnd_plot_merge_dash <- ggplot(
  subset(rlnd_plot_4, GeneType != "Background"), 
  aes(x = Sex.Bias.Category, y = q50, color = GeneType)
) +
  # Points and error bars
  geom_pointrange(aes(ymin = q05, ymax = q95),
                  position = position_dodge(width = 0.5),
                  size = 1.2) +
  geom_point(position = position_dodge(width = 0.5), size = 2) +
  geom_errorbar(aes(ymin = q05, ymax = q95),
                position = position_dodge(width = 0.5),
                width = 0.4,
                size = 1) +
  # Sample sizes
  geom_text(aes(y = q95 + 0.006, label = SampleSize, group = GeneType),
            position = position_dodge(width = 0.5),
            color = "black", size = 4) +
  
  # Background gene as black dashed line across bins
  geom_segment(
    data = subset(rlnd_plot_4, GeneType == "Background"),
    aes(x = as.numeric(Sex.Bias.Category) - 0.4,
        xend = as.numeric(Sex.Bias.Category) + 0.4,
        y = q50, yend = q50),
    color = "black", linetype = "longdash", size = 1.2,
    inherit.aes = FALSE
  ) +
  
  # Manual colors for the others
  scale_color_manual(values = c("Disease" = "darkgreen",
                                "Candidate" = "purple",
                                "Ortholog" = "orange")) +
  
  # Labels, theme, etc (keep as in your code)
  labs(x = "Sex Bias Category",
       y = expression(tilde(pi)[italic(r)*N]),
       color = "Gene Type") +
  scale_x_discrete(labels = c("All", "Highly FB", "Female", "Unbiased", "Male", "Highly MB")) +
  geom_vline(xintercept = 1.5, color = "black", size = 1.2) +
  geom_vline(xintercept = 2.5:5.5, color = "grey60", size = 0.5) +
  theme_classic(base_size = 16) +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 14),
    axis.text.y = element_text(size = 14),
    axis.title = element_text(size = 16, face = "bold"),
    legend.position = "top",
    legend.box = "horizontal",
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14, face = "bold"),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 1.2),
    panel.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(20, 20, 20, 20)
  )

print(rlnd_plot_merge_dash )

# Combining 4 genomic signatures plots (background dashed line) ----

(rMF_plot_merge_dash / DoS_plot_merge_dash /rlnd_plot_merge_dash  /tajD_plot_merge_dash) + plot_layout(guides = "collect") + plot_annotation(title = "My Combined Plot")

combined_plot_4_dash<- (rMF_plot_merge_dash  / DoS_plot_merge_dash / rlnd_plot_merge_dash  / tajD_plot_merge_dash)+ 
  plot_layout(guides = "collect") + 
  plot_annotation(tag_levels = "A") & 
  theme(
    legend.position = "top",
    plot.tag = element_text(face = "bold", size = 14),
    plot.tag.position = c(0, 1)
  )
print(combined_plot_4_dash)


ggsave("combined_plot_4_dash .png", combined_plot_4_dash, width = 17, height = 24, dpi = 300)
