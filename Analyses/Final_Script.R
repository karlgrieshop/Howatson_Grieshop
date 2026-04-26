#__________________________----
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

pub_theme <- function(base_size = 14) {
  theme_classic(base_size = base_size) +
    theme(
      axis.text        = element_text(size = base_size, colour = "black"),
      axis.title       = element_text(size = base_size + 2, face = "bold"),
      legend.text      = element_text(size = base_size),
      legend.title     = element_text(size = base_size + 1, face = "bold"),
      legend.position  = "top",
      axis.line        = element_blank(),
      panel.border     = element_rect(colour = "black", fill = NA, linewidth = 0.5),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      plot.margin      = margin(12, 12, 12, 12)
    )
}

#__________________________----
#IMPORTING DATA (GENOMIC METRICS) ----

setwd("/Users/elliehowatson/Final_Year_Project/R/data/Sex specific diseases")

# FlyBase disease name search----

# Imported list of flybase gene ID text files from sex specific diseases

sexspecific_disease_files <- list.files(path = "/Users/elliehowatson/Final_Year_Project/R/data/Sex specific diseases", pattern = "*.txt")

# Compiling list of sex specific disease flybase ids

sexspecific_disease_data <- lapply(sexspecific_disease_files, function(file) {
  disease_genes <- read_lines(file)
  disease_name <- tools::file_path_sans_ext(basename(file))  # Extract disease name from the file name
  data.frame(FlyBaseID = disease_genes, disease = disease_name)
})

# Combining sex specific disease flybase geneID tsv into one dataframe

all_sexspecific_diseases_data <- bind_rows(sexspecific_disease_data)
colnames(all_sexspecific_diseases_data)[1] <- "FlyBaseID"
all_sexspecific_diseases_data %>% summarise(count = n_distinct(FlyBaseID))

# Other diseases (literature search) ----

PCOS_orthologs <- read.csv("/Users/elliehowatson/Final_Year_Project/R/data/diopt_results_PCOS_genes.csv")
colnames(PCOS_orthologs)[8] <- "FlyBaseID"

human_SA_orthologs <- read.csv("/Users/elliehowatson/Final_Year_Project/R/data/diopt_results.csv")
colnames(human_SA_orthologs)[8] <- "FlyBaseID"

CAD_orthologs <- read.csv("/Users/elliehowatson/Final_Year_Project/R/data/coronary_artery_disease_orthologs.csv")
colnames(CAD_orthologs)[8] <- "FlyBaseID"

osteoporosis_orthologs <- read.csv("/Users/elliehowatson/Final_Year_Project/R/data/osteoporosis_orthologs.csv")
colnames(osteoporosis_orthologs)[8] <- "FlyBaseID"

AD_orthologs <- read.csv("/Users/elliehowatson/Final_Year_Project/R/data/alzheimer_disease_orthologs.csv")
colnames(AD_orthologs)[8] <- "FlyBaseID"

parkinsons_orthologs <- read.csv("/Users/elliehowatson/Final_Year_Project/R/data/parkinsons_orthologs.csv")
colnames(parkinsons_orthologs)[8] <- "FlyBaseID"

# Imported dataset for Mishra et al. 2024 ----

# SBGE contrasts

SBG_whole_bodies <-read.csv("/Users/elliehowatson/Final_Year_Project/R/data/Sex.Biased.Genes.Whole.csv")
colnames(SBG_whole_bodies)[1] <- "FlyBaseID"

SBG_whole_bodies <- SBG_whole_bodies %>% 
  mutate(Tissue = "Whole Body")

# Select and rename desired columns of SBGE 
# NOTE: use select()+rename() not cbind(), which coerces numeric log2FoldChange to character
SBG_whole_bodies <- SBG_whole_bodies %>%
  select(FlyBaseID, Sex.Bias.Category, log2FoldChange, Tissue) %>%
  rename(exp_SBGE_ase = log2FoldChange)

# Imported dataset from Fraisse et al. 2018 File S7 -----

Fra <- read.delim(file="/Users/elliehowatson/Final_Year_Project/R/data/FileS7.txt", header=TRUE)
colnames(Fra)[2] <- "FlyBaseID"
Fra[sapply(Fra, is.infinite)] <- NA

# Imported dataset from Singh and Agrawal 2023 ----

Singh <- read.csv(file="/Users/elliehowatson/Final_Year_Project/R/data/SBGEandSSSdataForMBE.csv", header=TRUE)
colnames(Singh)[2] <- "FlyBaseID"

# Add disease name column

PCOS_orthologs$disease <- rep(c("PCOS"), dim(PCOS_orthologs)[1])
human_SA_orthologs$disease <- rep(c("Human_SA"), dim(human_SA_orthologs)[1])
CAD_orthologs$disease <- rep(c("CAD"), dim(CAD_orthologs)[1])
osteoporosis_orthologs$disease <- rep(c("Osteoporosis"), dim(osteoporosis_orthologs)[1])
AD_orthologs$disease <- rep(c("Alzheimer"), dim(AD_orthologs)[1])
parkinsons_orthologs$disease <- rep(c("Parkinsons"), dim(parkinsons_orthologs)[1])

# rbind those all together

candidate_sexspecific_disease_orthologs <- rbind(PCOS_orthologs, human_SA_orthologs, CAD_orthologs, osteoporosis_orthologs, AD_orthologs, parkinsons_orthologs)

# Check number of uniq FlyBaseIDs
# all_orthologs %>% summarise(count = n_distinct(FlyBaseID))

#IMPORTING DATA (ENRICHMENT TESTS) ----
#Imported Innocenti and Morrow 2010 ----
InnoMorrow  <- read.csv(file="/Users/elliehowatson/Final_Year_Project/R/data/Innocenti_and_Morrow_ANTAGONISTIC.csv", header=TRUE)

#Imported Ruzicka et al. 2019 ----
Ruzicka <- read.csv(file="/Users/elliehowatson/Final_Year_Project/R/data/Ruzicka_et_al_pbio.3000244.s008.csv", header=TRUE)

#Imported Mishra et al. AI data 2024 ----

MishraAI <- read.csv(file="/Users/elliehowatson/Final_Year_Project/R/data/Mishra_et_al_Genes.with.Sex.Reversed.AI .csv", header=TRUE)
colnames(MishraAI )[1] <- "FlyBaseID"

#Imported Dataset Grieshop et al. 2025 ----
Grieshop <- read.csv(file="/Users/elliehowatson/Final_Year_Project/R/data/DE_candidates.csv", header=TRUE)

#Imported Dataset SA candidates ----
SA_genes <- read.csv(file="/Users/elliehowatson/Final_Year_Project/R/data/gene_hits_snps_90.csv", header=TRUE)

#__________________________----
#CLEAN DATA ----

# Filtering for the most accurate orthologs ----

filtered_orthologs <-candidate_sexspecific_disease_orthologs%>%
  filter(
    `Rank` == "high",
    `Best.Score` == "Yes",
    `Best.Score.Reverse` == "Yes",
    `DIOPT.Score` >= 3)

sexspecific_disease_orthologs <- filtered_orthologs[,c("Human.Symbol","Fly.Symbol", "FlyBaseID", "disease")]

# Assigning gene name to FlyBase geneIDs ----

ensembl <- useMart("ENSEMBL_MART_ENSEMBL", dataset = "dmelanogaster_gene_ensembl")

FlyBaseID <-all_sexspecific_data$FlyBaseID

gene_info <- getBM(
  attributes = c("flybase_gene_id","external_gene_name"),
  filters = "flybase_gene_id",
  values = FlyBaseID,
  mart = ensembl
)

# Renaming variables to match ortholog dataset ----

names(gene_info )[names(gene_info ) == "external_gene_name"] <- "Fly.Symbol"
names(gene_info )[names(gene_info ) == "flybase_gene_id"] <- "FlyBaseID"
gene_info <- merge(gene_info, all_sexspecific_diseases_data, by = "FlyBaseID", all = TRUE) 

sexspecific_disease_orthologs <- sexspecific_disease_orthologs[, c("Fly.Symbol", "FlyBaseID", "disease")]
names(gene_info)[names(gene_info) == "old_name"] <- "new_name"


sexspecific_disease_gene_data <- rbind(sexspecific_disease_orthologs, gene_info)

# Retrieving FlyBaseIDs for Innocenti and Morrow (2010) because only has the gene name----

flybase <- useMart("metazoa_mart", dataset = "dmelanogaster_eg_gene", host = "https://metazoa.ensembl.org")

gene_symbols <-InnoMorrow$Symbol

# Get the mapping (FlyBaseID)

mapping <- getBM(attributes = c("external_gene_name", "flybase_gene_id"),
                 filters = "external_gene_name",
                 values = gene_symbols,
                 mart = flybase)

# Merge into Innocenti and Morrow (2010)

InnoMorrow <- InnoMorrow%>%
  left_join(mapping , by = c("Symbol" = "external_gene_name"))

# 'Symbol' equal to 'flybase_gene_id'
InnoMorrow$Symbol <- InnoMorrow$flybase_gene_id

# Remove the original 'flybase_gene_id' column
InnoMorrow$flybase_gene_id <- NULL

# Rename 'Symbol' to 'FlyBaseID'
names(InnoMorrow)[names(InnoMorrow) == "Symbol"] <- "FlyBaseID"


#________________________----
#DATA ANALYSIS----

#Chromosome locations from reference genome ----

all.genes <- read.delim(file="~/Final_Year_Project/R/data/GeneLists/all.genes.tsv", sep="\t", header=FALSE)
colnames(all.genes) = c("FlyBaseID")

Xchr <- read.delim(file="~/Final_Year_Project/R/data/GeneLists/X.chromosome.genes.tsv", sep="\t", header=FALSE)
colnames(Xchr) = c("FlyBaseID")
Xchr$Chr <- rep("X", dim(Xchr)[1])

Ychr <- read.delim(file="~/Final_Year_Project/R/data/GeneLists/Y.chromosome.genes.tsv", sep="\t", header=FALSE)
colnames(Ychr) = c("FlyBaseID")
Ychr$Chr <- rep("Y", dim(Ychr)[1])

chr2L <- read.delim(file="~/Final_Year_Project/R/data/GeneLists/2L.chromosome.genes.tsv", sep="\t", header=FALSE)
colnames(chr2L) = c("FlyBaseID")
chr2R <- read.delim(file="~/Final_Year_Project/R/data/GeneLists/2R.chromosome.genes.tsv", sep="\t", header=FALSE)
colnames(chr2R) = c("FlyBaseID")

#
chr2 <- rbind(chr2L, chr2R)
chr2$Chr <- rep("2", dim(chr2)[1])

chr3L <- read.delim(file="~/Final_Year_Project/R/data/GeneLists/3L.chromosome.genes.tsv", sep="\t", header=FALSE)
colnames(chr3L) = c("FlyBaseID")
chr3R <- read.delim(file="~/Final_Year_Project/R/data/GeneLists/3R.chromosome.genes.tsv", sep="\t", header=FALSE)
colnames(chr3R) = c("FlyBaseID")

#
chr3 <- rbind(chr3L, chr3R)
chr3$Chr <- rep("3", dim(chr3)[1])

Chrs <- rbind(Xchr, Ychr, chr2, chr3) # Not all genes; just X, Y, 2, and 3.
Chrs$Chr <- as.factor(Chrs$Chr)


Chrs_sexspecific_disease <- merge(Chrs,sexspecific_disease_gene_data)

#GO enrichment for candidate DE genes----

# Target set
write_clip(Chrs_sexspecific_disease $FlyBaseID)# You decide which genes to use as target set
# Background set
write_clip(Chrs$FlyBaseID) 

# Background set
# These were used here: https://biit.cs.ut.ee/gprofiler/gost
# Alternatively, use: https://bioinformatics.sdstate.edu/go/

# MERGING DISEASE DATASETS ----

# Merging disease df with Fraisse et al. (2023) ----

Fra_master_data <- merge(Fra, sexspecific_disease_gene_data, by = "FlyBaseID")

# Reorder columns to put 'disease' first
cols <- c("disease", setdiff(names(Fra_master_data), "disease"))
Fra_master_data <- Fra_master_data[, cols]

# Removes sex specific disease name to have individual gene
Fra_master_data  <- Fra_master_data  %>% distinct()

# Merging sex specific disease df with Mishra et al. (2024) ----

SBGE_disease_data <- SBG_whole_bodies %>%
  inner_join(sexspecific_disease_gene_data, by = "FlyBaseID")

# Rename columns 
SBGE_disease_data<- SBGE_disease_data%>% distinct()
colnames(SBGE_disease_data) <- c("exp_SBGE_ase","Sex.Bias.Category","FlyBaseID","Tissue","Fly.Symbol","Disease")

#Remove disease name column
SBGE_disease_data<- SBGE_disease_data[, names(SBGE_disease_data) != "Disease"]

#Remove disease name column
SBGE_disease_data<- SBGE_disease_data[, names(SBGE_disease_data) != "Fly.Symbol"]

#Remove duplicates
SBGE_disease_data  <- SBGE_disease_data  %>% distinct()

#saving df
write.csv(SBGE_disease_data, 
          "/Users/elliehowatson/Final_Year_Project/R/data/SBGE_disease_data.csv", 
          row.names = FALSE)

# Merging Mishra  et al. (2024) Sex biased genes and Fraisse et al. (2018) ----

# Inner join by "FlyBaseID"
mish_fra_SBGE<- merge(SBG_whole_bodies,Fra_master_data, by = "FlyBaseID")

# Removes disease name to have individual gene
mish_fra_SBGE$disease <- NULL

#saving df
write.csv(mish_fra_SBGE, 
          "/Users/elliehowatson/Final_Year_Project/R/data/mish_fra_SBGE.csv", 
          row.names = FALSE)

# Merging disease genes df with Singh and Agrawal (2023) ----

singh_merge<- sexspecific_disease_gene_data %>%
  inner_join(Singh, by = "FlyBaseID")

# Merge Mishra  et al. (2024) with Singh and Agrawal (2023) ----

mish_singh_SBG <- SBG_whole_bodies %>%
  inner_join(Singh, by = "FlyBaseID")

#saving df
write.csv(mish_singh_SBG, 
          "/Users/elliehowatson/Final_Year_Project/R/data/mish_singh_SBG.csv", 
          row.names = FALSE)

# GENERATING CONTROL GENES ----

# Create a control set by excluding disease genes

sexspecific_disease_gene_ids <- unique(Chrs_sexspecific_disease$FlyBaseID)

control_genes <- Chrs%>%
  filter(!(FlyBaseID %in% sexspecific_disease_gene_ids))

# Select FlyBaseID and add GeneType column for background genes
# BUG FIX: exclude disease genes so they are not double-counted as Background
background_genes <- data.frame(
  FlyBaseID = Fra$FlyBaseID[!Fra$FlyBaseID %in% sexspecific_disease_gene_ids],
  GeneType = "Background",
  stringsAsFactors = FALSE
)

# Get unique FlyBaseID and add GeneType column for disease genes
sexspecific_disease_genes <- data.frame(
  FlyBaseID = unique(Chrs_sexspecific_disease$FlyBaseID),
  GeneType = "Candidate",
  stringsAsFactors = FALSE
)

background_and_sexspecific_disease<- bind_rows(background_genes , sexspecific_disease_genes)

# Merging control and sex specific disease genes with Mishra  et al. (2024)
mishra_merged<- left_join(SBG_whole_bodies, background_and_sexspecific_disease, by = "FlyBaseID")

# Merging Control Genes with Singh and Agrawal ----

# Filter Singh genes that are in Fra
Singh_filtered <- Singh[Singh$FlyBaseID %in% Fra$FlyBaseID, ]

# Merge control genes with Singh
control_singh <- merge(control_genes, Singh_filtered, by = "FlyBaseID")

# Add control label and rename columns
SBGE_control_data <- SBG_whole_bodies
SBGE_control_data$GeneType <- "Background"
colnames(SBGE_control_data) <- c("exp_SBGE_ase", "Sex.Bias.Category", "FlyBaseID", "Tissue", "GeneType")

# Add disease label to disease dataset and remove fly symbol column
SBGE_disease_data$GeneType <- "Disease"

# Bind control and disease rows into one data frame
SBGE_Mishra_data <- rbind(SBGE_control_data,SBGE_disease_data)

# Merge Control SBGE Data (Mishra) with Control Genes ----

control_mish_SBGE <- merge(SBGE_control_data, control_genes, by = "FlyBaseID")
control_mish_SBGE_wb <- merge(SBG_whole_bodies, control_genes, by = "FlyBaseID")

# Combine Control Singh and Mishra Data ----

control_singh_mish_SBGE_wb <- merge(control_mish_SBGE_wb, control_singh, by = "FlyBaseID")

# Label and format the control dataset
control_singh_mish_SBGE_wb$GeneType <- "Background"
names(control_singh_mish_SBGE_wb)[names(control_singh_mish_SBGE_wb) == "log2FoldChange"] <- "exp_SBGE_ase"

# Remove duplicated chromosome columns if they exist
control_singh_mish_SBGE_wb <- control_singh_mish_SBGE_wb[, !names(control_singh_mish_SBGE_wb) %in% c("Chr.x", "Chr.y")]

# Align column order and bind with disease data
mish_singh_SBG<- mish_singh_SBG[ , c("FlyBaseID",setdiff(names(mish_singh_SBG),"FlyBaseID"))]
mish_singh_SBG$GeneType <- "Disease"
mish_singh_SBG$FlyBaseID <- mish_singh_SBG$FlyBaseID  # Ensure column order for bind
control_set_genes_mish_singh_SBGE_wb <- rbind(control_singh_mish_SBGE_wb, mish_singh_SBG)


# Merging Control Genes with Fraisse and Mishra ----

# Merge control genes with Fra
control_Fra <- merge(control_genes, Fra, by = "FlyBaseID")

# IMPORTANT: use whole-body SBGE only for tissue-specificity background.
# Using all 3 tissues (gonads/heads/whole bodies) here would inflate Background ~3x
# because each gene would appear once per tissue with a different exp_SBGE_ase value,
# surviving the downstream distinct() call and tripling the effective gene count.
control_mish_SBGE_wb <- merge(
  SBG_whole_bodies[SBG_whole_bodies$Tissue == "Whole Body", ],
  control_genes,
  by = "FlyBaseID"
)

# Merge Mishra and Fraisse control genes
control_Fra_mish_SBGE_wb <- merge(control_mish_SBGE_wb, control_Fra, by = "FlyBaseID")

# Validation: assert one row per gene (catch future multi-tissue contamination)
stopifnot(
  "control_Fra_mish_SBGE_wb has duplicate FlyBaseIDs — check that only whole-body SBGE was used" =
    nrow(control_Fra_mish_SBGE_wb) == length(unique(control_Fra_mish_SBGE_wb$FlyBaseID))
)

# Format and label the final Fra+Mishra control dataset
control_Fra_mish_SBGE_wb$GeneType <- "Background"
names(control_Fra_mish_SBGE_wb)[names(control_Fra_mish_SBGE_wb) == "log2FoldChange"] <- "exp_SBGE_ase"

# Remove duplicate chromosome columns if they exist
control_Fra_mish_SBGE_wb <- control_Fra_mish_SBGE_wb[, !names(control_Fra_mish_SBGE_wb) %in% c("Chr.x", "Chr.y")]

# Align and bind with disease dataset
mish_fra_SBGE <- mish_fra_SBGE[, names(mish_fra_SBGE) != "Fly.Symbol"]
mish_fra_SBGE$GeneType <- "Disease"
mish_fra_SBGE$FlyBaseID <- mish_fra_SBGE$FlyBaseID  # Ensure order for bind
# Deduplicate disease genes (multiple diseases per gene create duplicate rows)
control_set_genes_mish_Fra_SBGE_wb <- rbind(control_Fra_mish_SBGE_wb, mish_fra_SBGE)


# Tajimas D Summary table ----

# Change into numerical 
library(dplyr)
SBGE_tajDS_compare<- control_set_genes_mish_singh_SBGE_wb%>%
  select(FlyBaseID,Sex.Bias.Category,exp_SBGE_ase,tajD.S,GeneType)%>% 
  distinct() # Remove disease variable to avoid repeats

SBGE_tajDS_compare<- SBGE_tajDS_compare %>%
  mutate(exp_SBGE_ase = as.numeric(exp_SBGE_ase))

# Creating Bins

SBGE_tajDS_bins_compare<-SBGE_tajDS_compare %>%
  mutate(Sex.Bias.Category = case_when(
    exp_SBGE_ase >= 2 ~ "Highly male-biased",
    exp_SBGE_ase >= 0.5 & exp_SBGE_ase < 2 ~ "Moderately male-biased",
    exp_SBGE_ase >= -0.5 & exp_SBGE_ase < 0.5 ~ "Unbiased",
    exp_SBGE_ase >= -2 & exp_SBGE_ase < -0.5 ~ "Moderately female-biased",
    exp_SBGE_ase < -2 ~ "Highly female-biased"
  ))

SBGE_tajDS_bins_compare <- SBGE_tajDS_bins_compare[!is.na(SBGE_tajDS_bins_compare$tajD.S), ]

# Compute Mean and Standard Error (SE) 

# Calculate mean
mean_tajD <- aggregate(tajD.S ~ Sex.Bias.Category + GeneType, 
                       data = SBGE_tajDS_bins_compare, 
                       FUN = function(x) mean(x, na.rm = TRUE))

# Calculate standard error
se_tajD <- aggregate(tajD.S ~ Sex.Bias.Category + GeneType, 
                     data = SBGE_tajDS_bins_compare, 
                     FUN = function(x) sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x))))

# Combine into one data frame
summary_tajDS_data_comparative <- merge(mean_tajD, se_tajD, 
                                        by = c("Sex.Bias.Category", "GeneType"))

# Rename columns
colnames(summary_tajDS_data_comparative) <- c("Sex.Bias.Category", "GeneType", 
                                              "Mean_tajD.S", "SE_tajD.S")

sample_sizes_tadj <- as.data.frame(table(SBGE_tajDS_bins_compare$Sex.Bias.Category,
                                         SBGE_tajDS_bins_compare$GeneType))

# Rename columns
colnames(sample_sizes_tadj) <- c("Sex.Bias.Category", "GeneType", "SampleSize")

summary_tajDS_data_comparative <- left_join(summary_tajDS_data_comparative, sample_sizes_tadj, by = c("Sex.Bias.Category", "GeneType"))

# Then compute the "All" category stats

split_data<- split(summary_tajDS_data_comparative, summary_tajDS_data_comparative$GeneType)

all_category <- do.call(rbind, lapply(split_data, function(df) {
  data.frame(
    GeneType = unique(df$GeneType),
    Mean_tajD.S = mean(df$Mean_tajD.S, na.rm = TRUE),
    SE_tajD.S = sqrt(sum(df$SE_tajD.S^2, na.rm = TRUE) / nrow(df)),
    SampleSize = sum(df$SampleSize, na.rm = TRUE),
    Sex.Bias.Category = "All"
  )
}))

summary_tajDS_data_comparative <- bind_rows(summary_tajDS_data_comparative, all_category)

# rMF Summary table ----

# Change into numerical 

SBGE_rMF_compare<- control_set_genes_mish_singh_SBGE_wb%>%
  select(FlyBaseID,Sex.Bias.Category,exp_SBGE_ase,Tissue,rmf,GeneType)%>% 
  distinct() # Remove disease variable to avoid repeats

SBGE_rMF_compare<- SBGE_rMF_compare%>%
  mutate(exp_SBGE_ase = as.numeric(exp_SBGE_ase))

# Creating Bins

SBGE_rMF_bins_compare<-SBGE_rMF_compare%>%
  mutate(Sex.Bias.Category = case_when(
    exp_SBGE_ase >= 2 ~ "Highly male-biased",
    exp_SBGE_ase >= 0.5 & exp_SBGE_ase < 2 ~ "Moderately male-biased",
    exp_SBGE_ase >= -0.5 & exp_SBGE_ase < 0.5 ~ "Unbiased",
    exp_SBGE_ase >= -2 & exp_SBGE_ase < -0.5 ~ "Moderately female-biased",
    exp_SBGE_ase < -2 ~ "Highly female-biased"
  ))

SBGE_rMF_bins_compare <- SBGE_rMF_bins_compare[!is.na(SBGE_rMF_bins_compare$rmf), ]

# Compute Mean and Standard Error (SE) 

# Compute mean
mean_rmf <- aggregate(rmf ~ Sex.Bias.Category + GeneType, 
                      data = SBGE_rMF_bins_compare, 
                      FUN = function(x) mean(x, na.rm = TRUE))

# Compute standard error
se_rmf <- aggregate(rmf ~ Sex.Bias.Category + GeneType, 
                    data = SBGE_rMF_bins_compare, 
                    FUN = function(x) sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x))))

# Merge the two summaries
summary_rMF_data_comparative <- merge(mean_rmf, se_rmf, 
                                      by = c("Sex.Bias.Category", "GeneType"))

# Rename columns
colnames(summary_rMF_data_comparative) <- c("Sex.Bias.Category", "GeneType", 
                                            "Mean_rmf", "SE_rmf")

rMF_sample_sizes <- as.data.frame(table(SBGE_rMF_bins_compare$Sex.Bias.Category,
                                        SBGE_rMF_bins_compare$GeneType))

# Rename columns 
colnames(rMF_sample_sizes) <- c("Sex.Bias.Category", "GeneType", "SampleSize")
summary_rMF_data_comparative <- left_join(summary_rMF_data_comparative, rMF_sample_sizes, by = c("Sex.Bias.Category", "GeneType"))

# Split by GeneType
split_data_rmf <- split(summary_rMF_data_comparative, summary_rMF_data_comparative$GeneType)

# For each GeneType, compute "All" summary row
all_category_rmf <- do.call(rbind, lapply(split_data_rmf , function(df) {
  data.frame(
    GeneType = unique(df$GeneType),
    Mean_rmf = mean(df$Mean_rmf, na.rm = TRUE),
    SE_rmf = sqrt(sum(df$SE_rmf^2, na.rm = TRUE) / nrow(df)),
    SampleSize = sum(df$SampleSize, na.rm = TRUE),
    Sex.Bias.Category = "All"
  )
}))

summary_rMF_data_comparative <- bind_rows(summary_rMF_data_comparative,all_category_rmf)

# DoS Summary table ----
# Change into numerical 

SBGE_DoS_compare<- control_set_genes_mish_Fra_SBGE_wb%>%
  select(FlyBaseID,Sex.Bias.Category,exp_SBGE_ase,DoS,GeneType)%>% 
  distinct() # Remove disease variable to avoid repeats

SBGE_DoS_compare<- SBGE_DoS_compare%>%
  mutate(exp_SBGE_ase = as.numeric(exp_SBGE_ase))

# Creating Bins

SBGE_DoS_bins_compare<-SBGE_DoS_compare%>%
  mutate(Sex.Bias.Category = case_when(
    exp_SBGE_ase >= 2 ~ "Highly male-biased",
    exp_SBGE_ase >= 0.5 & exp_SBGE_ase < 2 ~ "Moderately male-biased",
    exp_SBGE_ase >= -0.5 & exp_SBGE_ase < 0.5 ~ "Unbiased",
    exp_SBGE_ase >= -2 & exp_SBGE_ase < -0.5 ~ "Moderately female-biased",
    exp_SBGE_ase < -2 ~ "Highly female-biased"
  ))

# Compute Mean and Standard Error (SE) 

# Calculate mean
mean_DoS <- aggregate(DoS ~ Sex.Bias.Category + GeneType, 
                      data = SBGE_DoS_bins_compare, 
                      FUN = function(x) mean(x, na.rm = TRUE))

# Calculate standard error
se_DoS <- aggregate(DoS ~ Sex.Bias.Category + GeneType, 
                    data = SBGE_DoS_bins_compare, 
                    FUN = function(x) sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x))))

# Merge both summaries
summary_DoS_data_comparative <- merge(mean_DoS, se_DoS,
                                      by = c("Sex.Bias.Category", "GeneType"))

# Rename columns
colnames(summary_DoS_data_comparative) <- c("Sex.Bias.Category", "GeneType", "Mean_DoS", "SE_DoS")


DoS_sample_sizes <- as.data.frame(table(SBGE_DoS_bins_compare$Sex.Bias.Category,
                                        SBGE_DoS_bins_compare$GeneType))

# Rename columns
colnames(DoS_sample_sizes) <- c("Sex.Bias.Category", "GeneType", "SampleSize")

summary_DoS_data_comparative <- left_join(summary_DoS_data_comparative, DoS_sample_sizes, by = c("Sex.Bias.Category", "GeneType"))


# Split by GeneType
split_data_DoS <- split(summary_DoS_data_comparative, summary_DoS_data_comparative$GeneType)

# For each GeneType, compute "All" summary row
all_category_DoS <- do.call(rbind, lapply(split_data_DoS , function(df) {
  data.frame(
    GeneType = unique(df$GeneType),
    Mean_DoS = mean(df$Mean_DoS, na.rm = TRUE),
    SE_DoS = sqrt(sum(df$SE_DoS^2, na.rm = TRUE) / nrow(df)),
    SampleSize = sum(df$SampleSize, na.rm = TRUE),
    Sex.Bias.Category = "All"
  )
}))

# Reorder columns to match the expected format
all_category_DoS<- all_category_DoS [, c("Sex.Bias.Category", "GeneType", "Mean_DoS", "SE_DoS", "SampleSize")]


summary_DoS_data_comparative <- bind_rows(summary_DoS_data_comparative,all_category_DoS)

# Tissue Specificity Summary table  ----

# Change into numerical 

SBGE_tissuespec_compare <- control_set_genes_mish_Fra_SBGE_wb %>%
  select(FlyBaseID, Sex.Bias.Category, exp_SBGE_ase, tissueSpecificityAdults, GeneType) %>%
  distinct() # Collapse any remaining per-disease duplicates in disease genes

SBGE_tissuespec_compare<- SBGE_tissuespec_compare%>%
  mutate(exp_SBGE_ase = as.numeric(exp_SBGE_ase))

# Creating Bins

SBGE_tissuespec_bins_compare<-SBGE_tissuespec_compare%>%
  mutate(Sex.Bias.Category = case_when(
    exp_SBGE_ase >= 2 ~ "Highly male-biased",
    exp_SBGE_ase >= 0.5 & exp_SBGE_ase < 2 ~ "Moderately male-biased",
    exp_SBGE_ase >= -0.5 & exp_SBGE_ase < 0.5 ~ "Unbiased",
    exp_SBGE_ase >= -2 & exp_SBGE_ase < -0.5 ~ "Moderately female-biased",
    exp_SBGE_ase < -2 ~ "Highly female-biased"
  ))

SBGE_tissuespec_bins_compare <- SBGE_tissuespec_bins_compare[!is.na(SBGE_tissuespec_bins_compare$tissueSpecificityAdults), ]

# Compute Mean and Standard Error (SE) 

# Calculate mean
mean_tissueSpecificityAdults<- aggregate(tissueSpecificityAdults~ Sex.Bias.Category + GeneType, 
                                         data = SBGE_tissuespec_bins_compare, 
                                         FUN = function(x) mean(x, na.rm = TRUE))

# Calculate standard error
se_tissueSpecificityAdults<- aggregate(tissueSpecificityAdults ~ Sex.Bias.Category + GeneType, 
                                       data = SBGE_tissuespec_bins_compare, 
                                       FUN = function(x) sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x))))

# Merge both summaries
summary_tissuespec_data_comparative <- merge(mean_tissueSpecificityAdults,se_tissueSpecificityAdults,
                                             by = c("Sex.Bias.Category", "GeneType"))

# Rename columns
colnames(summary_tissuespec_data_comparative) <- c("Sex.Bias.Category", "GeneType", "Mean_tissueSpecificityAdults", "SE_tissueSpecificityAdults")


tissuespec_sample_sizes <- as.data.frame(table(SBGE_tissuespec_bins_compare$Sex.Bias.Category,
                                               SBGE_tissuespec_bins_compare$GeneType))

# Rename columns
colnames(tissuespec_sample_sizes) <- c("Sex.Bias.Category", "GeneType", "SampleSize")

summary_tissuespec_data_comparative <- left_join(summary_tissuespec_data_comparative, tissuespec_sample_sizes, by = c("Sex.Bias.Category", "GeneType"))


# Split by GeneType
split_data_tissuespec <- split(summary_tissuespec_data_comparative, summary_tissuespec_data_comparative$GeneType)

# For each GeneType, compute "All" summary row
all_category_tissuespec <- do.call(rbind, lapply(split_data_tissuespec , function(df) {
  data.frame(
    GeneType = unique(df$GeneType),
    Mean_tissueSpecificityAdults = mean(df$Mean_tissueSpecificityAdults, na.rm = TRUE),
    SE_tissueSpecificityAdults = sqrt(sum(df$SE_tissueSpecificityAdults^2, na.rm = TRUE) / nrow(df)),
    SampleSize = sum(df$SampleSize, na.rm = TRUE),
    Sex.Bias.Category = "All"
  )
}))

# Reorder columns to match the expected format
all_category_tissuespec <- all_category_tissuespec[, c("Sex.Bias.Category", "GeneType", "Mean_tissueSpecificityAdults", "SE_tissueSpecificityAdults", "SampleSize")]


summary_tissuespec_data_comparative <- bind_rows(summary_tissuespec_data_comparative,all_category_tissuespec)

# Relative Level of Nonsynonymous Diversity Summary table  ----

# Change into numerical 

SBGE_rlnd_compare<- control_set_genes_mish_singh_SBGE_wb%>%
  select(FlyBaseID,Sex.Bias.Category,exp_SBGE_ase,piNfract,GeneType)%>% 
  distinct() # Remove disease variable to avoid repeats

SBGE_rlnd_compare<- SBGE_rlnd_compare%>%
  mutate(exp_SBGE_ase = as.numeric(exp_SBGE_ase))

# Creating Bins

SBGE_rlnd_bins_compare<-SBGE_rlnd_compare%>%
  mutate(Sex.Bias.Category = case_when(
    exp_SBGE_ase >= 2 ~ "Highly male-biased",
    exp_SBGE_ase >= 0.5 & exp_SBGE_ase < 2 ~ "Moderately male-biased",
    exp_SBGE_ase >= -0.5 & exp_SBGE_ase < 0.5 ~ "Unbiased",
    exp_SBGE_ase >= -2 & exp_SBGE_ase < -0.5 ~ "Moderately female-biased",
    exp_SBGE_ase < -2 ~ "Highly female-biased"
  ))

SBGE_rlnd_bins_compare <- SBGE_rlnd_bins_compare[!is.na(SBGE_rlnd_bins_compare$piNfract), ]

# Compute Mean and Standard Error (SE) 

# Calculate mean
mean_rlnd<- aggregate(piNfract~ Sex.Bias.Category + GeneType, 
                      data = SBGE_rlnd_bins_compare, 
                      FUN = function(x) mean(x, na.rm = TRUE))

# Calculate standard error
se_rlnd<- aggregate(piNfract ~ Sex.Bias.Category + GeneType, 
                    data = SBGE_rlnd_bins_compare, 
                    FUN = function(x) sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x))))

# Merge both summaries
summary_rlnd_data_comparative <- merge(mean_rlnd,se_rlnd,
                                       by = c("Sex.Bias.Category", "GeneType"))

# Rename columns
colnames(summary_rlnd_data_comparative) <- c("Sex.Bias.Category", "GeneType", "Mean_rlnd", "SE_rlnd")


rlnd_sample_sizes <- as.data.frame(table(SBGE_rlnd_bins_compare$Sex.Bias.Category,
                                         SBGE_rlnd_bins_compare$GeneType))

# Rename columns
colnames(rlnd_sample_sizes) <- c("Sex.Bias.Category", "GeneType", "SampleSize")

summary_rlnd_data_comparative <- left_join(summary_rlnd_data_comparative, rlnd_sample_sizes, by = c("Sex.Bias.Category", "GeneType"))


# Split by GeneType
split_data_rlnd <- split(summary_rlnd_data_comparative, summary_rlnd_data_comparative$GeneType)

# For each GeneType, compute "All" summary row
all_category_rlnd <- do.call(rbind, lapply(split_data_rlnd , function(df) {
  data.frame(
    GeneType = unique(df$GeneType),
    Mean_rlnd = mean(df$Mean_rlnd, na.rm = TRUE),
    SE_rlnd = sqrt(sum(df$SE_rlnd^2, na.rm = TRUE) / nrow(df)),
    SampleSize = sum(df$SampleSize, na.rm = TRUE),
    Sex.Bias.Category = "All"
  )
}))

# Reorder columns to match the expected format

all_category_rlnd<- all_category_rlnd[, c("Sex.Bias.Category", "GeneType", "Mean_rlnd", "SE_rlnd", "SampleSize")]

summary_rlnd_data_comparative <- bind_rows(summary_rlnd_data_comparative,all_category_rlnd)


#___________________________-----
# STATISTICAL ANALYSES----

# Permutation tests----

# Adding logical variable enables function to work

summary_DoS_data_comparative$GeneType_logical <- summary_DoS_data_comparative$GeneType == "Background"
str(summary_DoS_data_comparative)

SBGE_DoS_bins_compare$GeneType_logical <- SBGE_DoS_bins_compare$GeneType == "Background"
str(SBGE_DoS_bins_compare)

summary_rMF_data_comparative$GeneType_logical <- summary_rMF_data_comparative$GeneType == "Background"
str(summary_rMF_data_comparative)

SBGE_rMF_bins_compare$GeneType_logical <- SBGE_rMF_bins_compare$GeneType == "Background"
str(SBGE_rMF_bins_compare)

summary_tajDS_data_comparative$GeneType_logical <- summary_tajDS_data_comparative$GeneType == "Background"
str(summary_tajDS_data_comparative)

SBGE_tajDS_bins_compare$GeneType_logical <- SBGE_tajDS_bins_compare$GeneType == "Background"
str(SBGE_tajDS_bins_compare)

summary_tissuespec_data_comparative$GeneType_logical <- summary_tissuespec_data_comparative$GeneType == "Background"
str(summary_tissuespec_data_comparative)

SBGE_tissuespec_bins_compare$GeneType_logical <- SBGE_tissuespec_bins_compare$GeneType == "Background"
str(SBGE_tissuespec_bins_compare)

summary_rlnd_data_comparative$GeneType_logical <- summary_rlnd_data_comparative$GeneType == "Background"
str(summary_rlnd_data_comparative)

SBGE_rlnd_bins_compare$GeneType_logical <- SBGE_rlnd_bins_compare$GeneType == "Background"
str(SBGE_rlnd_bins_compare)


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


#rMF perm----
# Now run TwoPerm
result_rmf<- TwoPerm(
  perm_dat = SBGE_rMF_bins_compare,
  x_col = "rmf",
  groupBy = "GeneType_logical"
)

print(result_rmf)

# Now run TwoPerm_SBGE
result_rmf_SBGE_cat <- TwoPerm_SBGE(
  perm_dat = SBGE_rMF_bins_compare,
  x_col = "rmf",
  groupBy = "GeneType_logical",
  SBGE_cat = "Sex.Bias.Category"
)

result_rmf_all <- data.frame(
  Sex.Bias.Category = "All",
  obs.diff = as.numeric(result_rmf[1]),
  n_TRUE = as.numeric(result_rmf[2]),
  n_FALSE = as.numeric(result_rmf[3]),
  pval = as.character(result_rmf[4]),
  Sig = as.character(result_rmf[5])
)

rmf_perm_results<- bind_rows(result_rmf_all, result_rmf_SBGE_cat)

#TajD perm----

# Now run TwoPerm
result_tajD <- TwoPerm(
  perm_dat = SBGE_tajDS_bins_compare,
  x_col = "tajD.S",
  groupBy = "GeneType_logical"
)

# Now run TwoPerm_SBGE() — it will include "All"
result_tajD_SBGE_cat <- TwoPerm_SBGE(
  perm_dat = SBGE_tajDS_bins_compare,
  x_col = "tajD.S",
  groupBy = "GeneType_logical",
  SBGE_cat = "Sex.Bias.Category"
)

result_tajD_all <- data.frame(
  Sex.Bias.Category = "All",
  obs.diff = as.numeric(result_tajD[1]),
  n_TRUE = as.numeric(result_tajD[2]),
  n_FALSE = as.numeric(result_tajD[3]),
  pval = as.character(result_tajD[4]),
  Sig = as.character(result_tajD[5])
)

tajD_perm_results<- bind_rows(result_tajD_all, result_tajD_SBGE_cat )

#DoS perm----

# Now run TwoPerm
result_DoS<- TwoPerm(
  perm_dat = SBGE_DoS_bins_compare,
  x_col = "DoS",
  groupBy = "GeneType_logical"
)

# Now run TwoPerm_SBGE() — it will include "All"
result_DoS_SBGE_cat <- TwoPerm_SBGE(
  perm_dat = SBGE_DoS_bins_compare,
  x_col = "DoS",
  groupBy = "GeneType_logical",
  SBGE_cat = "Sex.Bias.Category"
)

result_DoS_all <- data.frame(
  Sex.Bias.Category = "All",
  obs.diff = as.numeric(result_DoS[1]),
  n_TRUE = as.numeric(result_DoS[2]),
  n_FALSE = as.numeric(result_DoS[3]),
  pval = as.character(result_DoS[4]),
  Sig = as.character(result_DoS[5])
)

DoS_perm_results<- bind_rows(result_DoS_all, result_DoS_SBGE_cat )
print(result_DoS)

# Tissue Specificity perm----
# Now run TwoPerm
result_TissueSpecificity<- TwoPerm(
  perm_dat = SBGE_tissuespec_bins_compare,
  x_col = "tissueSpecificityAdults",
  groupBy = "GeneType_logical"
)

# Now run TwoPerm_SBGE() — it will include "All"
result_TissueSpecificity_SBGE_cat <- TwoPerm_SBGE(
  perm_dat = SBGE_tissuespec_bins_compare,
  x_col = "tissueSpecificityAdults",
  groupBy = "GeneType_logical",
  SBGE_cat = "Sex.Bias.Category"
)

result_TissueSpecificity_all <- data.frame(
  Sex.Bias.Category = "All",
  obs.diff = as.numeric(result_TissueSpecificity[1]),
  n_TRUE = as.numeric(result_TissueSpecificity[2]),
  n_FALSE = as.numeric(result_TissueSpecificity[3]),
  pval = as.character(result_TissueSpecificity[4]),
  Sig = as.character(result_TissueSpecificity[5])
)

TissueSpecificity_perm_results<- bind_rows(result_TissueSpecificity_all, result_TissueSpecificity_SBGE_cat )

# Relative level of Nonsynonymous Diversity----

# Now run TwoPerm
result_rlnd<- TwoPerm(
  perm_dat = SBGE_rlnd_bins_compare,
  x_col = "piNfract",
  groupBy = "GeneType_logical"
)

# Now run TwoPerm_SBGE() — it will include "All"
result_rlnd_SBGE_cat <- TwoPerm_SBGE(
  perm_dat = SBGE_rlnd_bins_compare,
  x_col = "piNfract",
  groupBy = "GeneType_logical",
  SBGE_cat = "Sex.Bias.Category"
)

result_rlnd_all <- data.frame(
  Sex.Bias.Category = "All",
  obs.diff = as.numeric(result_rlnd[1]),
  n_TRUE = as.numeric(result_rlnd[2]),
  n_FALSE = as.numeric(result_rlnd[3]),
  pval = as.character(result_rlnd[4]),
  Sig = as.character(result_rlnd[5])
)

rlnd_perm_results<- bind_rows(result_rlnd_all , result_rlnd_SBGE_cat )

# Bootstrapping functions ----

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




# testing

#rMF boot----
rmf_boot<-TwoBoot(SBGE_rMF_bins_compare, x_col = "rmf", groupBy = "GeneType", boot_n = 10)
rmf_boot_SBGE <-TwoBoot_SBGE(SBGE_rMF_bins_compare, x_col = "rmf", groupBy = "GeneType",SBGE_cat = "Sex.Bias.Category")

# Add "All" category to first dataframe
rmf_boot$Sex.Bias.Category <- "All"

names(rmf_boot_SBGE)[names(rmf_boot_SBGE) == "Sig"] <- "GeneType"


rmf_boot_results<- bind_rows(rmf_boot, rmf_boot_SBGE )

#DoS boot----
DoS_boot<-TwoBoot(SBGE_DoS_bins_compare, x_col = "DoS", groupBy = "GeneType", boot_n = 10)
DoS_boot_SBGE <-TwoBoot_SBGE(SBGE_DoS_bins_compare, x_col = "DoS", groupBy = "GeneType",SBGE_cat = "Sex.Bias.Category")

# Add "All" category to first dataframe
DoS_boot$Sex.Bias.Category <- "All"
names(DoS_boot_SBGE)[names(DoS_boot_SBGE) == "Sig"] <- "GeneType"

DoS_boot_results<- bind_rows(DoS_boot, DoS_boot_SBGE )

#tajD boot----
tajD_boot<-TwoBoot(SBGE_tajDS_bins_compare, x_col = "tajD.S", groupBy = "GeneType", boot_n = 10)
tajD_boot_SBGE <-TwoBoot_SBGE(SBGE_tajDS_bins_compare, x_col = "tajD.S", groupBy = "GeneType",SBGE_cat = "Sex.Bias.Category")

# Add "All" category to first dataframe
tajD_boot$Sex.Bias.Category <- "All"
names(tajD_boot_SBGE)[names(tajD_boot_SBGE) == "Sig"] <- "GeneType"

tajD_boot_results<- bind_rows(tajD_boot, tajD_boot_SBGE )

#Tissue Specificity boot----
TissueSpecificity_boot<-TwoBoot(SBGE_tissuespec_bins_compare, x_col = "tissueSpecificityAdults", groupBy = "GeneType", boot_n = 10)
TissueSpecificity_boot_SBGE <-TwoBoot_SBGE(SBGE_tissuespec_bins_compare, x_col = "tissueSpecificityAdults", groupBy = "GeneType",SBGE_cat = "Sex.Bias.Category")

# Add "All" category to first dataframe
TissueSpecificity_boot$Sex.Bias.Category <- "All"
names(TissueSpecificity_boot_SBGE)[names(TissueSpecificity_boot_SBGE) == "Sig"] <- "GeneType"

TissueSpecificity_boot_results<- bind_rows(TissueSpecificity_boot, TissueSpecificity_boot_SBGE )

#Relative level of Nonsynonymous Diversity boot----
rlnd_boot<-TwoBoot(SBGE_rlnd_bins_compare, x_col = "piNfract", groupBy = "GeneType", boot_n = 10)
rlnd_boot_SBGE <-TwoBoot_SBGE(SBGE_rlnd_bins_compare, x_col = "piNfract", groupBy = "GeneType",SBGE_cat = "Sex.Bias.Category")

# Add "All" category to first dataframe
rlnd_boot$Sex.Bias.Category <- "All"
names(rlnd_boot_SBGE)[names(rlnd_boot_SBGE) == "Sig"] <- "GeneType"

rlnd_boot_results<- bind_rows(rlnd_boot, rlnd_boot_SBGE )

# Enrichment Test for Candidate Sexually Antagonistic Genes----

library(biomaRt)

# Making background set using all ca 1300 protein coding genes in Drosophila
# Connect to Ensembl Metazoa and use the correct dataset
flymart <- useMart("metazoa_mart",
                   dataset = "dmelanogaster_eg_gene",
                   host = "https://metazoa.ensembl.org")

# Get all protein-coding genes with FlyBase ID, gene symbol, and chromosome
coding_genes <- getBM(attributes = c("flybase_gene_id", "external_gene_name", "chromosome_name"),
                      filters = "biotype",
                      values = "protein_coding",
                      mart = flymart)

coding_genes <- dplyr::rename(coding_genes,
                              FlyBaseID = flybase_gene_id,
                              GeneName = external_gene_name,
                              Chromosome = chromosome_name)


# Disease candidates

# Removing disease name and removing repeats
disease_chrs <- Chrs_disease[ , !(names(Chrs_disease) %in% c("disease", "Fly.Symbol"))]
disease_chrs <- unique(disease_chrs)

# Prepare dataset

disease_chrs<-disease_chrs %>%
  dplyr::mutate(Sig = TRUE)


coding_genes<-coding_genes%>%
  dplyr::mutate(Sig = FALSE)
test_df <- merge(disease_chrs[,c("FlyBaseID", "Sig")], 
                 coding_genes[,c("FlyBaseID", "Sig")], 
                 by ="FlyBaseID", all = T)

str(test_df)
test_df <- test_df[!(is.na(test_df$Sig.x) & is.na(test_df$Sig.y)),]
test_df <- test_df %>%
  dplyr::mutate(Sig = ifelse((!is.na(Sig.x) & Sig.x), TRUE,
                             ifelse(!is.na(Sig.y) & Sig.y, TRUE, FALSE)))

colnames(test_df) <- c("FlyBaseID", "Sig.Disease", "Sig.Control", "Sig")


# Then

# Get Ruzicka et al. (2019) SA candidates

colnames(Ruzicka) <- c("FlyBaseID", "FlyBaseName", "Chr", "SA_miss", # change accordingly
                       "SA_nonmiss", "NSA_miss", "NSA_nonmiss") # change accordingly
test_df <- test_df%>% 
  mutate(IsRuz = FlyBaseID %in% Ruzicka$FlyBaseID)



# Then 
# Get Innocenti & Morrow (2010) SA candidates
# combine InnoMorrow status with disease.genes dataset
test_df <- test_df %>%
  dplyr::mutate(IsInnoMorr = FlyBaseID %in% InnoMorrow$FlyBaseID) 



# Then
# Get Grieshop et al. (2025) SA candidates
# combine Grieshop status with disease.genes dataset
test_df <- test_df %>%
  dplyr::mutate(IsGrieshop = FlyBaseID %in% Grieshop$FlyBaseID) 


# Then
# Get Mishra et al. (2024) dominance-reversed genes
# combine Mishra status with disease.genes dataset
test_df  <- test_df  %>%
  dplyr::mutate(IsMishraAI = FlyBaseID %in% MishraAI$FlyBaseID) 


# Then
# Get New Sexually Antagonistic gene set 

colnames(SA_genes) <- c("chromosome_name","start_position","end_position","FlyBaseID","external_gene_name","gene_biotype") 

test_df  <- test_df  %>%
  dplyr::mutate(ISSA_genes = FlyBaseID %in% SA_genes$FlyBaseID) 


TestEnrichment <- function(data, df1 = "Sig", df2, xlab = "DE", ylab){
  test <- fisher.test(x = data[[df1]], y = data[[df2]])
  print(test)
  
  # kinda impractical but can't figure out how to 
  #     get ggbarstats to work with strings... alas
  colnames(data)[colnames(data) == paste0(df1)] = "df1"
  colnames(data)[colnames(data) == paste0(df2)] = "df2"
  
  mos_plot <- ggbarstats(
    data, df2, df1,
    results.subtitle = FALSE, # include fisher's exact test results
    subtitle = paste0(
      "Fisher's exact test", ", p-value = ",
      ifelse(test$p.value < 0.001, "< 0.001", round(test$p.value, 3))
    ),
    label = "both", # use both percentage and count as labels, 
    perc.k = 2, # include 2 decimal points in perc
    xlab = NULL
  ) +
    # some plot and theme settings
    scale_fill_manual(labels = c(paste0(ylab), paste0("not ",ylab)),
                      values = c("purple3", "darkgrey")) + 
    scale_x_discrete(labels = c(paste0("not",xlab), paste0(xlab))) +
    theme(plot.title.position = c("panel"),
          legend.title = element_blank(),
          legend.position = c("bottom"),
          legend.box.background = element_rect(),
          legend.text = element_text(size = 20, color = "black"),
          axis.text.x = element_text(size=20, margin = margin(5,0,0,0), color = "black"),
          axis.text.y = element_text(size=10, margin = margin(0,5,0,0), color = "black"),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size=40, margin = margin(0,10,0,0), color = "black"),
          plot.title = element_text(size=40, margin = margin(0,0,0,0), color = "black"),
          plot.margin = margin(6,6,6,6)
    )
  
  # get the contingency table
  plot.data <- mos_plot$data[,1:4]
  # rename columns appropriately
  colnames(plot.data)[colnames(plot.data) == "df1"] = paste0(xlab)
  colnames(plot.data)[colnames(plot.data) == "df2"] = paste0(df2)
  
  # assign new object in parent env
  assign(paste0(xlab,".vs.",ylab,".data"), plot.data, envir = parent.frame())
  
  # return the plot
  return(mos_plot)
}


# Run the enrichment tests----

# Testing for enrichment of disease genes in Innocenti Morrow candidates
disease.vs.Inn <- TestEnrichment(test_df, df2 = "IsInnoMorr", xlab = "Disease", ylab = "InnoMorr")
print(disease.vs.Inn)

# Testing for enrichment of disease genes in Ruzicka candidates
disease.vs.Ruz <- TestEnrichment(test_df, df2 = "IsRuz", xlab = "Disease", ylab = "Ruz")
print(disease.vs.Ruz)

# Testing for enrichment of disease genes in Grieshop candidates
disease.vs.Grieshop <- TestEnrichment(test_df, df2 = "IsGrieshop", xlab = "Disease", ylab = "Grieshop")
print(disease.vs.Grieshop)

# Testing for enrichment of disease genes in Mishra candidates
disease.vs.Mishra <- TestEnrichment(test_df, df2 = "IsMishraAI", xlab = "Disease", ylab = "MishraAI")
print(disease.vs.Mishra)

# Testing for enrichment of disease genes in ISSA candidates
disease.vs.SA_genes<- TestEnrichment(test_df, df2 = "ISSA_genes", xlab = "Disease", ylab = "SA_genes")
print(disease.vs.SA_genes)
library(ggstatsplot)
#___________________----
#DATA VISUALS----

#Tajimas D Plot----
summary_tajDS_data_comparative <- left_join(summary_tajDS_data_comparative, tajD_boot_results, by = c("Sex.Bias.Category", "GeneType"))

summary_tajDS_data_comparative$Sex.Bias.Category <- factor(summary_tajDS_data_comparative$Sex.Bias.Category, 
                                                           levels = c("All","Highly female-biased", "Moderately female-biased", 
                                                                      "Unbiased", "Moderately male-biased", "Highly male-biased"))

# --- Tajimas D Plot---
tajD_plot <- ggplot(summary_tajDS_data_comparative, 
                    aes(x = Sex.Bias.Category, y = q50, color = GeneType)) + 
  geom_pointrange(aes(ymin = q05, ymax = q95), 
                  position = position_dodge(width = 0.5), 
                  linewidth = 1.2) +
  # Sample size labels
  geom_text(aes(y = q95+0.002, label = SampleSize, group = GeneType), 
            position = position_dodge(width = 0.5) ,
            color = "black", size = 4) +
  
  # Manual colors & shapes
  scale_color_manual(values = c("Control" = "#999999", "Disease" = "#0072B2")) +
  
  # Axis labels
  labs(x = "Sex Bias Category", y = expression("Tajima's " * italic(D)), color = "Gene Type") +
  
  # X-axis category labels
  scale_x_discrete(labels = c("All", "Highly FB", "Female", "Unbiased", "Male", "Highly MB")) +
  
  # Vertical dividers
  geom_vline(xintercept = 1.5, color = "black", linewidth = 1.2) +
  geom_vline(xintercept = 2.5:5.5, color = "grey60", linewidth = 0.5) +
  
  # Theme
  pub_theme() +
  theme(legend.position = "top")

ggsave("tajD_plot.png", width = 13, height = 6, dpi = 300)

# --- rMF Plot----
summary_rMF_data_comparative <- left_join(summary_rMF_data_comparative, rmf_boot_results, by = c("Sex.Bias.Category", "GeneType"))

summary_rMF_data_comparative$Sex.Bias.Category <- factor(summary_rMF_data_comparative$Sex.Bias.Category, 
                                                         levels = c("All","Highly female-biased", "Moderately female-biased", 
                                                                    "Unbiased", "Moderately male-biased", "Highly male-biased"))
rmf_plot <- ggplot(summary_rMF_data_comparative, 
                   aes(x = Sex.Bias.Category, y = q50, color = GeneType)) + 
  geom_pointrange(aes(ymin = q05, ymax = q95), 
                  position = position_dodge(width = 0.5), 
                  linewidth = 1.2) +
  # Sample size labels
  geom_text(aes(y = q95+0.02, label = SampleSize, group = GeneType), 
            position = position_dodge(width = 0.5) ,
            color = "black", size = 4) +
  
  # Manual colors & shapes
  scale_color_manual(values = c("Control" = "#999999", "Disease" = "#0072B2")) +
  
  # Axis labels
  labs(x = "Sex Bias Category", y = expression(italic(r)[mf]), color = "Gene Type") +
  
  # X-axis category labels
  scale_x_discrete(labels = c("All", "Highly FB", "Female", "Unbiased", "Male", "Highly MB")) +
  
  # Vertical dividers
  geom_vline(xintercept = 1.5, color = "black", linewidth = 1.2) +
  geom_vline(xintercept = 2.5:5.5, color = "grey60", linewidth = 0.5) +
  
  # Theme
  pub_theme() +
  theme(legend.position = "top")

ggsave("rmf_plot.png", width = 13, height = 6, dpi = 300)

# DoS Plot----

summary_DoS_data_comparative <- left_join(summary_DoS_data_comparative, DoS_boot_results, by = c("Sex.Bias.Category", "GeneType"))

summary_DoS_data_comparative$Sex.Bias.Category <- factor(summary_DoS_data_comparative$Sex.Bias.Category, 
                                                         levels = c("All","Highly female-biased", "Moderately female-biased", 
                                                                    "Unbiased", "Moderately male-biased", "Highly male-biased"))


DoS_plot <- ggplot(summary_DoS_data_comparative, 
                   aes(x = Sex.Bias.Category, y = q50, color = GeneType)) + 
  geom_pointrange(aes(ymin = q05, ymax = q95), 
                  position = position_dodge(width = 0.5), 
                  linewidth = 1.2) +
  # Sample size labels
  geom_text(aes(y = q95+0.01, label = SampleSize, group = GeneType), 
            position = position_dodge(width = 0.5) ,
            color = "black", size = 4) +
  
  # Manual colors & shapes
  scale_color_manual(values = c("Control" = "#999999", "Disease" = "#0072B2")) +
  
  # Axis labels
  labs(x = "Sex Bias Category", y = expression(italic(DoS)), color = "Gene Type") +
  
  # X-axis category labels
  scale_x_discrete(labels = c("All", "Highly FB", "Female", "Unbiased", "Male", "Highly MB")) +
  
  # Vertical dividers
  geom_vline(xintercept = 1.5, color = "black", linewidth = 1.2) +
  geom_vline(xintercept = 2.5:5.5, color = "grey60", linewidth = 0.5) +
  
  # Theme
  pub_theme() +
  theme(legend.position = "top")

ggsave("DoS_plot.png", width = 13, height = 6, dpi = 300)

# Tissue Specificity Plot ----
summary_tissuespec_data_comparative<- left_join(summary_tissuespec_data_comparative, TissueSpecificity_boot_results, by = c("Sex.Bias.Category", "GeneType"))

summary_tissuespec_data_comparative$Sex.Bias.Category <- factor(summary_tissuespec_data_comparative$Sex.Bias.Category, 
                                                                levels = c("All","Highly female-biased", "Moderately female-biased", 
                                                                           "Unbiased", "Moderately male-biased", "Highly male-biased"))

tissuespec <- ggplot(summary_tissuespec_data_comparative, 
                     aes(x = Sex.Bias.Category, y = q50, color = GeneType)) + 
  geom_pointrange(aes(ymin = q05, ymax = q95), 
                  position = position_dodge(width = 0.5), 
                  linewidth = 1.2) +
  # Sample size labels
  geom_text(aes(y = q95+0.02, label = SampleSize, group = GeneType), 
            position = position_dodge(width = 0.5) ,
            color = "black", size = 4) +
  
  # Manual colors & shapes
  scale_color_manual(values = c("Control" = "#999999", "Disease" = "#0072B2")) +
  
  # Axis labels
  labs(x = "Sex Bias Category", y = "Tissue Specificity", color = "Gene Type") +
  
  # X-axis category labels
  scale_x_discrete(labels = c("All", "Highly FB", "Female", "Unbiased", "Male", "Highly MB")) +
  
  # Vertical dividers
  geom_vline(xintercept = 1.5, color = "black", linewidth = 1.2) +
  geom_vline(xintercept = 2.5:5.5, color = "grey60", linewidth = 0.5) +
  
  # Theme
  pub_theme() +
  theme(legend.position = "top")

ggsave("tissuespec.png", width = 13, height = 6, dpi = 300)


#Relative Level of Nonsynonymous Diveristy Plot ----

summary_rlnd_data_comparative <- left_join(summary_rlnd_data_comparative, rlnd_boot_results, by = c("Sex.Bias.Category", "GeneType"))

summary_rlnd_data_comparative$Sex.Bias.Category <- factor(summary_rlnd_data_comparative$Sex.Bias.Category, 
                                                          levels = c("All","Highly female-biased", "Moderately female-biased", 
                                                                     "Unbiased", "Moderately male-biased", "Highly male-biased"))

rlnd_plot <- ggplot(summary_rlnd_data_comparative, 
                    aes(x = Sex.Bias.Category, y = q50, color = GeneType)) + 
  geom_pointrange(aes(ymin = q05, ymax = q95), 
                  position = position_dodge(width = 0.5), 
                  linewidth = 1.2) +
  # Sample size labels
  geom_text(aes(y = q95+0.006, label = SampleSize, group = GeneType), 
            position = position_dodge(width = 0.5) ,
            color = "black", size = 4) +
  
  # Manual colors & shapes
  scale_color_manual(values = c("Control" = "#999999", "Disease" = "#0072B2")) +
  
  # Axis labels
  labs(x = "Sex Bias Category", y = expression(tilde(pi)[italic(r)*N]), color = "Gene Type") +
  
  # X-axis category labels
  scale_x_discrete(labels = c("All", "Highly FB", "Female", "Unbiased", "Male", "Highly MB")) +
  
  # Vertical dividers
  geom_vline(xintercept = 1.5, color = "black", linewidth = 1.2) +
  geom_vline(xintercept = 2.5:5.5, color = "grey60", linewidth = 0.5) +
  
  # Theme
  pub_theme() +
  theme(legend.position = "top")

ggsave("rlnd_plot.png", width = 13, height = 6, dpi = 300)


#Combining 4 genomic signatures plots ----

(rmf_plot / DoS_plot / rlnd_plot/ tajD_plot) + plot_layout(guides = "collect") + plot_annotation(title = "My Combined Plot")

combined_plot <- (rmf_plot / DoS_plot / rlnd_plot/ tajD_plot)+ 
  plot_layout(guides = "collect") + 
  plot_annotation(tag_levels = "A") & 
  theme(
    legend.position = "top",
    plot.tag = element_text(face = "bold", size = 14),
    plot.tag.position = c(0, 1)
  )
print(combined_plot)


ggsave("combined_plot.png", combined_plot, width = 10, height = 17, dpi = 300)
