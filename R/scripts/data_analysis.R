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
#__________________________----
#IMPORT DATA ----

setwd("/Users/elliehowatson/Final_Year_Project/R/data/Sex specific diseases")

# Imported list of flybase gene ID text files from sex specific diseases

disease_files <- list.files(path = "/Users/elliehowatson/Final_Year_Project/R/data/Sex specific diseases", pattern = "*.txt")

# Compiling list of sex specific disease flybase ids

disease_data <- lapply(disease_files, function(file) {
  disease_genes <- read_lines(file)
  disease_name <- tools::file_path_sans_ext(basename(file))  # Extract disease name from the file name
  data.frame(FlyBaseID = disease_genes, disease = disease_name)
})

# Combining sex specific disease flybase geneID tsv into one dataframe

all_diseases_data <- bind_rows(disease_data)
colnames(all_diseases_data)[1] <- "FlyBaseID"
all_diseases_data %>% summarise(count = n_distinct(FlyBaseID))


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

SBG_gonads <-read.csv("/Users/elliehowatson/Final_Year_Project/R/data/Sex.Biased.Genes.Gonads.csv")
colnames(SBG_gonads)[1] <- "FlyBaseID"

SBG_gonads <- SBG_gonads %>% 
  mutate(Tissue = "Gonads")

SBG_heads <-read.csv("/Users/elliehowatson/Final_Year_Project/R/data/Sex.Biased.Genes.Heads.csv")
colnames(SBG_heads)[1] <- "FlyBaseID"

SBG_heads <- SBG_heads%>% 
  mutate(Tissue = "Heads")

# Bind SBG Dataframe Together

SBG_data <- bind_rows(SBG_gonads, SBG_heads, SBG_whole_bodies)

# Select and rename desired columns of SBGE 

SBG_data <- data.frame(cbind(SBG_data$log2FoldChange, SBG_data$Sex.Bias.Category, SBG_data$FlyBaseID, SBG_data$Tissue))
colnames(SBG_data) <- c("exp_SBGE_ase","Sex.Bias.Category","FlyBaseID","Tissue")

# # Allele specific contrasts
# 
# DGE_whole_bodies <- read.csv("/Users/elliehowatson/Final_Year_Project/R/data/DGE_whole_bodies.csv")
# 
# DGE_whole_bodies <- DGE_whole_bodies%>% 
#   mutate(Tissue = "Whole Body")
# 
# DGE_gonads <- read.csv("/Users/elliehowatson/Final_Year_Project/R/data/DifferentialGeneExpressionGonads.csv")
# 
# DGE_gonads <- DGE_gonads%>% 
#   mutate(Tissue = "Gonads")
# 
# DGE_heads <- read.csv("/Users/elliehowatson/Final_Year_Project/R/data/DifferentialGeneExpression.csv")
# 
# DGE_heads <- DGE_heads%>% 
#   mutate(Tissue = "Heads")
# 
# # Bind DGE Dataframe Together
# 
# DGE_data<- rbind(DGE_whole_bodies, DGE_gonads, DGE_heads)
# 
# # Select and rename desired columns of Differential gene expression
# 
# DGE_data <- data.frame(cbind(DGE_data$log2FoldChange, DGE_data$lfcSE, DGE_data$FlyBaseID, DGE_data$Tissue))
# colnames(DGE_data) <- c("exp_SBGE_ase", "se_SBGE_ase", "FlyBaseID", "Tissue")
# 
# # Convert columns to numeric
# DGE_data$exp_SBGE_ase <- as.numeric(DGE_data$exp_SBGE_ase)
# DGE_data$se_SBGE_ase <- as.numeric(DGE_data$se_SBGE_ase)
# str(DGE_data)

# Imported dataset from Fraisse et al. 2018 -----
Fra <- read.delim(file="/Users/elliehowatson/Final_Year_Project/R/data/FileS7.txt", header=TRUE)
colnames(Fra)[2] <- "FlyBaseID"
Fra[sapply(Fra, is.infinite)] <- NA

# Imported dataset from Singh and Agrawal 2023 ----
Singh <- read.csv(file="/Users/elliehowatson/Final_Year_Project/R/data/SBGEandSSSdataForMBE.csv", header=TRUE)
colnames(Singh)[2] <- "FlyBaseID"


# Imported dataset from Innocenti and Morrow 2010
# ant_genes <- read.csv("/Users/elliehowatson/Final_Year_Project/R/data/Innocenti Antagonistic genes.csv")
# female_fitness_genes<- read.csv("/Users/elliehowatson/Final_Year_Project/R/data/Innocenti female fitness ass. genes.csv")
# male_fitness_genes<- read.csv("/Users/elliehowatson/Final_Year_Project/R/data/Innocenti male fitness ass. genes.csv")


# Add disease column to each of 

PCOS_orthologs$disease <- rep(c("PCOS"), dim(PCOS_orthologs)[1])
human_SA_orthologs$disease <- rep(c("Human_SA"), dim(human_SA_orthologs)[1])
CAD_orthologs$disease <- rep(c("CAD"), dim(CAD_orthologs)[1])
osteoporosis_orthologs$disease <- rep(c("Osteoporosis"), dim(osteoporosis_orthologs)[1])
AD_orthologs$disease <- rep(c("Alzheimer"), dim(AD_orthologs)[1])
parkinsons_orthologs$disease <- rep(c("Parkinsons"), dim(parkinsons_orthologs)[1])

# rbind those all together

all_orthologs <- rbind(PCOS_orthologs, human_SA_orthologs, CAD_orthologs, osteoporosis_orthologs, AD_orthologs, parkinsons_orthologs)

# Check number of uniq FlyBaseIDs
# all_orthologs %>% summarise(count = n_distinct(FlyBaseID))

# Imported Innocenti and Morrow 2010 ----
InnoMorrow  <- read.csv(file="/Users/elliehowatson/Final_Year_Project/R/data/Innocenti_and_Morrow_ANTAGONISTIC.csv", header=TRUE)

#Imported Ruzicka et al. 2019 ----
Ruzicka <- read.csv(file="/Users/elliehowatson/Final_Year_Project/R/data/Ruzicka_et_al_pbio.3000244.s008.csv", header=TRUE)

#Imported Mishra et al. AI data 2024 ----

MishraAI <- read.csv(file="/Users/elliehowatson/Final_Year_Project/R/data/Mishra_et_al_Genes.with.Sex.Reversed.AI .csv", header=TRUE)
colnames(MishraAI )[1] <- "FlyBaseID"


#Imported Dataset Greishop 2025 et al.----
Grieshop <- read.csv(file="/Users/elliehowatson/Final_Year_Project/R/data/DE_candidates.csv", header=TRUE)



#__________________________----
#CLEAN DATA ----

# Filtering for the most accurate orthologs ----

filtered_orthologs <-all_orthologs %>%
  filter(
    `Rank` == "high",
    `Best.Score` == "Yes",
    `Best.Score.Reverse` == "Yes",
    `DIOPT.Score` >= 3)

orthologs <- filtered_orthologs[,c("Human.Symbol","Fly.Symbol", "FlyBaseID", "disease")]

# Assigning gene name to FlyBase geneIDs ----

ensembl <- useMart("ENSEMBL_MART_ENSEMBL", dataset = "dmelanogaster_gene_ensembl")

FlyBaseID <-all_diseases_data$FlyBaseID

gene_info <- getBM(
  attributes = c("flybase_gene_id","external_gene_name"),
  filters = "flybase_gene_id",
  values = FlyBaseID,
  mart = ensembl
)

# Renaming variables to match ortholog dataset ----

names(gene_info )[names(gene_info ) == "external_gene_name"] <- "Fly.Symbol"
names(gene_info )[names(gene_info ) == "flybase_gene_id"] <- "FlyBaseID"

gene_info <- merge(gene_info, all_diseases_data, by = "FlyBaseID", all = TRUE) 

# Combining flybase search genes with ortholog genes ----

gene_info_inc_human<- gene_info %>%
  mutate("Human.Symbol" = NA) %>%  # Add missing column with NA
  select("Human.Symbol" , "Fly.Symbol", "FlyBaseID", "disease")

final_data_inc_human_genes<- bind_rows(orthologs, gene_info_inc_human)

final_data<- bind_rows(orthologs, gene_info)%>%
  select( "Fly.Symbol", "FlyBaseID", "disease")

# Retrieving FlyBaseIDs for Innocenti and Morrow (2010) because only has the gene name----

flybase <- useMart("metazoa_mart", dataset = "dmelanogaster_eg_gene", host = "https://metazoa.ensembl.org")

# Assume your gene symbols are in a vector called gene_symbols
gene_symbols <-InnoMorrow$Symbol

# Get the mapping
mapping <- getBM(attributes = c("external_gene_name", "flybase_gene_id"),
                 filters = "external_gene_name",
                 values = gene_symbols,
                 mart = flybase)

# Merge into your original data
InnoMorrow <- InnoMorrow%>%
  left_join(mapping , by = c("Symbol" = "external_gene_name"))

InnoMorrow <- InnoMorrow %>%
  mutate(Symbol = flybase_gene_id) %>%
  select(-flybase_gene_id)  # optionally drop the extra column
names(InnoMorrow )[names(InnoMorrow ) == "Symbol"] <- "FlyBaseID"
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


Chrs_disease <- merge(Chrs,final_data)


#GO enrichment for candidate differentially expressed genes----

# Target set
write_clip(Chrs_disease $FlyBaseID)# You decide which genes to use as target set
# Background set
write_clip(Chrs$FlyBaseID) 

# Background set
# These were used here: https://biit.cs.ut.ee/gprofiler/gost
# Alternatively, use: https://bioinformatics.sdstate.edu/go/

# Chromosomal distribution ----

# Count occurrences of each chromosome

chr_counts <- as.data.frame(table(Chrs_disease$Chr))
colnames(chr_counts) <- c("Chr", "count")

# Merging data frames with Fraisse et al.2023 ----

Fra_master_data <- Fra %>%
  inner_join(final_data, by = "FlyBaseID")%>%
  select(disease, everything())  # Reorder columns

# Remove exact duplicates
Fra_master_data  <- Fra_master_data  %>% distinct()

# Merging data frames with Mishra et al. 2024 ----

# DGE_master_data <- DGE_data%>%
#   inner_join(final_data, by = "FlyBaseID")

SBGE_disease_data <- SBG_data %>%
  inner_join(final_data, by = "FlyBaseID")

SBGE_disease_data<- SBGE_disease_data%>% distinct()
colnames(SBGE_disease_data) <- c("exp_SBGE_ase","Sex.Bias.Category","FlyBaseID","Tissue","Fly.Symbol","Disease")

SBG_wb_data <- SBG_whole_bodies %>%
  inner_join(final_data, by = "FlyBaseID")

SBG_wb_data <- SBG_wb_data %>% distinct()
colnames(SBG_wb_data) <- c("FlyBaseID","Sex.Bias.Category","exp_SBGE_ase","Tissue","Fly.Symbol","disease")


# Removing repeats where disease listed twice on DGE AND SBGE datasets 

SBG_no_disease_column <- SBG_wb_data %>%
  select(FlyBaseID,Sex.Bias.Category,exp_SBGE_ase,Tissue)# Reorder columns

SBG_no_disease_column <- SBG_no_disease_column %>% distinct() # Remove repeats

SBGE_disease_data_no_disease_column <- SBGE_disease_data %>%
  select(FlyBaseID,Sex.Bias.Category,exp_SBGE_ase,Tissue)# Reorder columns

SBGE_disease_data_no_disease_column<- SBGE_disease_data_no_disease_column %>% distinct() # Remove repeats

SBGE_disease_data_no_disease_column$exp_SBGE_ase <- as.numeric(SBGE_disease_data_no_disease_column$exp_SBGE_ase)
# DGE_no_disease_column <- DGE_master_data  %>%
#   select(FlyBaseID,exp_SBGE_ase,Tissue,Fly.Symbol)# Reorder columns
# 
# DGE_no_disease_column <- DGE_no_disease_column %>% distinct() # Remove repeats

# Merged Mishra Sex biased genes and Fraisse et al. 2018 ----
str(SBG_wb_data)

mish_fra_SBG_wb<- SBG_wb_data %>%
  inner_join(Fra, by = "FlyBaseID")%>%
  select(disease, everything())  # Reorder columns

# Removing disease column and repeats on mishra Sex biased genes and Fraisse et al. 2018

mish_fra_SBG_no_disease<- mish_fra_SBG_wb  %>%
  select(-disease, everything()) # Reorder columns


# Merging Data Frames with Singh and Agrawal 2023 ----

singh_merge_inc_human<- final_data_inc_human_genes %>%
  inner_join(Singh, by = "FlyBaseID")

singh_merge <- final_data %>%
  inner_join(Singh, by = "FlyBaseID")

# Merge Mishra with Singh and Agrawal 2023 ----

mish_singh_SBG_wb <- SBG_wb_data %>%
  inner_join(Singh, by = "FlyBaseID")

# GENERATING CONTROL GENES ----

# Create a control set by excluding disease genes

disease_gene_ids <- unique(Chrs_disease$FlyBaseID)

control_genes <- Chrs%>%
  filter(!(FlyBaseID %in% disease_gene_ids))

# Background----

# BUG FIX: exclude disease genes so Background and Disease sets are mutually exclusive
background_genes <- Fra %>%
  filter(!FlyBaseID %in% disease_gene_ids) %>%
  select(FlyBaseID) %>%
  mutate(GeneType = "Background")

disease_genes <- Chrs_disease %>%
  distinct(FlyBaseID) %>%             
  mutate(GeneType = "Disease")

background_and_disease <- bind_rows(background_genes , disease_genes)

mishra_merged <- left_join(SBG_whole_bodies, background_and_disease, by = "FlyBaseID")

# Remove all duplicate rows
clean_control_genes <- unique(control_genes)

# Merging Control genes with Mishra, Singh and Agrawal ----

# Merging Singh and Agrawal

Singh_filtered <- Singh[(Singh$FlyBaseID %in% Fra$FlyBaseID),]

control_singh<- control_genes %>%
  inner_join(Singh_filtered, by = "FlyBaseID")

# Merged control and Mishra

SBG_data<- rbind(SBG_whole_bodies, SBG_gonads, SBG_heads)

SBGE_control_data<- SBG_data%>%
  mutate(GeneType = "Control")
colnames(SBGE_control_data) <- c("FlyBaseID","Sex.Bias.Category","exp_SBGE_ase","Tissue","GeneType")

SBGE_disease_data_no_disease_column<- SBGE_disease_data_no_disease_column%>%
  mutate(GeneType = "Disease")


SBGE_Mishra_data <- bind_rows(SBGE_control_data,SBGE_disease_data_no_disease_column)


# Merged control and Mishra

control_mish_SBGE<-SBGE_control_data%>%
  inner_join( control_genes, by = "FlyBaseID")

control_mish_SBGE_wb<-SBG_whole_bodies%>%
  inner_join( control_genes, by = "FlyBaseID")

control_set_mish_SBGE<- control_mish_SBGE


# Merging control Mishra genes and control Singh genes

control_singh_mish_SBGE_wb<-control_mish_SBGE_wb%>%
  inner_join(control_singh, by = "FlyBaseID")
# Add a label to distinguish between control and disease genes

control_singh_mish_SBGE_wb$GeneType<- "Control"
mish_singh_SBG_wb$GeneType <- "Disease"

colnames(control_singh_mish_SBGE_wb)[colnames(control_singh_mish_SBGE_wb) == "log2FoldChange"] <- "exp_SBGE_ase" #changing control dataset so its the same as disease for binding

control_singh_mish_SBGE_wb<- control_singh_mish_SBGE_wb  %>%
  select(-Chr.x,-Chr.y) # Remove 'Chr.x' and 'Chr.y'

# Add missing columns: Fly.Symbol and disease (fill with NA)

control_singh_mish_SBGE_wb$Fly.Symbol <- NA

control_singh_mish_SBGE_wb$disease <- NA

# Bind the datasets together

mish_singh_SBG_wb <- mish_singh_SBG_wb%>% relocate(FlyBaseID)

control_set_genes_mish_singh_SBGE_wb <- rbind(control_singh_mish_SBGE_wb,mish_singh_SBG_wb)

# Merging Control genes with Mishra and Fraisse ----

control_Fra<- control_genes %>%
  inner_join(Fra, by = "FlyBaseID")

# IMPORTANT: use whole-body SBGE only for tissue-specificity background.
# Multi-tissue data would inflate Background ~3x because each gene would appear
# once per tissue with a different exp_SBGE_ase value, surviving distinct().
control_mish_SBGE_wb <- SBG_whole_bodies %>%
  filter(Tissue == "Whole Body") %>%
  inner_join(control_genes, by = "FlyBaseID")

# Merging control Mishra genes and control Fraisse genes

control_Fra_mish_SBGE_wb<-control_mish_SBGE_wb%>%
  inner_join(control_Fra, by = "FlyBaseID")

# Validation: assert one row per gene (catches future multi-tissue contamination)
stopifnot(
  "control_Fra_mish_SBGE_wb has duplicate FlyBaseIDs" =
    nrow(control_Fra_mish_SBGE_wb) == length(unique(control_Fra_mish_SBGE_wb$FlyBaseID))
)

# Add a label to distinguish between control and disease genes

control_Fra_mish_SBGE_wb$GeneType<- "Control"
mish_fra_SBG_wb$GeneType <- "Disease"

colnames(control_Fra_mish_SBGE_wb)[colnames(control_Fra_mish_SBGE_wb) == "log2FoldChange"] <- "exp_SBGE_ase" #changing control dataset so its the same as disease for binding

control_Fra_mish_SBGE_wb<- control_Fra_mish_SBGE_wb  %>%
  select(-Chr.x,-Chr.y) # Remove 'Chr.x' and 'Chr.y'

# Add missing columns: Fly.Symbol and disease (fill with NA)

control_Fra_mish_SBGE_wb$Fly.Symbol <- NA

control_Fra_mish_SBGE_wb$disease <- NA

# Bind the datasets together

mish_fra_SBG_wb <- mish_fra_SBG_wb%>% relocate(FlyBaseID)

control_set_genes_mish_Fra_SBGE_wb <- rbind(control_Fra_mish_SBGE_wb,mish_fra_SBG_wb)

# Tajimas D Summary table ----

# Change into numerical 

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
  distinct() # Collapse per-disease duplicates; do not include 'disease' column here

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
  distinct() # Collapse per-disease duplicates; do not include 'disease' column here

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
  distinct() # Collapse per-disease duplicates; do not include 'disease' column here

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

# Making mishra dataset with control and disease genes wide ----
SBGE_Mishra_data_wide <- SBGE_Mishra_data%>%
  select(FlyBaseID, Tissue, exp_SBGE_ase, GeneType, Sex.Bias.Category) %>%
  pivot_wider(
    names_from = Tissue,
    values_from = exp_SBGE_ase,
    names_prefix = "log2FC_"
  )

#TISSUE----

#Gonads

# Change into numerical 

SBGE_gonads_compare<- SBGE_Mishra_data_wide %>%
  select(FlyBaseID,Sex.Bias.Category,`log2FC_Whole Body`,log2FC_Gonads,GeneType)%>% 
  distinct() # Remove disease variable to avoid repeats

SBGE_gonads_compare<- SBGE_gonads_compare%>%
  mutate(`log2FC_Whole Body` = as.numeric(`log2FC_Whole Body`),
         log2FC_Gonads= as.numeric(`log2FC_Gonads`))

# Creating Bins

SBGE_gonads_bins_compare<-SBGE_gonads_compare%>%
  mutate(Sex.Bias.Category = case_when(
    `log2FC_Whole Body` >= 2 ~ "Highly male-biased",
    `log2FC_Whole Body` >= 0.5 & `log2FC_Whole Body` < 2 ~ "Moderately male-biased",
    `log2FC_Whole Body`>= -0.5 & `log2FC_Whole Body` < 0.5 ~ "Unbiased",
    `log2FC_Whole Body` >= -2 & `log2FC_Whole Body` < -0.5 ~ "Moderately female-biased",
    `log2FC_Whole Body` < -2 ~ "Highly female-biased"
  ))

SBGE_gonads_bins_compare <- SBGE_gonads_bins_compare[!is.na(SBGE_gonads_bins_compare$log2FC_Gonads), ]
SBGE_gonads_bins_compare <- SBGE_gonads_bins_compare[!is.na(SBGE_gonads_bins_compare$`log2FC_Whole Body`), ]
SBGE_gonads_bins_compare <- SBGE_gonads_bins_compare[!is.na(SBGE_gonads_bins_compare$Sex.Bias.Category), ]
# Compute Mean and Standard Error (SE) 

# Calculate mean
mean_gonads<- aggregate(log2FC_Gonads~ Sex.Bias.Category + GeneType, 
                      data = SBGE_gonads_bins_compare, 
                      FUN = function(x) mean(x, na.rm = TRUE))

# Calculate standard error
se_gonads<- aggregate(log2FC_Gonads ~ Sex.Bias.Category + GeneType, 
                    data = SBGE_gonads_bins_compare, 
                    FUN = function(x) sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x))))

# Merge both summaries
summary_gonads_data_comparative <- merge(mean_gonads,se_gonads,
                                       by = c("Sex.Bias.Category", "GeneType"))

# Rename columns
colnames(summary_gonads_data_comparative) <- c("Sex.Bias.Category", "GeneType", "Mean_lFC", "SE_lFC")

gonads_sample_sizes <- as.data.frame(table(SBGE_gonads_bins_compare$Sex.Bias.Category,
                                         SBGE_gonads_bins_compare$GeneType))

# Rename columns
colnames(gonads_sample_sizes) <- c("Sex.Bias.Category", "GeneType", "SampleSize")

summary_gonads_data_comparative  <- left_join(summary_gonads_data_comparative , gonads_sample_sizes, by = c("Sex.Bias.Category", "GeneType"))


# Split by GeneType
split_data_gonads<- split(summary_gonads_data_comparative, summary_gonads_data_comparative$GeneType)

# For each GeneType, compute "All" summary row
all_category_gonads <- do.call(rbind, lapply(split_data_gonads, function(df) {
  data.frame(
    GeneType = unique(df$GeneType),
    Mean_lFC = mean(df$Mean_lFC, na.rm = TRUE),
    SE_lFC = sqrt(sum(df$SE_lFC^2, na.rm = TRUE) / nrow(df)),
    SampleSize = sum(df$SampleSize, na.rm = TRUE),
    Sex.Bias.Category = "All"
  )
}))

# Reorder columns to match the expected format

all_category_gonads<- all_category_gonads[, c("Sex.Bias.Category", "GeneType", "Mean_lFC", "SE_lFC", "SampleSize")]

summary_gonads_data_comparative <- bind_rows(summary_gonads_data_comparative,all_category_gonads)

# HEADS

# Change into numerical 

SBGE_heads_compare<- SBGE_Mishra_data_wide %>%
  select(FlyBaseID,Sex.Bias.Category,`log2FC_Whole Body`,log2FC_Heads,GeneType)%>% 
  distinct() # Remove disease variable to avoid repeats

SBGE_heads_compare<- SBGE_heads_compare%>%
  mutate(`log2FC_Whole Body` = as.numeric(`log2FC_Whole Body`),
         log2FC_Heads= as.numeric(`log2FC_Heads`))

# Creating Bins

SBGE_heads_bins_compare<-SBGE_heads_compare%>%
  mutate(Sex.Bias.Category = case_when(
    `log2FC_Whole Body` >= 2 ~ "Highly male-biased",
    `log2FC_Whole Body` >= 0.5 & `log2FC_Whole Body` < 2 ~ "Moderately male-biased",
    `log2FC_Whole Body`>= -0.5 & `log2FC_Whole Body` < 0.5 ~ "Unbiased",
    `log2FC_Whole Body` >= -2 & `log2FC_Whole Body` < -0.5 ~ "Moderately female-biased",
    `log2FC_Whole Body` < -2 ~ "Highly female-biased"
  ))

SBGE_heads_bins_compare <- SBGE_heads_bins_compare[!is.na(SBGE_heads_bins_compare$log2FC_Heads), ]
SBGE_heads_bins_compare <- SBGE_heads_bins_compare[!is.na(SBGE_heads_bins_compare$`log2FC_Whole Body`), ]
SBGE_heads_bins_compare<- SBGE_heads_bins_compare[!is.na(SBGE_heads_bins_compare$Sex.Bias.Category), ]

# Compute Mean and Standard Error (SE) 

# Calculate mean
mean_heads<- aggregate(log2FC_Heads~ Sex.Bias.Category + GeneType, 
                        data = SBGE_heads_bins_compare, 
                        FUN = function(x) mean(x, na.rm = TRUE))

# Calculate standard error
se_heads<- aggregate(log2FC_Heads ~ Sex.Bias.Category + GeneType, 
                      data = SBGE_heads_bins_compare, 
                      FUN = function(x) sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x))))

# Merge both summaries
summary_heads_data_comparative <- merge(mean_heads,se_heads,
                                         by = c("Sex.Bias.Category", "GeneType"))

# Rename columns
colnames(summary_heads_data_comparative) <- c("Sex.Bias.Category", "GeneType", "Mean_lFC", "SE_lFC")

heads_sample_sizes <- as.data.frame(table(SBGE_heads_bins_compare$Sex.Bias.Category,
                                           SBGE_heads_bins_compare$GeneType))

# Rename columns
colnames(heads_sample_sizes) <- c("Sex.Bias.Category", "GeneType", "SampleSize")

summary_heads_data_comparative  <- left_join(summary_heads_data_comparative , heads_sample_sizes, by = c("Sex.Bias.Category", "GeneType"))


# Split by GeneType
split_data_heads<- split(summary_heads_data_comparative, summary_heads_data_comparative$GeneType)

# For each GeneType, compute "All" summary row
all_category_heads <- do.call(rbind, lapply(split_data_heads, function(df) {
  data.frame(
    GeneType = unique(df$GeneType),
    Mean_lFC = mean(df$Mean_lFC, na.rm = TRUE),
    SE_lFC = sqrt(sum(df$SE_lFC^2, na.rm = TRUE) / nrow(df)),
    SampleSize = sum(df$SampleSize, na.rm = TRUE),
    Sex.Bias.Category = "All"
  )
}))

# Reorder columns to match the expected format

all_category_heads<- all_category_heads[, c("Sex.Bias.Category", "GeneType", "Mean_lFC", "SE_lFC", "SampleSize")]

summary_heads_data_comparative <- bind_rows(summary_heads_data_comparative,all_category_heads)


# # #__________________________----
# #EXPLORATORY ANALYSIS ----
# #  
# # # Removing repeats where disease listed twice 
# # 
# # SBG_no_disease_column <- mish_fra_SBG  %>%
# #   select(FlyBaseID,Sex.Bias.Category,exp_SBGE_ase,Tissue,DoS)  # Reorder columns
# # 
# # SBG_no_disease_column <- SBG_no_disease_column%>% distinct()
# # 
# # # Change into numerical 
# # str(SBG_no_disease_column$exp_SBGE_ase)
# # SBG_no_disease_column <- SBG_no_disease_column %>%
# #   mutate(exp_SBGE_ase = as.numeric(exp_SBGE_ase))
# # 
# # #Creating Bins 
# # bins <- SBG_no_disease_column %>%
# #   mutate(Sex.Bias.Category = case_when(
# #     exp_SBGE_ase >= 2 ~ "Highly male-biased",
# #     exp_SBGE_ase >= 0.5 & exp_SBGE_ase < 2 ~ "Moderately male-biased",
# #     exp_SBGE_ase >= -0.5 & exp_SBGE_ase < 0.5 ~ "Unbiased",
# #     exp_SBGE_ase > -2 & exp_SBGE_ase <= -0.5 ~ "Moderately female-biased",
# #     exp_SBGE_ase <= -2 ~ "Highly female-biased"
# #   ))
# # 
# # 
# # color_map <- c("Highly female-biased" = "red", "Moderately female-biased" = "pink", 
# #                "Unbiased" = "black", "Moderately male-biased" = "lightblue", "Highly male-biased" = "blue")
# # 
# # bins$Sex.Bias.Category <- factor(bins $Sex.Bias.Category, 
# #                                  levels = c("Highly female-biased", "Moderately female-biased", 
# #                                             "Unbiased", "Moderately male-biased", "Highly male-biased"))
# # 
# # 
# # # Compute Mean and Standard Error (SE) for DoS per Sex Bias Category
# # summary_data <- bins %>%
# #   group_by(Sex.Bias.Category) %>%
# #   summarise(
# #     Mean_DoS = mean(DoS, na.rm = TRUE),
# #     SE_DoS = sd(DoS, na.rm = TRUE) / sqrt(n())  # Standard Error
# #   )
# # 
# # # Creating Summary table for DoS tissue data
# # summary_tissue_data <- bins %>%
# #   group_by(Sex.Bias.Category, Tissue) %>%
# #   summarise(
# #     Mean_DoS = mean(DoS, na.rm = TRUE),
# #     SE_DoS = sd(DoS, na.rm = TRUE) / sqrt(n()),  # Standard Error
# #     .groups = "drop"
# #   )
# # 
# # # Sex Bias/DoS 
# # 
# # ggplot(summary_data, aes(x = Sex.Bias.Category, y = Mean_DoS, color = Sex.Bias.Category)) +
# #   geom_point(size = 4) +  # Plot mean DoS
# #   geom_errorbar(aes(ymin = Mean_DoS - SE_DoS, ymax = Mean_DoS + SE_DoS), width = 0.2) +  # Error bars
# #   scale_color_manual(values = color_map) +  # Custom colors
# #   labs(title = "", x = "", y = expression(italic(DoS)))+
# #   
# #   theme_bw()+
# #   theme(
# #     axis.text.x = element_text(size = 12),
# #     axis.title.y = element_text(size = 14, face = "italic"),
# #     panel.grid = element_blank(),
# #     legend.position = "right"
# #   )+
# #   
# #   # Add bins for each category
# #   geom_vline(xintercept = 1.5, linetype = "solid", color = "black") +
# #   geom_vline(xintercept = 2.5, linetype = "solid", color = "black") +
# #   geom_vline(xintercept = 3.5, linetype = "solid", color = "black") +
# #   geom_vline(xintercept = 4.5, linetype = "solid", color = "black")
# # 
# # 
# # summary_tissue_data$Sex.Bias.Category <- factor(summary_tissue_data$Sex.Bias.Category, 
# #                                          levels = c("Highly female-biased", "Moderately female-biased", 
# #                                                      "Unbiased", "Moderately male-biased", "Highly male-biased"))
# # 
# # tissue_colors <- c("Heads" = "purple", "Heads" = "black", "Whole Body" = "gray")
# # 
# # # Sex Bias/DoS -Tissue type 
# # 
# # summary_tissue_data$Sex.Bias.Category <- factor(summary_tissue_data$Sex.Bias.Category, 
# #                                                 levels = c("Highly female-biased", "Moderately female-biased", 
# #                                                            "Unbiased", "Moderately male-biased", "Highly male-biased"))
# # 
# # tissue_colors <- c("Gonads" = "purple", "Heads" = "black", "Whole Body" = "gray")
# # 
# # 
# # dos_tissue_plot<- ggplot(summary_tissue_data, aes(x = Sex.Bias.Category, y = Mean_DoS, color = Tissue)) +
# #   geom_point(position = position_dodge(width = 0.6),
# #             size = 3) +  # Plot mean DoS without dodge
# #   geom_errorbar(aes(ymin = Mean_DoS - SE_DoS, ymax = Mean_DoS + SE_DoS), width = 0.2, position = position_dodge(width = 0.6)) +  # Error bars
# #   scale_color_manual(values = tissue_colors) +  # Colors for tissues
# #   theme_bw() +
# #   labs(title = "", x = "", y = expression(italic(DoS)), color = "Tissue Type") +
# #   geom_abline(intercept = 0, slope = 0,  size = 0.5, linetype= "solid", color = "black") +
# #   theme_classic() +
# #   theme(plot.title.position = c("panel"),
# #         legend.title = element_blank(),
# #         legend.position = c("bottom"),
# #         legend.box.background = element_rect(),
# #         legend.text = element_text(size = 20, color = "black"),
# #         axis.text.x = element_text(size=20, margin = margin(5,0,0,0), color = "black"),
# #         axis.text.y = element_text(size=20, margin = margin(0,5,0,0), color = "black"), 
# #         axis.title.x = element_text(size=40, margin = margin(10,0,0,0), color = "black"),
# #         axis.title.y = element_text(size=40, margin = margin(0,10,0,0), color = "black"),
# #         plot.title = element_text(size=40, margin = margin(0,0,0,0), color = "black"),
# #         plot.margin = margin(6,6,6,6)
# #   )
# #   
# # print(dos_tissue_plot)
# # 
# # # Expression Level Whole body
# # 
# # SBGE_expression<- SBG_whole_bodies%>%
# #   inner_join(Fra_master_data, by = "FlyBaseID")%>%  
# #   rename("exp_SBGE_ase" = "log2FoldChange")%>%
# #   select(FlyBaseID,Sex.Bias.Category,exp_SBGE_ase,exprLevelWB)%>% 
# #   distinct() # Remove disease variable to avoid repeats
# # 
# # # Creating Bins
# # bins_expr <-SBGE_expression %>%
# #   mutate(Sex.Bias.Category = case_when(
# #     exp_SBGE_ase >= 2 ~ "Highly male-biased",
# #     exp_SBGE_ase >= 0.5 & exp_SBGE_ase < 2 ~ "Moderately male-biased",
# #     exp_SBGE_ase >= -0.5 & exp_SBGE_ase < 0.5 ~ "Unbiased",
# #     exp_SBGE_ase >= -2 & exp_SBGE_ase < -0.5 ~ "Moderately female-biased",
# #     exp_SBGE_ase < -2 ~ "Highly female-biased"
# #   ))
# # 
# # 
# # summary_expr_data <- bins_expr %>%
# #   group_by(Sex.Bias.Category) %>%
# #   summarise(
# #     Mean_expr = mean(exprLevelWB, na.rm = TRUE),
# #     SE_expr = sd(exprLevelWB, na.rm = TRUE) / sqrt(n())  # Standard Error
# #   )
# # 
# # summary_expr_data$Sex.Bias.Category <- factor(summary_expr_data$Sex.Bias.Category, 
# #                                                 levels = c("Highly female-biased", "Moderately female-biased", 
# #                                                            "Unbiased", "Moderately male-biased", "Highly male-biased"))
# # # Average Expression Level Whole Body 
# # 
# # expr_plot<- ggplot(summary_expr_data, aes(x = Sex.Bias.Category, y = Mean_expr)) +
# #   geom_point(position = position_dodge(width = 0.6),
# #              size = 3) +  # Plot mean DoS without dodge
# #   geom_errorbar(aes(ymin = Mean_expr - SE_expr, ymax = Mean_expr+ SE_expr), width = 0.2, position = position_dodge(width = 0.6)) +  # Error bars
# #   theme_bw() +
# #   labs(title = "", x = "", y = "Level of Whole body expression") +
# #   geom_abline(intercept = 0, slope = 0,  size = 0.5, linetype= "solid", color = "black") +
# #   theme_classic() +
# #   theme(plot.title.position = c("panel"),
# #         legend.title = element_blank(),
# #         legend.position = c("bottom"),
# #         legend.box.background = element_rect(),
# #         legend.text = element_text(size = 20, color = "black"),
# #         axis.text.x = element_text(size=20, margin = margin(5,0,0,0), color = "black"),
# #         axis.text.y = element_text(size=20, margin = margin(0,5,0,0), color = "black"), 
# #         axis.title.x = element_text(size=40, margin = margin(10,0,0,0), color = "black"),
# #         axis.title.y = element_text(size=20, margin = margin(0,10,0,0), color = "black"),
# #         plot.title = element_text(size=40, margin = margin(0,0,0,0), color = "black"),
# #         plot.margin = margin(6,6,6,6)
# #   )
# # 
# # print(expr_plot)
# # 
# # # Expression Level Whole Body 
# # 
# # 
# # hist(SBGE_expression$exprLevelWB)
# # 
# # SBGE_expression$Sex.Bias.Category <- factor(SBGE_expression$Sex.Bias.Category, 
# #                                  levels = c("Highly female-biased", "Moderately female-biased", 
# #                                             "Unbiased", "Moderately male-biased", "Highly male-biased"))
# # 
# # colour_map <- c("Highly female-biased" = "red", "Moderately female-biased" = "pink", 
# #                "Unbiased" = "black", "Moderately male-biased" = "lightblue", "Highly male-biased" = "blue")
# # 
# # ggplot(SBGE_expression, aes(x = Sex.Bias.Category, y = exprLevelWB, color = Sex.Bias.Category)) +
# #   geom_point(size = 4, alpha = 0.7) +  # Scatter points
# #   geom_smooth(method = "lm", se = FALSE, linetype = "dashed") +  # Add regression line
# #   scale_color_manual(values = colour_map) +  # Apply custom color scheme
# #   labs(
# #     title = "Sex Bias Gene expression of whole body",
# #     x = "",
# #     y = "Expression Level (WB)",
# #     color = "Sex Bias Category"
# #   ) +
# #   theme_minimal()
# # 
# # ggplot(SBGE_expression, aes(x = Sex.Bias.Category, y = exprLevelWB, fill = Sex.Bias.Category)) +
# #   geom_boxplot(alpha = 0.7, outlier.shape = NA, width = 0.6) + # Wider box plots
# #   geom_jitter(aes(color = Sex.Bias.Category), size = 2, alpha = 0.7, width = 0.15) + # Less horizontal jitter
# #   scale_fill_manual(values = colour_map) +
# #   scale_color_manual(values = colour_map) +
# #   scale_x_discrete(expand = c(0.2, 0.2)) + # Increase spacing
# #   labs(
# #     title = "Expression Level (WB) Across Sex Bias Categories",
# #     x = "Sex Bias Category",
# #     y = "Expression Level (WB)",
# #     fill = "Sex.Bias.Category"
# #   ) +
# #   theme_minimal() +
# #   theme(axis.text.x = element_text(angle = 45, hjust = 1)) # Rotate x-axis labels
# # 
# # ggplot(SBGE_expression, aes(x = Sex.Bias.Category, y = exprLevelWB, fill = Sex.Bias.Category)) +
# #   geom_violin(alpha = 0.5, width = 0.8, draw_quantiles = c(0.25, 0.5, 0.75)) + # Adds quartiles
# #   geom_jitter(aes(color = Sex.Bias.Category), size = 1.5, width = 0.15, alpha = 0.7) + # Keeps jittered points
# #   scale_fill_manual(values = colour_map) +
# #   scale_color_manual(values = colour_map) +
# #   labs(
# #     title = "Expression Level (WB) Across Sex Bias Categories",
# #     x = "Sex Bias Category",
# #     y = "Expression Level (WB)"
# #   ) +
# #   theme_minimal() +
# #   theme(axis.text.x = element_text(angle = 45, hjust = 1))
# # 
# # ggplot(SBGE_expression, aes(x = Sex.Bias.Category, y = exprLevelWB, fill = Sex.Bias.Category)) +
# #   geom_boxplot(alpha = 0.7, outlier.shape = NA, width = 0.6) +
# #   geom_jitter(aes(color = Sex.Bias.Category), size = 2, width = 0.15, alpha = 0.7) +
# #   scale_y_log10() + # Apply log scale
# #   scale_fill_manual(values = colour_map) +
# #   scale_color_manual(values = colour_map) +
# #   labs(
# #     title = "Expression Level (WB) Across Sex Bias Categories (Log Scale)",
# #     x = "Sex Bias Category",
# #     y = "Expression Level (WB) (log scale)"
# #   ) +
# #   theme_minimal() +
# #   theme(axis.text.x = element_text(angle = 45, hjust = 1))
# # 
# # # Tajimas D 
# # 
# # # Change into numerical 
# # 
# # SBGE_tajDN<- mish_singh_SBG%>%
# #   select(FlyBaseID,Sex.Bias.Category,exp_SBGE_ase,tajD.N,Tissue)%>% 
# #   distinct() # Remove disease variable to avoid repeats
# # 
# # SBGE_tajDN<- SBGE_tajDN %>%
# #   mutate(exp_SBGE_ase = as.numeric(exp_SBGE_ase))
# # 
# # SBGE_tajDS<- mish_singh_SBG%>%
# #   select(FlyBaseID,Sex.Bias.Category,exp_SBGE_ase,tajD.S,Tissue)%>% 
# #   distinct() # Remove disease variable to avoid repeats
# # 
# # SBGE_tajDS<- SBGE_tajDS %>%
# #   mutate(exp_SBGE_ase = as.numeric(exp_SBGE_ase))
# # 
# # # Creating Bins
# # 
# # SBGE_tajDN_bins<-SBGE_tajDN %>%
# #   mutate(Sex.Bias.Category = case_when(
# #     exp_SBGE_ase >= 2 ~ "Highly male-biased",
# #     exp_SBGE_ase >= 0.5 & exp_SBGE_ase < 2 ~ "Moderately male-biased",
# #     exp_SBGE_ase >= -0.5 & exp_SBGE_ase < 0.5 ~ "Unbiased",
# #     exp_SBGE_ase >= -2 & exp_SBGE_ase < -0.5 ~ "Moderately female-biased",
# #     exp_SBGE_ase < -2 ~ "Highly female-biased"
# #   ))
# # 
# # SBGE_tajDS_bins<-SBGE_tajDS %>%
# #   mutate(Sex.Bias.Category = case_when(
# #     exp_SBGE_ase >= 2 ~ "Highly male-biased",
# #     exp_SBGE_ase >= 0.5 & exp_SBGE_ase < 2 ~ "Moderately male-biased",
# #     exp_SBGE_ase >= -0.5 & exp_SBGE_ase < 0.5 ~ "Unbiased",
# #     exp_SBGE_ase >= -2 & exp_SBGE_ase < -0.5 ~ "Moderately female-biased",
# #     exp_SBGE_ase < -2 ~ "Highly female-biased"
# #   ))
# # 
# # 
# # # Compute Mean and Standard Error (SE) for tajDN and tajDS per Sex Bias Category
# # 
# # summary_tajDN_data <-SBGE_tajDN_bins %>%
# #   group_by(Sex.Bias.Category) %>%
# #   summarise(
# #     Mean_tajD.N = mean(tajD.N, na.rm = TRUE),
# #     SE_tajD.N = sd(tajD.N, na.rm = TRUE) / sqrt(n())  # Standard Error
# #   )
# # 
# # 
# # summary_tajDS_data <-SBGE_tajDS_bins %>%
# #   group_by(Sex.Bias.Category) %>%
# #   summarise(
# #     Mean_tajD.S = mean(tajD.S, na.rm = TRUE),
# #     SE_tajD.S = sd(tajD.S, na.rm = TRUE) / sqrt(n())  # Standard Error
# #   )
# # 
# # # Merge summary tajD data
# # summary_tajD_data <- merge(summary_tajDN_data,summary_tajDS_data)
# # 
# # 
# # 
# # # Creating Summary table for tajD tissue data
# # 
# # SBGE_tajD_tissue<- mish_singh_SBG%>%
# #   select(FlyBaseID,Sex.Bias.Category,exp_SBGE_ase,tajD.S,tajD.N,Tissue)%>% 
# #   distinct() # Remove disease variable to avoid repeats
# # 
# # summary_tissue_tajDS_data <- SBGE_tajDS_bins %>%
# #   group_by(Sex.Bias.Category, Tissue) %>%
# #   summarise( 
# #     Mean_tajD.S = mean(tajD.S, na.rm = TRUE),
# #     SE_tajD.S = sd(tajD.S, na.rm = TRUE) / sqrt(n()),  # Standard Error
# #     .groups = "drop"
# #   )
# # 
# # summary_tissue_tajDN_data <- SBGE_tajDN_bins %>%
# #   group_by(Sex.Bias.Category, Tissue) %>%
# #   summarise( 
# #     Mean_tajD.N = mean(tajD.N, na.rm = TRUE),
# #     SE_tajD.N = sd(tajD.N, na.rm = TRUE) / sqrt(n()),  # Standard Error
# #     .groups = "drop"
# #   )
# # 
# # summary_tajD_tissue_data <- merge(summary_tissue_tajDN_data,summary_tissue_tajDS_data)
# # 
# # # Plot tajD non-synonymous
# # 
# # summary_tajDN_data$Sex.Bias.Category <- factor(summary_tajDN_data$Sex.Bias.Category, 
# #                                              levels = c("Highly female-biased", "Moderately female-biased", 
# #                                                         "Unbiased", "Moderately male-biased", "Highly male-biased"))
# # 
# # ggplot(summary_tajDN_data, aes(x = Sex.Bias.Category, y = Mean_tajD.N, color = Sex.Bias.Category)) +
# #   geom_point(size = 4) +  # Plot mean DoS
# #   geom_errorbar(aes(ymin = Mean_tajD.N - SE_tajD.N, ymax = Mean_tajD.N + SE_tajD.N), width = 0.2) +  # Error bars
# #   scale_color_manual(values = color_map) +  # Custom colors
# #   labs(title = "", x = "", y = "Tajimas D")+
# #   
# #   theme_bw()+
# #   theme(
# #     axis.text.x = element_text(size = 12),
# #     axis.title.y = element_text(size = 14, face = "italic"),
# #     panel.grid = element_blank(),
# #     legend.position = "right"
# #   )+
# #   
# #   # Add bins for each category
# #   geom_vline(xintercept = 1.5, linetype = "solid", color = "black") +
# #   geom_vline(xintercept = 2.5, linetype = "solid", color = "black") +
# #   geom_vline(xintercept = 3.5, linetype = "solid", color = "black") +
# #   geom_vline(xintercept = 4.5, linetype = "solid", color = "black")
# # 
# # # Plot tajD 
# # 
# # summary_tajDS_data$Sex.Bias.Category <- factor(summary_tajDS_data$Sex.Bias.Category, 
# #                                                levels = c("Highly female-biased", "Moderately female-biased", 
# #                                                           "Unbiased", "Moderately male-biased", "Highly male-biased"))
# # 
# # ggplot(summary_tajDS_data, aes(x = Sex.Bias.Category, y = Mean_tajD.S, color = Sex.Bias.Category)) +
# #   geom_point(size = 4) +  # Plot mean DoS
# #   geom_errorbar(aes(ymin = Mean_tajD.S - SE_tajD.S, ymax = Mean_tajD.S + SE_tajD.S), width = 0.2) +  # Error bars
# #   scale_color_manual(values = color_map) +  # Custom colors
# #   labs(title = "", x = "", y = "Tajimas D")+
# #   
# #   theme_bw()+
# #   theme(
# #     axis.text.x = element_text(size = 12),
# #     axis.title.y = element_text(size = 14, face = "italic"),
# #     panel.grid = element_blank(),
# #     legend.position = "right"
# #   )+
# #   
# #   # Add bins for each category
# #   geom_vline(xintercept = 1.5, linetype = "solid", color = "black") +
# #   geom_vline(xintercept = 2.5, linetype = "solid", color = "black") +
# #   geom_vline(xintercept = 3.5, linetype = "solid", color = "black") +
# #   geom_vline(xintercept = 4.5, linetype = "solid", color = "black")
# # 
# # # Plot tajD tissue
# # 
# # summary_tajD_tissue_data$Sex.Bias.Category <- factor(summary_tajD_tissue_data$Sex.Bias.Category, 
# #                                             levels = c("Highly female-biased", "Moderately female-biased", 
# #                                                        "Unbiased", "Moderately male-biased", "Highly male-biased"))
# # 
# # tajD_tissue_plot<- ggplot(summary_tajD_tissue_data, aes(x = Sex.Bias.Category, y =  Mean_tajD.S, color = Tissue)) +
# #   geom_point(position = position_dodge(width = 0.6),
# #              size = 3) +  # Plot mean DoS without dodge
# #   geom_errorbar(aes(ymin =  Mean_tajD.S - SE_tajD.S, ymax = Mean_tajD.S+ SE_tajD.S), width = 0.2, position = position_dodge(width = 0.6)) +  # Error bars
# #   scale_color_manual(values = tissue_colors) +  # Colors for tissues
# #   theme_bw() +
# #   labs(title = "", x = "", y = "Tajimas D", color = "Tissue Type") +
# #   geom_abline(intercept = 0, slope = 0,  size = 0.5, linetype= "solid", color = "black") +
# #   theme_classic() +
# #   theme(plot.title.position = c("panel"),
# #         legend.title = element_blank(),
# #         legend.position = c("right"),
# #         legend.box.background = element_rect(),
# #         legend.text = element_text(size = 20, color = "black"),
# #         axis.text.x = element_text(angle = 20, size=15, hjust = 1, margin = margin(5,0,0,0), color = "black"),
# #         axis.text.y = element_text(size=20, margin = margin(0,5,0,0), color = "black"), 
# #         axis.title.x = element_text(size=20, margin = margin(10,0,0,0), color = "black"),
# #         axis.title.y = element_text(size=20, margin = margin(0,10,0,0), color = "black"),
# #         plot.title = element_text(size=20, margin = margin(0,0,0,0), color = "black"),
# #         plot.margin = margin(6,6,6,6)
# #   )+
# #   
# #   # Add bins for each category
# #   geom_vline(xintercept = 1.5, linetype = "solid", color = "black") +
# #   geom_vline(xintercept = 2.5, linetype = "solid", color = "black") +
# #   geom_vline(xintercept = 3.5, linetype = "solid", color = "black") +
# #   geom_vline(xintercept = 4.5, linetype = "solid", color = "black")
# # 
# # print(tajD_tissue_plot)
# # 
# # 
# # # Chromosomal distribution 
# # 
# # # Convert Chr to factor with an order: X, 2, then 3
# # chr_counts$Chr <- factor(chr_counts$Chr, levels = c("X", "2", "3"))
# # 
# # # Create the bar plot 
# # ggplot(chr_counts, aes(x = Chr, y = count, fill = Chr)) +
# #   geom_bar(stat = "identity") +
# #   labs(title = "Distribution of Chromosomes in Disease Genes",
# #        x = "Chromosome",
# #        y = "Number of Disease-Associated Genes") +
# #   theme_minimal() + 
# #   theme(panel.grid.major = element_blank(),  # Remove major gridlines
# #         panel.grid.minor = element_blank())  # Remove minor gridlines
# # 
# # # rMF 
# # 
# # # rMF data 
# # 
# # SBGE_rMF<- mish_singh_SBG%>%
# #   select(FlyBaseID,Sex.Bias.Category,exp_SBGE_ase,DoS,Tissue)%>% 
# #   distinct() # Remove disease variable to avoid repeats
# # 
# # # Change into numerical 
# # str(SBGE_rMF$exp_SBGE_ase)
# # SBGE_rMF<- SBGE_rMF %>%
# #   mutate(exp_SBGE_ase = as.numeric(exp_SBGE_ase))
# # 
# # # Creating Bins
# # SBGE_rMF_bins<-SBGE_rMF %>%
# #   mutate(Sex.Bias.Category = case_when(
# #     exp_SBGE_ase >= 2 ~ "Highly male-biased",
# #     exp_SBGE_ase >= 0.5 & exp_SBGE_ase < 2 ~ "Moderately male-biased",
# #     exp_SBGE_ase >= -0.5 & exp_SBGE_ase < 0.5 ~ "Unbiased",
# #     exp_SBGE_ase >= -2 & exp_SBGE_ase < -0.5 ~ "Moderately female-biased",
# #     exp_SBGE_ase < -2 ~ "Highly female-biased"
# #   ))
# # 
# # # Compute Mean and Standard Error (SE) for rMF per Sex Bias Category
# # 
# # summary_rMF_data <-SBGE_rMF_bins %>%
# #   group_by(Sex.Bias.Category) %>%
# #   summarise(
# #     Mean_rmf = mean(rmf, na.rm = TRUE),
# #     SE_rmf = sd(rmf, na.rm = TRUE) / sqrt(n())  # Standard Error
# #   )
# # 
# # 
# # # Creating Summary table for rMF tissue data
# # 
# # SBGE_rMF_tissue<- mish_singh_SBG%>%
# #   select(FlyBaseID,Sex.Bias.Category,exp_SBGE_ase,rmf,Tissue)%>% 
# #   distinct() # Remove disease variable to avoid repeats
# # 
# # summary_tissue_rMF_data <- SBGE_rMF_bins %>%
# #   group_by(Sex.Bias.Category, Tissue) %>%
# #   summarise(
# #     Mean_rmf = mean(rmf, na.rm = TRUE),
# #     SE_rmf = sd(rmf, na.rm = TRUE) / sqrt(n()),  # Standard Error
# #     .groups = "drop"
# #   )
# # 
# # 
# # # rMF no tissue plot
# # 
# # summary_rMF_data$Sex.Bias.Category <- factor(summary_rMF_data$Sex.Bias.Category, 
# #                                                     levels = c("Highly female-biased", "Moderately female-biased", 
# #                                                                "Unbiased", "Moderately male-biased", "Highly male-biased"))
# # 
# # ggplot(summary_rMF_data, aes(x = Sex.Bias.Category, y = Mean_rmf, color = Sex.Bias.Category)) +
# #   geom_point(size = 4) +  # Plot mean DoS
# #   geom_errorbar(aes(ymin = Mean_rmf - SE_rmf, ymax = Mean_rmf + SE_rmf), width = 0.2) +  # Error bars
# #   scale_color_manual(values = color_map) +  # Custom colors
# #   labs(title = "", x = "", y = "rmf")+
# #   
# #   theme_bw()+
# #   theme(
# #     axis.text.x = element_text(size = 12),
# #     axis.title.y = element_text(size = 14, face = "italic"),
# #     panel.grid = element_blank(),
# #     legend.position = "right"
# #   )+
# #   
# #   # Add bins for each category
# #   geom_vline(xintercept = 1.5, linetype = "solid", color = "black") +
# #   geom_vline(xintercept = 2.5, linetype = "solid", color = "black") +
# #   geom_vline(xintercept = 3.5, linetype = "solid", color = "black") +
# #   geom_vline(xintercept = 4.5, linetype = "solid", color = "black")
# # 
# # summary_tissue_rMF_data$Sex.Bias.Category <- factor(summary_tissue_rMF_data$Sex.Bias.Category, 
# #                                                 levels = c("Highly female-biased", "Moderately female-biased", 
# #                                                            "Unbiased", "Moderately male-biased", "Highly male-biased"))
# # 
# # # rMF tissue plot
# # 
# # rMF_tissue_plot<- ggplot(summary_tissue_rMF_data, aes(x = Sex.Bias.Category, y =  Mean_rmf, color = Tissue)) +
# #   geom_point(position = position_dodge(width = 0.6),
# #              size = 3) +  # Plot mean DoS without dodge
# #   geom_errorbar(aes(ymin =  Mean_rmf - SE_rmf, ymax = Mean_rmf + SE_rmf), width = 0.2, position = position_dodge(width = 0.6)) +  # Error bars
# #   scale_color_manual(values = tissue_colors) +  # Colors for tissues
# #   theme_bw() +
# #   labs(title = "", x = "", y = "rMF", color = "Tissue Type") +
# #   geom_abline(intercept = 0, slope = 0,  size = 0.5, linetype= "solid", color = "black") +
# #   theme_classic() +
# #   theme(plot.title.position = c("panel"),
# #         legend.title = element_blank(),
# #         legend.position = c("right"),
# #         legend.box.background = element_rect(),
# #         legend.text = element_text(size = 20, color = "black"),
# #         axis.text.x = element_text(angle = 45, size=15, hjust = 1, margin = margin(5,0,0,0), color = "black"),
# #         axis.text.y = element_text(size=20, margin = margin(0,5,0,0), color = "black"), 
# #         axis.title.x = element_text(size=40, margin = margin(10,0,0,0), color = "black"),
# #         axis.title.y = element_text(size=40, margin = margin(0,10,0,0), color = "black"),
# #         plot.title = element_text(size=40, margin = margin(0,0,0,0), color = "black"),
# #         plot.margin = margin(6,6,6,6)
# #   )+
# #   
# #   # Add bins for each category
# #   geom_vline(xintercept = 1.5, linetype = "solid", color = "black") +
# #   geom_vline(xintercept = 2.5, linetype = "solid", color = "black") +
# #   geom_vline(xintercept = 3.5, linetype = "solid", color = "black") +
# #   geom_vline(xintercept = 4.5, linetype = "solid", color = "black")
# # 
# # print(rMF_tissue_plot)
# # 
# # #COMPARATIVE PLOTS WITH CONTROL GENEs
# 
# # Tajimas D 
# 
# # Change into numerical 
# 
# SBGE_tajDS_compare<- control_set_genes_mish_singh_SBGE_wb%>%
#   select(FlyBaseID,Sex.Bias.Category,exp_SBGE_ase,tajD.S,GeneType)%>% 
#   distinct() # Remove disease variable to avoid repeats
# 
# SBGE_tajDS_compare<- SBGE_tajDS_compare %>%
#   mutate(exp_SBGE_ase = as.numeric(exp_SBGE_ase))
# 
# # Creating Bins
# 
# SBGE_tajDS_bins_compare<-SBGE_tajDS_compare %>%
#   mutate(Sex.Bias.Category = case_when(
#     exp_SBGE_ase >= 2 ~ "Highly male-biased",
#     exp_SBGE_ase >= 0.5 & exp_SBGE_ase < 2 ~ "Moderately male-biased",
#     exp_SBGE_ase >= -0.5 & exp_SBGE_ase < 0.5 ~ "Unbiased",
#     exp_SBGE_ase >= -2 & exp_SBGE_ase < -0.5 ~ "Moderately female-biased",
#     exp_SBGE_ase < -2 ~ "Highly female-biased"
#   ))
# 
# SBGE_tajDS_bins_compare <- SBGE_tajDS_bins_compare[!is.na(SBGE_tajDS_bins_compare$tajD.S), ]
# 
# # Compute Mean and Standard Error (SE) for tajDN and tajDS per Sex Bias Category
# 
# # Calculate mean
# mean_tajD <- aggregate(tajD.S ~ Sex.Bias.Category + GeneType, 
#                        data = SBGE_tajDS_bins_compare, 
#                        FUN = function(x) mean(x, na.rm = TRUE))
# 
# # Calculate standard error
# se_tajD <- aggregate(tajD.S ~ Sex.Bias.Category + GeneType, 
#                      data = SBGE_tajDS_bins_compare, 
#                      FUN = function(x) sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x))))
# 
# # Combine into one data frame
# summary_tajDS_data_comparative <- merge(mean_tajD, se_tajD, 
#                                         by = c("Sex.Bias.Category", "GeneType"))
# 
# # Rename columns
# colnames(summary_tajDS_data_comparative) <- c("Sex.Bias.Category", "GeneType", 
#                                               "Mean_tajD.S", "SE_tajD.S")
# 
# 
# # summary_tissue_tajDS_data_comparative <-SBGE_tajDS_bins_compare%>%
# #   group_by(Sex.Bias.Category, Tissue, GeneType) %>%
# #   summarise( 
# #     Mean_tajD.S = mean(tajD.S, na.rm = TRUE),
# #     SE_tajD.S = sd(tajD.S, na.rm = TRUE) / sqrt(n()),  # Standard Error
# #     .groups = "drop"
# #   )
# 
# 
# # # Plot using a single x-axis with grouped categories
# # ggplot(summary_tajDS_data_comparative, aes(x = Sex.Bias.Category, y = Mean_tajD.S, color = GeneType)) +
# #   geom_pointrange(aes(ymin = Mean_tajD.S - SE_tajD.S, ymax = Mean_tajD.S + SE_tajD.S),
# #                   position = position_dodge(width = 0.5),
# #                   size = 1) +  # Adjust error bar size
# #   geom_point(position = position_dodge(width = 0.5), size = 3) +  # Set size of points
# #   scale_color_manual(values = c("Control" = "blue", "Disease" = "red")) +  # Color points by Gene Type
# #   labs(x = "Sex Bias Category", y = "Tajima's D") +
# #   theme_minimal(base_size = 14) +  # Remove default gridlines
# #   theme(
# #     panel.grid.major = element_blank(),  # Remove major gridlines
# #     panel.grid.minor = element_blank(),  # Remove minor gridlines
# #     panel.border = element_rect(color = "black", fill = NA, size = 1),  # Add border around plot
# #     axis.text.x = element_text(angle = 45, hjust = 1),
# #     panel.background = element_blank()
# #   ) +
# #   geom_vline(xintercept = seq(0.5, length(unique(summary_data$Sex.Bias.Category)) + 0.5, 1),
# #              linetype = "solid", color = "black", size = 0.8)  # Add bin separators
# 
# 
# #summary_tissue_tajDS_data_comparative$Sex.Bias.Category <- factor(summary_tissue_tajDS_data_comparative$Sex.Bias.Category,
#                                                      #levels = c("Highly female-biased", "Moderately female-biased",
#                                                                 #"Unbiased", "Moderately male-biased", "Highly male-biased"))
# 
# # 
# # ggplot(summary_tissue_tajDS_data_comparative, aes(x = Sex.Bias.Category, y = Mean_tajD.S, color = GeneType)) +
# #   geom_pointrange(aes(ymin = Mean_tajD.S - SE_tajD.S, ymax = Mean_tajD.S + SE_tajD.S),
# #                   position = position_dodge(width = 0.5),
# #                   size = 1) +  # Adjust error bar size
# #   geom_point(position = position_dodge(width = 0.5), size = 3) +  # Set size of points
# #   scale_color_manual(values = c("Control" = "blue", "Disease" = "red")) +  # Color points by Gene Type
# #   labs(x = "Sex Bias Category", y = "Tajima's D") +
# #   theme_minimal(base_size = 14) +  # Remove default gridlines
# #   theme(
# #     panel.grid.major = element_blank(),  # Remove major gridlines
# #     panel.grid.minor = element_blank(),  # Remove minor gridlines
# #     panel.border = element_rect(color = "black", fill = NA, size = 1),  # Add border around plot
# #     axis.text.x = element_text(angle = 45, hjust = 1),
# #     panel.background = element_blank()
# #   ) +
# #   geom_vline(xintercept = seq(0.5, length(unique(summary_data$Sex.Bias.Category)) + 0.5, 1),
# #              linetype = "solid", color = "black", size = 0.8) +  # Add bin separators
# #   facet_wrap(~Tissue, ncol = 1, )  # Facet by Tissue, one column
# # 
# 
# sample_sizes_tadj <- as.data.frame(table(SBGE_tajDS_bins_compare$Sex.Bias.Category,
#                                          SBGE_tajDS_bins_compare$GeneType))
# 
# # Rename columns
# colnames(sample_sizes_tadj) <- c("Sex.Bias.Category", "GeneType", "SampleSize")
# 
# summary_tajDS_data_comparative <- left_join(summary_tajDS_data_comparative, sample_sizes_tadj, by = c("Sex.Bias.Category", "GeneType"))
# 
# # Then compute the "All" category stats
# 
# split_data<- split(summary_tajDS_data_comparative, summary_tajDS_data_comparative$GeneType)
# 
# all_category <- do.call(rbind, lapply(split_data, function(df) {
#   data.frame(
#     GeneType = unique(df$GeneType),
#     Mean_tajD.S = mean(df$Mean_tajD.S, na.rm = TRUE),
#     SE_tajD.S = sqrt(sum(df$SE_tajD.S^2, na.rm = TRUE) / nrow(df)),
#     SampleSize = sum(df$SampleSize, na.rm = TRUE),
#     Sex.Bias.Category = "All"
#   )
# }))
# 
# summary_tajDS_data_comparative <- bind_rows(summary_tajDS_data_comparative, all_category)
# 
# summary_tajDS_data_comparative <- left_join(summary_tajDS_data_comparative, tajD_boot_results, by = c("Sex.Bias.Category", "GeneType"))
# 
# summary_tajDS_data_comparative$Sex.Bias.Category <- factor(summary_tajDS_data_comparative$Sex.Bias.Category, 
#                                                      levels = c("All","Highly female-biased", "Moderately female-biased", 
#                                                                 "Unbiased", "Moderately male-biased", "Highly male-biased"))
# 
# # Plot using a single x-axis with grouped categories
# 
# ggplot(summary_tajDS_data_comparative, 
#        aes(x = Sex.Bias.Category, y = q50, color = GeneType)) + 
#   geom_pointrange(aes(ymin = q05, ymax = q95),
#                   position = position_dodge(width = 0.5),
#                   size = 1) +  # Adjust error bar thickness
#   geom_point(position = position_dodge(width = 0.5), size = 3) +  # Set size of points
#   geom_text(aes(label = SampleSize),  
#             position = position_dodge(width = 0.5),  
#             color = "black", size = 3) + # Adjust label position and size
#   scale_color_manual(values = c("Control" = "darkgreen", "Disease" = "purple")) +  # Match colors
#   labs(x = "Sex Bias Category", y = "Tajima’s D",color = "Gene Type") +
#   guides(color = guide_legend(override.aes = list(size = 1))) +
#   theme_minimal(base_size = 14) +  # Remove default gridlines
#   theme(
#     panel.grid.major = element_blank(),  # Remove major gridlines
#     panel.grid.minor = element_blank(),  # Remove minor gridlines
#     panel.border = element_rect(color = "black", fill = NA, size = 1),  # Add border
#     panel.background = element_blank()
#   ) +
#   geom_vline(xintercept = 1.5, color = "black", size = 1.2) +  # Bold divider
#   geom_vline(xintercept = 2.5:5.5, color = "grey60", size = 0.5)+  # Thin grey lines
#   scale_x_discrete(labels = c("All", "Highly FB", "Female-Biased", "Unbiased", "Male-Biased", "Highly MB")) 
# 
# 
# 
# # rMF
# 
# # Change into numerical 
# 
# SBGE_rMF_compare<- control_set_genes_mish_singh_SBGE_wb%>%
#   select(FlyBaseID,Sex.Bias.Category,exp_SBGE_ase,Tissue,rmf,GeneType)%>% 
#   distinct() # Remove disease variable to avoid repeats
# 
# SBGE_rMF_compare<- SBGE_rMF_compare%>%
#   mutate(exp_SBGE_ase = as.numeric(exp_SBGE_ase))
# 
# # Creating Bins
# 
# SBGE_rMF_bins_compare<-SBGE_rMF_compare%>%
#   mutate(Sex.Bias.Category = case_when(
#     exp_SBGE_ase >= 2 ~ "Highly male-biased",
#     exp_SBGE_ase >= 0.5 & exp_SBGE_ase < 2 ~ "Moderately male-biased",
#     exp_SBGE_ase >= -0.5 & exp_SBGE_ase < 0.5 ~ "Unbiased",
#     exp_SBGE_ase >= -2 & exp_SBGE_ase < -0.5 ~ "Moderately female-biased",
#     exp_SBGE_ase < -2 ~ "Highly female-biased"
#   ))
# 
# SBGE_rMF_bins_compare <- SBGE_rMF_bins_compare[!is.na(SBGE_rMF_bins_compare$rmf), ]
# 
# # Compute Mean and Standard Error (SE) for tajDN and tajDS per Sex Bias Category
# 
# # Compute mean
# mean_rmf <- aggregate(rmf ~ Sex.Bias.Category + GeneType, 
#                       data = SBGE_rMF_bins_compare, 
#                       FUN = function(x) mean(x, na.rm = TRUE))
# 
# # Compute standard error
# se_rmf <- aggregate(rmf ~ Sex.Bias.Category + GeneType, 
#                     data = SBGE_rMF_bins_compare, 
#                     FUN = function(x) sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x))))
# 
# # Merge the two summaries
# summary_rMF_data_comparative <- merge(mean_rmf, se_rmf, 
#                                       by = c("Sex.Bias.Category", "GeneType"))
# 
# # Rename columns
# colnames(summary_rMF_data_comparative) <- c("Sex.Bias.Category", "GeneType", 
#                                             "Mean_rmf", "SE_rmf")
# 
# 
# # summary_tissue_rMF_data_comparative <-SBGE_rMF_bins_compare%>%
# #   group_by(Sex.Bias.Category, Tissue, GeneType) %>%
# #   summarise( 
# #     Mean_rmf = mean(rmf, na.rm = TRUE),
# #     SE_rmf = sd(rmf, na.rm = TRUE) / sqrt(n()),  # Standard Error
# #     .groups = "drop"
# #   )
# 
# rMF_sample_sizes <- as.data.frame(table(SBGE_rMF_bins_compare$Sex.Bias.Category,
#                                         SBGE_rMF_bins_compare$GeneType))
# 
# # Rename columns 
# colnames(rMF_sample_sizes) <- c("Sex.Bias.Category", "GeneType", "SampleSize")
# summary_rMF_data_comparative <- left_join(summary_rMF_data_comparative, rMF_sample_sizes, by = c("Sex.Bias.Category", "GeneType"))
# 
# # Split by GeneType
# split_data_rmf <- split(summary_rMF_data_comparative, summary_rMF_data_comparative$GeneType)
# 
# # For each GeneType, compute "All" summary row
# all_category_rmf <- do.call(rbind, lapply(split_data_rmf , function(df) {
#   data.frame(
#     GeneType = unique(df$GeneType),
#     Mean_rmf = mean(df$Mean_rmf, na.rm = TRUE),
#     SE_rmf = sqrt(sum(df$SE_rmf^2, na.rm = TRUE) / nrow(df)),
#     SampleSize = sum(df$SampleSize, na.rm = TRUE),
#     Sex.Bias.Category = "All"
#   )
# }))
# 
# summary_rMF_data_comparative <- bind_rows(summary_rMF_data_comparative,all_category_rmf)
# 
# summary_rMF_data_comparative <- left_join(summary_rMF_data_comparative, rmf_boot_results, by = c("Sex.Bias.Category", "GeneType"))
# 
# summary_rMF_data_comparative$Sex.Bias.Category <- factor(summary_rMF_data_comparative$Sex.Bias.Category, 
#                                                            levels = c("All","Highly female-biased", "Moderately female-biased", 
#                                                                       "Unbiased", "Moderately male-biased", "Highly male-biased"))
# # rmf Plot
# 
# ggplot(summary_rMF_data_comparative, 
#        aes(x = Sex.Bias.Category, y = q50, color = GeneType)) + 
#   geom_pointrange(aes(ymin = q05, ymax = q95),
#                   position = position_dodge(width = 0.5),
#                   size = 1) +  # Adjust error bar thickness
#   geom_point(position = position_dodge(width = 0.5), size = 2) +  # Set size of points
#   geom_text(aes(label = SampleSize),  
#             position = position_dodge(width = 0.5),
#             color = "black", size = 3) + # Adjust label position and size
#   scale_color_manual(values = c("Control" = "darkgreen", "Disease" = "purple")) +  # Match colors
#   labs(x = "Sex Bias Category", y = "rMF",color = "Gene Type") +
#   guides(color = guide_legend(override.aes = list(size = 1))) +
#   theme_minimal(base_size = 14) +  # Remove default gridlines
#   theme(
#     panel.grid.major = element_blank(),  # Remove major gridlines
#     panel.grid.minor = element_blank(),  # Remove minor gridlines
#     panel.border = element_rect(color = "black", fill = NA, size = 1),  # Add border
#     panel.background = element_blank()
#   )+
#   geom_vline(xintercept = 1.5, color = "black", size = 1.2) +  # Bold divider
#   geom_vline(xintercept = 2.5:5.5, color = "grey60", size = 0.5)+  # Thin grey lines
#   scale_x_discrete(labels = c("All", "Highly FB", "Female-Biased", "Unbiased", "Male-Biased", "Highly MB")) 
# 
# 
# 
# # DoS 
# # Change into numerical 
# 
# SBGE_DoS_compare<- control_set_genes_mish_Fra_SBGE_wb%>%
#   select(FlyBaseID,Sex.Bias.Category,exp_SBGE_ase,DoS,GeneType,disease)%>% 
#   distinct() # Remove disease variable to avoid repeats
# 
# SBGE_DoS_compare<- SBGE_DoS_compare%>%
#   mutate(exp_SBGE_ase = as.numeric(exp_SBGE_ase))
# 
# # Creating Bins
# 
# SBGE_DoS_bins_compare<-SBGE_DoS_compare%>%
#   mutate(Sex.Bias.Category = case_when(
#     exp_SBGE_ase >= 2 ~ "Highly male-biased",
#     exp_SBGE_ase >= 0.5 & exp_SBGE_ase < 2 ~ "Moderately male-biased",
#     exp_SBGE_ase >= -0.5 & exp_SBGE_ase < 0.5 ~ "Unbiased",
#     exp_SBGE_ase >= -2 & exp_SBGE_ase < -0.5 ~ "Moderately female-biased",
#     exp_SBGE_ase < -2 ~ "Highly female-biased"
#   ))
# 
# # Compute Mean and Standard Error (SE) for tajDN and tajDS per Sex Bias Category
# 
# # Calculate mean
# mean_DoS <- aggregate(DoS ~ Sex.Bias.Category + GeneType, 
#                       data = SBGE_DoS_bins_compare, 
#                       FUN = function(x) mean(x, na.rm = TRUE))
# 
# # Calculate standard error
# se_DoS <- aggregate(DoS ~ Sex.Bias.Category + GeneType, 
#                     data = SBGE_DoS_bins_compare, 
#                     FUN = function(x) sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x))))
# 
# # Merge both summaries
# summary_DoS_data_comparative <- merge(mean_DoS, se_DoS,
#                                       by = c("Sex.Bias.Category", "GeneType"))
# 
# # Rename columns
# colnames(summary_DoS_data_comparative) <- c("Sex.Bias.Category", "GeneType", "Mean_DoS", "SE_DoS")
# 
# 
# DoS_sample_sizes <- as.data.frame(table(SBGE_DoS_bins_compare$Sex.Bias.Category,
#                                         SBGE_DoS_bins_compare$GeneType))
# 
# # Rename columns
# colnames(DoS_sample_sizes) <- c("Sex.Bias.Category", "GeneType", "SampleSize")
# 
# summary_DoS_data_comparative <- left_join(summary_DoS_data_comparative, DoS_sample_sizes, by = c("Sex.Bias.Category", "GeneType"))
# 
# 
# # Split by GeneType
# split_data_DoS <- split(summary_DoS_data_comparative, summary_DoS_data_comparative$GeneType)
# 
# # For each GeneType, compute "All" summary row
# all_category_DoS <- do.call(rbind, lapply(split_data_DoS , function(df) {
#   data.frame(
#     GeneType = unique(df$GeneType),
#     Mean_DoS = mean(df$Mean_DoS, na.rm = TRUE),
#     SE_DoS = sqrt(sum(df$SE_DoS^2, na.rm = TRUE) / nrow(df)),
#     SampleSize = sum(df$SampleSize, na.rm = TRUE),
#     Sex.Bias.Category = "All"
#   )
# }))
# 
# # Reorder columns to match the expected format
# DoS_all_category <- DoS_all_category[, c("Sex.Bias.Category", "GeneType", "Mean_DoS", "SE_DoS", "SampleSize")]
# 
# 
# summary_DoS_data_comparative <- bind_rows(summary_DoS_data_comparative,all_category_DoS)
# 
# summary_DoS_data_comparative <- left_join(summary_DoS_data_comparative, DoS_boot_results, by = c("Sex.Bias.Category", "GeneType"))
# 
# summary_DoS_data_comparative$Sex.Bias.Category <- factor(summary_DoS_data_comparative$Sex.Bias.Category, 
#                                                          levels = c("All","Highly female-biased", "Moderately female-biased", 
#                                                                     "Unbiased", "Moderately male-biased", "Highly male-biased"))
# # DoS Plot
# 
# ggplot(summary_DoS_data_comparative, 
#        aes(x = Sex.Bias.Category, y = q50, color = GeneType)) + 
#   geom_pointrange(aes(ymin = q05, ymax = q95),
#                   position = position_dodge(width = 0.5),
#                   size = 1) +  # Adjust error bar thickness
#   geom_point(position = position_dodge(width = 0.5), size = 2) +  # Set size of points
#   geom_text(aes(label = SampleSize),  
#             position = position_dodge(width = 0.5),
#             color = "black", size = 3) + # Adjust label position and size
#   scale_color_manual(values = c("Control" = "darkgreen", "Disease" = "purple")) +  # Match colors
#   labs(x = "Sex Bias Category", y = "DoS",color = "Gene Type") +
#   guides(color = guide_legend(override.aes = list(size = 1))) +
#   theme_minimal(base_size = 14) +  # Remove default gridlines
#   theme(
#     panel.grid.major = element_blank(),  # Remove major gridlines
#     panel.grid.minor = element_blank(),  # Remove minor gridlines
#     panel.border = element_rect(color = "black", fill = NA, size = 1),  # Add border
#     panel.background = element_blank()
#   )+
#   geom_vline(xintercept = 1.5, color = "black", size = 1.2) +  # Bold divider
#   geom_vline(xintercept = 2.5:5.5, color = "grey60", size = 0.5)+  # Thin grey lines
#   scale_x_discrete(labels = c("All", "Highly FB", "Female-Biased", "Unbiased", "Male-Biased", "Highly MB")) 
# 
# # Tissue Specificity
# 
# # Change into numerical 
# 
# SBGE_tissuespec_compare<- control_set_genes_mish_Fra_SBGE_wb%>%
#   select(FlyBaseID,Sex.Bias.Category,exp_SBGE_ase,tissueSpecificityAdults,GeneType,disease)%>% 
#   distinct() # Remove disease variable to avoid repeats
# 
# SBGE_tissuespec_compare<- SBGE_tissuespec_compare%>%
#   mutate(exp_SBGE_ase = as.numeric(exp_SBGE_ase))
# 
# # Creating Bins
# 
# SBGE_tissuespec_bins_compare<-SBGE_tissuespec_compare%>%
#   mutate(Sex.Bias.Category = case_when(
#     exp_SBGE_ase >= 2 ~ "Highly male-biased",
#     exp_SBGE_ase >= 0.5 & exp_SBGE_ase < 2 ~ "Moderately male-biased",
#     exp_SBGE_ase >= -0.5 & exp_SBGE_ase < 0.5 ~ "Unbiased",
#     exp_SBGE_ase >= -2 & exp_SBGE_ase < -0.5 ~ "Moderately female-biased",
#     exp_SBGE_ase < -2 ~ "Highly female-biased"
#   ))
# 
# SBGE_tissuespec_bins_compare <- SBGE_tissuespec_bins_compare[!is.na(SBGE_tissuespec_bins_compare$tissueSpecificityAdults), ]
# 
# # Compute Mean and Standard Error (SE) for tajDN and tajDS per Sex Bias Category
# 
# # Calculate mean
# mean_tissueSpecificityAdults<- aggregate(tissueSpecificityAdults~ Sex.Bias.Category + GeneType, 
#                       data = SBGE_tissuespec_bins_compare, 
#                       FUN = function(x) mean(x, na.rm = TRUE))
# 
# # Calculate standard error
# se_tissueSpecificityAdults<- aggregate(tissueSpecificityAdults ~ Sex.Bias.Category + GeneType, 
#                     data = SBGE_tissuespec_bins_compare, 
#                     FUN = function(x) sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x))))
# 
# # Merge both summaries
# summary_tissuespec_data_comparative <- merge(mean_tissueSpecificityAdults,se_tissueSpecificityAdults,
#                                       by = c("Sex.Bias.Category", "GeneType"))
# 
# # Rename columns
# colnames(summary_tissuespec_data_comparative) <- c("Sex.Bias.Category", "GeneType", "Mean_tissueSpecificityAdults", "SE_tissueSpecificityAdults")
# 
# 
# tissuespec_sample_sizes <- as.data.frame(table(SBGE_tissuespec_bins_compare$Sex.Bias.Category,
#                                         SBGE_tissuespec_bins_compare$GeneType))
# 
# # Rename columns
# colnames(tissuespec_sample_sizes) <- c("Sex.Bias.Category", "GeneType", "SampleSize")
# 
# summary_tissuespec_data_comparative <- left_join(summary_tissuespec_data_comparative, tissuespec_sample_sizes, by = c("Sex.Bias.Category", "GeneType"))
# 
# 
# # Split by GeneType
# split_data_tissuespec <- split(summary_tissuespec_data_comparative, summary_tissuespec_data_comparative$GeneType)
# 
# # For each GeneType, compute "All" summary row
# all_category_tissuespec <- do.call(rbind, lapply(split_data_tissuespec , function(df) {
#   data.frame(
#     GeneType = unique(df$GeneType),
#     Mean_tissueSpecificityAdults = mean(df$Mean_tissueSpecificityAdults, na.rm = TRUE),
#     SE_tissueSpecificityAdults = sqrt(sum(df$SE_tissueSpecificityAdults^2, na.rm = TRUE) / nrow(df)),
#     SampleSize = sum(df$SampleSize, na.rm = TRUE),
#     Sex.Bias.Category = "All"
#   )
# }))
# 
# # Reorder columns to match the expected format
# all_category_tissuespec <- all_category_tissuespec[, c("Sex.Bias.Category", "GeneType", "Mean_tissueSpecificityAdults", "SE_tissueSpecificityAdults", "SampleSize")]
# 
# 
# summary_tissuespec_data_comparative <- bind_rows(summary_tissuespec_data_comparative,all_category_tissuespec)
# 
# summary_tissuespec_data_comparative<- left_join(summary_tissuespec_data_comparative, TissueSpecificity_boot_results, by = c("Sex.Bias.Category", "GeneType"))
# 
# summary_tissuespec_data_comparative$Sex.Bias.Category <- factor(summary_tissuespec_data_comparative$Sex.Bias.Category, 
#                                                          levels = c("All","Highly female-biased", "Moderately female-biased", 
#                                                                     "Unbiased", "Moderately male-biased", "Highly male-biased"))
# # Tissue Specificity Plot
# 
# tissueSpecificity_plot<-ggplot(summary_tissuespec_data_comparative, 
#        aes(x = Sex.Bias.Category, y = q50, color = GeneType)) + 
#   geom_pointrange(aes(ymin = q05, ymax = q95),
#                   position = position_dodge(width = 0.5),
#                   size = 0.7) +  # Adjust error bar thickness
#   geom_point(position = position_dodge(width = 0.5), size = 2.5) +  # Set size of points
#   geom_text(aes(label = SampleSize),  
#             position = position_dodge(width = 0.5),
#             color = "black", size = 4) + # Adjust label position and size
#   scale_color_manual(values = c("Control" = "darkgreen", "Disease" = "purple")) +  # Match colors
#   labs(x = "Sex Bias Category", y = "Tissue Specificity",color = "Gene Type") +
#   guides(color = guide_legend(override.aes = list(size = 1))) +
#   theme_minimal(base_size = 16) +  # Remove default gridlines
#   theme(axis.text = element_text(size = 14),
#         axis.title = element_text(size = 16, face = "bold"),
#         legend.text = element_text(size = 14),
#         legend.title = element_text(size = 16),
#         strip.text = element_text(size = 16),
#         plot.margin = margin(10, 10, 10, 10),
#     panel.grid.major = element_blank(),  # Remove major gridlines
#     panel.grid.minor = element_blank(),  # Remove minor gridlines
#     panel.border = element_rect(color = "black", fill = NA, size = 1),  # Add border
#     panel.background = element_blank()
#   )+
#   geom_vline(xintercept = 1.5, color = "black", size = 1.2) +  # Bold divider
#   geom_vline(xintercept = 2.5:5.5, color = "grey60", size = 0.5)+  # Thin grey lines
#   scale_x_discrete(labels = c("All", "Highly FB", "Female-Biased", "Unbiased", "Male-Biased", "Highly MB")) 
# 
# 
# ggsave("tissueSpecificity_plot.pdf", width = 12, height = 7, dpi = 300)
# 
# tissueSpecificity_plot2 <- ggplot(summary_tissuespec_data_comparative, 
#                                  aes(x = Sex.Bias.Category, y = q50, color = GeneType)) + 
#   geom_pointrange(aes(ymin = q05, ymax = q95),
#                   position = position_dodge(width = 0.5),
#                   size = 0.7) +
#   geom_point(position = position_dodge(width = 0.5), size = 2.5) +
#   # Sample Size Text Centered Above All Points
#   geom_text(aes(label = SampleSize),
#             position = position_nudge(y = 0.02),  # Adjust this number for space
#             vjust = 0,  # Text below bottom edge = 1, centered = 0.5, above = 0
#             hjust = 0.5,
#             color = "black", size = 4) +
#   
#   scale_color_manual(values = c("Control" = "darkgreen", "Disease" = "purple")) +
#   labs(x = "Sex Bias Category", y = "Tissue Specificity", color = "Gene Type") +
#   guides(color = guide_legend(override.aes = list(size = 1))) +
#   theme_minimal(base_size = 16) +
#   theme(axis.text = element_text(size = 14),
#         axis.title = element_text(size = 16, face = "bold"),
#         legend.text = element_text(size = 14),
#         legend.title = element_text(size = 16),
#         strip.text = element_text(size = 16),
#         plot.margin = margin(10, 10, 10, 10),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.border = element_rect(color = "black", fill = NA, size = 1),
#         panel.background = element_blank()) +
#   geom_vline(xintercept = 1.5, color = "black", size = 1.2) +
#   geom_vline(xintercept = 2.5:5.5, color = "grey60", size = 0.5) +
#   scale_x_discrete(labels = c("All", "Highly FB", "Female-Biased", "Unbiased", "Male-Biased", "Highly MB"))
# ggsave("tissueSpecificity_plot2.pdf", width = 12, height = 7, dpi = 300)

#___________________________-----
# STATISTICAL ANALYSES----

# Permutation tests----

# Adding logical variable enables function to work----
summary_DoS_data_comparative$GeneType_logical <- summary_DoS_data_comparative$GeneType == "Control"
str(summary_DoS_data_comparative)

SBGE_DoS_bins_compare$GeneType_logical <- SBGE_DoS_bins_compare$GeneType == "Control"
str(SBGE_DoS_bins_compare)

summary_rMF_data_comparative$GeneType_logical <- summary_rMF_data_comparative$GeneType == "Control"
str(summary_rMF_data_comparative)

SBGE_rMF_bins_compare$GeneType_logical <- SBGE_rMF_bins_compare$GeneType == "Control"
str(SBGE_rMF_bins_compare)

summary_tajDS_data_comparative$GeneType_logical <- summary_tajDS_data_comparative$GeneType == "Control"
str(summary_tajDS_data_comparative)

SBGE_tajDS_bins_compare$GeneType_logical <- SBGE_tajDS_bins_compare$GeneType == "Control"
str(SBGE_tajDS_bins_compare)

summary_tissuespec_data_comparative$GeneType_logical <- summary_tissuespec_data_comparative$GeneType == "Control"
str(summary_tissuespec_data_comparative)

SBGE_tissuespec_bins_compare$GeneType_logical <- SBGE_tissuespec_bins_compare$GeneType == "Control"
str(SBGE_tissuespec_bins_compare)

summary_rlnd_data_comparative$GeneType_logical <- summary_rlnd_data_comparative$GeneType == "Control"
str(summary_rlnd_data_comparative)

SBGE_rlnd_bins_compare$GeneType_logical <- SBGE_rlnd_bins_compare$GeneType == "Control"
str(SBGE_rlnd_bins_compare)

summary_gonads_data_comparative$GeneType_logical <- summary_gonads_data_comparative$GeneType == "Control"
str(summary_gonads_data_comparative)

SBGE_gonads_bins_compare$GeneType_logical <- SBGE_gonads_bins_compare$GeneType == "Control"
str(SBGE_gonads_bins_compare)

summary_heads_data_comparative$GeneType_logical <- summary_heads_data_comparative$GeneType == "Control"
str(summary_heads_data_comparative)

SBGE_heads_bins_compare$GeneType_logical <- SBGE_heads_bins_compare$GeneType == "Control"
str(SBGE_heads_bins_compare)

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

# Tissue Specificity----
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

# Tissue SBGE ----
# GONADS

# Now run TwoPerm
result_gonads<- TwoPerm(
  perm_dat = SBGE_gonads_bins_compare,
  x_col = "log2FC_Gonads",
  groupBy = "GeneType_logical"
)

# Now run TwoPerm_SBGE() — it will include "All"
result_gonads_SBGE_cat <- TwoPerm_SBGE(
  perm_dat = SBGE_gonads_bins_compare,
  x_col = "log2FC_Gonads",
  groupBy = "GeneType_logical",
  SBGE_cat = "Sex.Bias.Category"
)

result_gonads_all <- data.frame(
  Sex.Bias.Category = "All",
  obs.diff = as.numeric(result_gonads[1]),
  n_TRUE = as.numeric(result_gonads[2]),
  n_FALSE = as.numeric(result_gonads[3]),
  pval = as.character(result_gonads[4]),
  Sig = as.character(result_gonads[5])
)

gonads_perm_results<- bind_rows(result_gonads_all , result_gonads_SBGE_cat )

# HEADS

# Now run TwoPerm
result_heads<- TwoPerm(
  perm_dat = SBGE_heads_bins_compare,
  x_col = "log2FC_Heads",
  groupBy = "GeneType_logical"
)

# Now run TwoPerm_SBGE() — it will include "All"
result_heads_SBGE_cat <- TwoPerm_SBGE(
  perm_dat = SBGE_heads_bins_compare,
  x_col = "log2FC_Heads",
  groupBy = "GeneType_logical",
  SBGE_cat = "Sex.Bias.Category"
)

result_heads_all <- data.frame(
  Sex.Bias.Category = "All",
  obs.diff = as.numeric(result_heads[1]),
  n_TRUE = as.numeric(result_heads[2]),
  n_FALSE = as.numeric(result_heads[3]),
  pval = as.character(result_heads[4]),
  Sig = as.character(result_heads[5])
)

heads_perm_results<- bind_rows(result_heads_all , result_heads_SBGE_cat )

gonads_perm_results$Tissue <- "Gonads"
heads_perm_results$Tissue <- "Heads"
combined_perm_tissue_results<- rbind(gonads_perm_results, heads_perm_results)



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

#rmf boot----
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

#Tissue boot----
#Gonads
gonads_boot<-TwoBoot(SBGE_gonads_bins_compare, x_col = "log2FC_Gonads", groupBy = "GeneType", boot_n = 10)
gonads_boot_SBGE <-TwoBoot_SBGE(SBGE_gonads_bins_compare, x_col = "log2FC_Gonads", groupBy = "GeneType",SBGE_cat = "Sex.Bias.Category")

# Add "All" category to first dataframe
gonads_boot$Sex.Bias.Category <- "All"
names(gonads_boot_SBGE)[names(gonads_boot_SBGE) == "Sig"] <- "GeneType"

gonads_boot_results<- bind_rows(gonads_boot, gonads_boot_SBGE )

#Heads
heads_boot<-TwoBoot(SBGE_heads_bins_compare, x_col = "log2FC_Heads", groupBy = "GeneType", boot_n = 10)
heads_boot_SBGE <-TwoBoot_SBGE(SBGE_heads_bins_compare, x_col = "log2FC_Heads", groupBy = "GeneType",SBGE_cat = "Sex.Bias.Category")

# Add "All" category to first dataframe
heads_boot$Sex.Bias.Category <- "All"
names(heads_boot_SBGE)[names(heads_boot_SBGE) == "Sig"] <- "GeneType"

heads_boot_results<- bind_rows(heads_boot, heads_boot_SBGE )

# Enrichment Test----


library(biomaRt)

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
# Prepare dataset

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
  scale_color_manual(values = c("Control" = "darkgreen", "Disease" = "purple")) +
  
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

ggsave("tajD_plot.pdf", width = 13, height = 6, dpi = 300)

# --- rMF Plot----
summary_rMF_data_comparative <- left_join(summary_rMF_data_comparative, rmf_boot_results, by = c("Sex.Bias.Category", "GeneType"))

summary_rMF_data_comparative$Sex.Bias.Category <- factor(summary_rMF_data_comparative$Sex.Bias.Category, 
                                                         levels = c("All","Highly female-biased", "Moderately female-biased", 
                                                                    "Unbiased", "Moderately male-biased", "Highly male-biased"))
rmf_plot <- ggplot(summary_rMF_data_comparative, 
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
  scale_color_manual(values = c("Control" = "darkgreen", "Disease" = "purple")) +
  
  # Axis labels
  labs(x = "Sex Bias Category", y = expression(italic(r)[mf]), color = "Gene Type") +
  
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

ggsave("rmf_plot.pdf", width = 13, height = 6, dpi = 300)

# DoS Plot----

summary_DoS_data_comparative <- left_join(summary_DoS_data_comparative, DoS_boot_results, by = c("Sex.Bias.Category", "GeneType"))

summary_DoS_data_comparative$Sex.Bias.Category <- factor(summary_DoS_data_comparative$Sex.Bias.Category, 
                                                         levels = c("All","Highly female-biased", "Moderately female-biased", 
                                                                    "Unbiased", "Moderately male-biased", "Highly male-biased"))


DoS_plot <- ggplot(summary_DoS_data_comparative, 
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
  scale_color_manual(values = c("Control" = "darkgreen", "Disease" = "purple")) +
  
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

ggsave("DoS_plot.pdf", width = 13, height = 6, dpi = 300)

# Tissue Specificity Plot ----
summary_tissuespec_data_comparative<- left_join(summary_tissuespec_data_comparative, TissueSpecificity_boot_results, by = c("Sex.Bias.Category", "GeneType"))

summary_tissuespec_data_comparative$Sex.Bias.Category <- factor(summary_tissuespec_data_comparative$Sex.Bias.Category, 
                                                                levels = c("All","Highly female-biased", "Moderately female-biased", 
                                                                           "Unbiased", "Moderately male-biased", "Highly male-biased"))

tissuespec <- ggplot(summary_tissuespec_data_comparative, 
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
  scale_color_manual(values = c("Control" = "darkgreen", "Disease" = "purple")) +
  
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

ggsave("tissuespec.pdf", width = 13, height = 6, dpi = 300)


#Relative Level of Nonsynonymous Diveristy Plot ----

summary_rlnd_data_comparative <- left_join(summary_rlnd_data_comparative, rlnd_boot_results, by = c("Sex.Bias.Category", "GeneType"))

summary_rlnd_data_comparative$Sex.Bias.Category <- factor(summary_rlnd_data_comparative$Sex.Bias.Category, 
                                                         levels = c("All","Highly female-biased", "Moderately female-biased", 
                                                                    "Unbiased", "Moderately male-biased", "Highly male-biased"))

rlnd_plot <- ggplot(summary_rlnd_data_comparative, 
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
  scale_color_manual(values = c("Control" = "darkgreen", "Disease" = "purple")) +
  
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

ggsave("rlnd_plot.pdf", width = 13, height = 6, dpi = 300)


#TISSUE----

#Gonads
summary_gonads_data_comparative <- left_join(summary_gonads_data_comparative, gonads_boot_results, by = c("Sex.Bias.Category", "GeneType"))

summary_gonads_data_comparative$Sex.Bias.Category <- factor(summary_gonads_data_comparative$Sex.Bias.Category, 
                                                          levels = c("All","Highly female-biased", "Moderately female-biased", 
                                                                     "Unbiased", "Moderately male-biased", "Highly male-biased"))
# Plot
ggplot(summary_gonads_data_comparative, 
       aes(x = Sex.Bias.Category, y = q50, color = GeneType)) + 
  geom_pointrange(aes(ymin = q05, ymax = q95),
                  position = position_dodge(width = 0.5),
                  size = 1) +  # Adjust error bar thickness
  geom_point(position = position_dodge(width = 0.5), size = 2) +  # Set size of points
  geom_text(aes(label = SampleSize),  
            position = position_dodge(width = 0.5),
            color = "black", size = 4) + # Adjust label position and size
  scale_color_manual(values = c("Control" = "darkgreen", "Disease" = "purple")) +  # Match colors
  labs(x = "SBGE",
       y = "SBGE in Gonads") +
  guides(color = guide_legend(override.aes = list(size = 1))) +
  theme_minimal(base_size = 16) +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 16, face = "bold"),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16),
        strip.text = element_text(size = 16),
        plot.margin = margin(10, 10, 10, 10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        panel.background = element_blank()) +
  geom_vline(xintercept = 1.5, color = "black", size = 1.2) +
  geom_vline(xintercept = 2.5:5.5, color = "grey60", size = 0.5) +
  scale_x_discrete(labels = c("All", "Highly FB", "Female-Biased", "Unbiased", "Male-Biased", "Highly MB")) 

#Heads
summary_heads_data_comparative <- left_join(summary_heads_data_comparative, heads_boot_results, by = c("Sex.Bias.Category", "GeneType"))

summary_heads_data_comparative$Sex.Bias.Category <- factor(summary_heads_data_comparative$Sex.Bias.Category, 
                                                          levels = c("All","Highly female-biased", "Moderately female-biased", 
                                                                     "Unbiased", "Moderately male-biased", "Highly male-biased"))
# Plot
ggplot(summary_heads_data_comparative, 
       aes(x = Sex.Bias.Category, y = q50, color = GeneType)) + 
  geom_pointrange(aes(ymin = q05, ymax = q95),
                  position = position_dodge(width = 0.5),
                  size = 1) +  # Adjust error bar thickness
  geom_point(position = position_dodge(width = 0.5), size = 2) +  # Set size of points
  geom_text(aes(label = SampleSize),  
            position = position_dodge(width = 0.5),
            color = "black", size = 4) + # Adjust label position and size
  scale_color_manual(values = c("Control" = "darkgreen", "Disease" = "purple")) +  # Match colors
  labs(x = "SBGE",
       y = "SBGE in Heads") +
  guides(color = guide_legend(override.aes = list(size = 1))) +
  theme_minimal(base_size = 16) +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 16, face = "bold"),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16),
        strip.text = element_text(size = 16),
        plot.margin = margin(10, 10, 10, 10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        panel.background = element_blank()) +
  geom_vline(xintercept = 1.5, color = "black", size = 1.2) +
  geom_vline(xintercept = 2.5:5.5, color = "grey60", size = 0.5) +
  scale_x_discrete(labels = c("All", "Highly FB", "Female-Biased", "Unbiased", "Male-Biased", "Highly MB")) 

summary_gonads_data_comparative$Tissue <- "Gonads"
summary_heads_data_comparative$Tissue <- "Heads"
combined_data_tissue <- rbind(summary_gonads_data_comparative, summary_heads_data_comparative)


combined_data_tissue$Sex.Bias.Category <- factor(combined_data_tissue$Sex.Bias.Category, 
                                          levels = c("All", "Highly female-biased", "Moderately female-biased", 
                                                     "Unbiased", "Moderately male-biased", "Highly male-biased"))


ggplot(combined_data_tissue, aes(x = Sex.Bias.Category, y = q50, color = GeneType)) + 
  geom_pointrange(aes(ymin = q05, ymax = q95), position = position_dodge(width = 0.5), size = 1) +
  geom_point(position = position_dodge(width = 0.5), size = 2) +
  geom_text(aes(label = SampleSize), position = position_dodge(width = 0.5), color = "black", size = 4) +
  scale_color_manual(values = c("Control" = "darkgreen", "Disease" = "purple")) +
  labs(x = "Sex Bias Category", y = "SBGE (log2FC)", color = "Gene Type") +
  facet_wrap(~Tissue, nrow = 1, scales = "fixed") +  # One column, vertical layout
  theme_minimal(base_size = 16) +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 16, face = "bold"),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16),
        strip.text = element_text(size = 16),
        plot.margin = margin(10, 10, 10, 10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        panel.background = element_blank()) +
  geom_vline(xintercept = 1.5, color = "black", size = 1.2) +
  geom_vline(xintercept = 2.5:5.5, color = "grey60", size = 0.5) +
  scale_x_discrete(labels = c("All", "Highly FB", "Female-Biased", "Unbiased", "Male-Biased", "Highly MB"))


tissue_expr<-ggplot(combined_data_tissue, aes(x = Sex.Bias.Category, y = q50, color = GeneType, shape = GeneType)) + 
  # Points and error bars
  geom_pointrange(aes(ymin = q05, ymax = q95), 
                  position = position_dodge(width = 0.5), 
                  size = 1.2)+
  geom_point(aes(shape = GeneType),
             position = position_dodge(width = 0.5),
             size = 2) +  
  geom_errorbar(aes(ymin = q05, ymax = q95), 
                position = position_dodge(width = 0.5),
                width = 0.4, 
                size = 1) +  # Thicker error bars
  # Sample size labels
  geom_text(aes(y = q95 + 0.8, label = SampleSize), 
            position = position_dodge(width = 0.7),
            color = "black", size = 4) +
  
  # Manual colors & shapes
  scale_color_manual(values = c("Control" = "darkgreen", "Disease" = "purple")) +
  scale_shape_manual(values = c("Control" = 16, "Disease" = 16)) +
  
  # Axis labels
  labs(x = "Sex Bias Category", y = "SBGE (log2FC)", color = "Gene Type", shape = "Gene Type") +
  
  # Facet layout
  facet_wrap(~Tissue, nrow = 1, strip.position = "top") +
  
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
    strip.text = element_text(size = 18, face = "bold", color = "#222222"),
    strip.background = element_rect(fill = "grey90", color = "black"),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 1.2),
    panel.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(20, 20, 20, 20)
  )

ggsave("tissue_expr.png", width = 14, height = 7, dpi = 300)

#Combining 4 evolutionary signatures ----

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


ggsave("combined_plot.pdf", combined_plot, width = 10, height = 17, dpi = 300)
