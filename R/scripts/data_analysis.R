#__________________________----
#PACKAGES ----
library(tidyverse) # Tidy data packages
library(janitor) # Clean data names
library(dplyr) #Data manipulation.
library(stringr) # Simplifies string operations
library(tidyr) # Tidying Data
library(biomaRt)
library(pheatmap)
library(ggplot2)
library(igraph)
library(ggraph)
library(tidygraph)
library(clipr)
#__________________________----
#IMPORT DATA ----

# Imported list of flybase gene ID text files from sex specific diseases
disease_files <- list.files(path = "~/Documents/Omics/Final_Year_Project/R/data/Sex specific diseases", pattern = "*.txt")
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


# Other diseases (literature search)
PCOS_orthologs <- read.csv("~/Documents/Omics/Final_Year_Project/R/data/diopt_results_PCOS_genes.csv")
colnames(PCOS_orthologs)[8] <- "FlyBaseID"

human_SA_orthologs <- read.csv("~/Documents/Omics/Final_Year_Project/R/data/diopt_results.csv")
colnames(human_SA_orthologs)[8] <- "FlyBaseID"

CAD_orthologs <- read.csv("~/Documents/Omics/Final_Year_Project/R/data/coronary_artery_disease_orthologs.csv")
colnames(CAD_orthologs)[8] <- "FlyBaseID"

osteoporosis_orthologs <- read.csv("~/Documents/Omics/Final_Year_Project/R/data/osteoporosis_orthologs.csv")
colnames(osteoporosis_orthologs)[8] <- "FlyBaseID"

AD_orthologs <- read.csv("~/Documents/Omics/Final_Year_Project/R/data/alzheimer_disease_orthologs.csv")
colnames(AD_orthologs)[8] <- "FlyBaseID"

parkinsons_orthologs <- read.csv("~/Documents/Omics/Final_Year_Project/R/data/parkinsons_orthologs.csv")
colnames(parkinsons_orthologs)[8] <- "FlyBaseID"


# Imported dataset for Mishra et al. 2024
# SBGE contrasts
SBG_whole_bodies <-read.csv("~/Documents/Omics/Final_Year_Project/R/data/Sex.Biased.Genes.Whole.csv")
colnames(SBG_whole_bodies)[1] <- "FlyBaseID"

# SBG_gonads <-read.csv("~/Documents/Omics/Final_Year_Project/R/data/Sex.Biased.Genes.Gonads.csv")
# colnames(SBG_gonads)[1] <- "FlyBaseID"
# 
# SBG_heads <-read.csv("~/Documents/Omics/Final_Year_Project/R/data/Sex.Biased.Genes.Heads.csv")
# colnames(SBG_heads)[1] <- "FlyBaseID"

# ASE contrasts
# DGE_whole_bodies <- read.csv("/Documents/Omics/Final_Year_Project/R/data/DGE_whole_bodies.csv")
# DGE_gonads <- read.csv("/Documents/Omics/Final_Year_Project/R/data/DifferentialGeneExpressionGonads.csv")
# DGE_heads <- read.csv("/Documents/Omics/Final_Year_Project/R/data/DifferentialGeneExpression.csv")

# Imported dataset from Fraisse et al. 2018
Fra <- read.delim(file="~/Documents/Omics/Final_Year_Project/R/data/FileS7.txt", header=TRUE)
colnames(Fra)[2] <- "FlyBaseID"
Fra[sapply(Fra, is.infinite)] <- NA

# Imported dataset from Singh and Agrawal 2023
Singh <- read.csv(file="~/Documents/Omics/Final_Year_Project/R/data/SBGEandSSSdataForMBE.csv", header=TRUE)
colnames(Singh)[2] <- "FlyBaseID"


# Imported dataset from Innocenti and Morrow 2010
# ant_genes <- read.csv("~/Documents/Omics/Final_Year_Project/R/data/Innocenti Antagonistic genes.csv")
# female_fitness_genes<- read.csv("~/Documents/Omics/Final_Year_Project/R/data/Innocenti female fitness ass. genes.csv")
# male_fitness_genes<- read.csv("~/Documents/Omics/Final_Year_Project/R/data/Innocenti male fitness ass. genes.csv")


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




#__________________________----
#CLEAN_DATA ----

# Filtering for the most accurate orthologs ----
filtered_orthologs <-all_orthologs %>%
  filter(
    `Rank` == "high",
    `Best.Score` == "Yes",
    `Best.Score.Reverse` == "Yes",
    `DIOPT.Score` >= 3)

orthologs <- filtered_orthologs[,c("Human.Symbol","Fly.Symbol", "FlyBaseID", "disease")]



# Assigning gene name to flybase geneIDs ----
ensembl <- useMart("ENSEMBL_MART_ENSEMBL", dataset = "dmelanogaster_gene_ensembl")

FlyBaseID <-all_diseases_data$FlyBaseID

gene_info <- getBM(
  attributes = c("flybase_gene_id","external_gene_name"),
  filters = "flybase_gene_id",
  values = FlyBaseID,
  mart = ensembl
)


######################
# 
# Stopped here... You'll have to rename the gene_info columns to match the all_diseases dataset before merging
# 
# There may be some other column and data frame naming that needs updating or fixing
# 
######################

# This could just be overwrite instead of a new df
combined_data <- merge(gene_info, all_diseases_data, by = "FlyBaseID", all = TRUE) # Keep all genes, even if no disease info










# Now get chromosome locations from reference genome (will be useful later)
##########

all.genes <- read.delim(file="~/Documents/BIO UEA/Teaching/Module - Data Science and Bioinformatics/Session10/RNA-Seq_Tutorial/ExtractGenes/GeneLists/all.genes.tsv", sep="\t", header=FALSE)
colnames(all.genes) = c("FlyBaseID")

Xchr <- read.delim(file="~/Documents/BIO UEA/Teaching/Module - Data Science and Bioinformatics/Session10/RNA-Seq_Tutorial/ExtractGenes/GeneLists/X.chromosome.genes.tsv", sep="\t", header=TRUE)
colnames(Xchr) = c("FlyBaseID")
Xchr$Chr <- rep("X", dim(Xchr)[1])

Ychr <- read.delim(file="~/Documents/BIO UEA/Teaching/Module - Data Science and Bioinformatics/Session10/RNA-Seq_Tutorial/ExtractGenes/GeneLists/Y.chromosome.genes.tsv", sep="\t", header=TRUE)
colnames(Ychr) = c("FlyBaseID")
Ychr$Chr <- rep("Y", dim(Ychr)[1])

chr2L <- read.delim(file="~/Documents/BIO UEA/Teaching/Module - Data Science and Bioinformatics/Session10/RNA-Seq_Tutorial/ExtractGenes/GeneLists/2L.chromosome.genes.tsv", sep="\t", header=FALSE)
colnames(chr2L) = c("FlyBaseID")
chr2R <- read.delim(file="~/Documents/BIO UEA/Teaching/Module - Data Science and Bioinformatics/Session10/RNA-Seq_Tutorial/ExtractGenes/GeneLists/2R.chromosome.genes.tsv", sep="\t", header=FALSE)
colnames(chr2R) = c("FlyBaseID")
#
chr2 <- rbind(chr2L, chr2R)
chr2$Chr <- rep("2", dim(chr2)[1])

chr3L <- read.delim(file="~/Documents/BIO UEA/Teaching/Module - Data Science and Bioinformatics/Session10/RNA-Seq_Tutorial/ExtractGenes/GeneLists/3L.chromosome.genes.tsv", sep="\t", header=FALSE)
colnames(chr3L) = c("FlyBaseID")
chr3R <- read.delim(file="~/Documents/BIO UEA/Teaching/Module - Data Science and Bioinformatics/Session10/RNA-Seq_Tutorial/ExtractGenes/GeneLists/3R.chromosome.genes.tsv", sep="\t", header=FALSE)
colnames(chr3R) = c("FlyBaseID")
#
chr3 <- rbind(chr3L, chr3R)
chr3$Chr <- rep("3", dim(chr3)[1])

Chrs <- rbind(Xchr, Ychr, chr2, chr3) # Not all genes; just X, Y, 2, and 3.
Chrs$Chr <- as.factor(Chrs$Chr)
##########







##########
#
# GO enrichment for candidate differentially expressed genes
# 
##########

# Target set
write_clip(disease.ortholog.master.list$FlyBaseID) # You decide which genes to use as target set
# Background set
write_clip(Chrs$FlyBaseID) 
# Background set
# These were used here: https://biit.cs.ut.ee/gprofiler/gost
# Alternatively, use: https://bioinformatics.sdstate.edu/go/
##########



















# Renaming variables to match ortholog dataset ----

combined_data_renamed <- combined_data %>%
  rename("Fly.Symbol" = "external_gene_name") %>%
  mutate("Human.Symbol" = NA) %>%  # Add missing column with NA
  select("Human.Symbol" , "Fly.Symbol", "FlyBaseID", "disease")


final_dataset <- bind_rows(orthologs, combined_data_renamed)%>%
  rename("FlyBaseID" = "FlyBaseID")


# Merging data frames with Fraisse et al.2023 ----

Fra_master_data <- Fra %>%
  inner_join(final_dataset, by = "FlyBaseID")%>%
  select(disease, everything())  # Reorder columns

# Remove exact duplicates
Fra_master_data  <- Fra_master_data  %>% distinct()

# Merging data frames with Mishra et al. 2024 ----

SBG_gonads_data <- SBG_gonads %>% mutate(Tissue = "Gonads")

SBG_heads_data <- SBG_heads%>% mutate(Tissue = "Heads")

SBG_wb_data <- SBG_whole_bodies %>% mutate(Tissue = "Whole Body")

SBG_data <- bind_rows(SBG_gonads_data, SBG_heads_data, SBG_wb_data)

DGE_wb_master_data <- DGE_whole_bodies %>%
  inner_join(final_dataset, by = "FlyBaseID")%>%
  select(disease, everything())  # Reorder columns

SBG_master_data <- SBG_data %>%
  inner_join(final_dataset, by = "FlyBaseID")%>%
  select(disease, everything())  # Reorder columns

SBG_master_data <- SBG_master_data%>% distinct()

# Merged mishra Sex biased genes and Fraisse et al. 2018

mish_fra_SBG<- SBG_master_data  %>%
  inner_join(Fra, by = "FlyBaseID")%>%
  select(disease, everything())  # Reorder columns

SBG_no_disease<- mish_fra_SBG  %>%
  select(FlyBaseID,Sex.Bias.Category,log2FoldChange,Tissue,DoS)  # Reorder columns

SBG_no_disease <- SBG_no_disease%>% distinct()
#__________________________----
#DATA VISUALS ----
 
# Removing repeats where disease listed twice 

SBG_no_disease_column <- mish_fra_SBG  %>%
  select(FlyBaseID,Sex.Bias.Category,log2FoldChange,Tissue,DoS)  # Reorder columns

SBG_no_disease_column <- SBG_no_disease%>% distinct()

# Creating Bins
bins <- SBG_no_disease %>%
  mutate(Sex.Bias.Category = case_when(
    log2FoldChange >= 2 ~ "Highly male-biased",
    log2FoldChange >= 0.5 & log2FoldChange < 2 ~ "Moderately male-biased",
    log2FoldChange >= -0.5 & log2FoldChange < 0.5 ~ "Unbiased",
    log2FoldChange >= -2 & log2FoldChange < -0.5 ~ "Moderately female-biased",
    log2FoldChange < -2 ~ "Highly female-biased"
  ))

color_map <- c("Highly female-biased" = "red", "Moderately female-biased" = "pink", 
               "Unbiased" = "black", "Moderately male-biased" = "lightblue", "Highly male-biased" = "blue")

bins$Sex.Bias.Category <- factor(bins $Sex.Bias.Category, 
                                 levels = c("Highly female-biased", "Moderately female-biased", 
                                            "Unbiased", "Moderately male-biased", "Highly male-biased"))


# Compute Mean and Standard Error (SE) for DoS per Sex Bias Category
summary_data <- bins %>%
  group_by(Sex.Bias.Category) %>%
  summarise(
    Mean_DoS = mean(DoS, na.rm = TRUE),
    SE_DoS = sd(DoS, na.rm = TRUE) / sqrt(n())  # Standard Error
  )

# Creating Summary table for tissue data
summary_tissue_data <- bins %>%
  group_by(Sex.Bias.Category, Tissue) %>%
  summarise(
    Mean_DoS = mean(DoS, na.rm = TRUE),
    SE_DoS = sd(DoS, na.rm = TRUE) / sqrt(n()),  # Standard Error
    .groups = "drop"
  )

# Sex Bias/DoS ----

ggplot(summary_data, aes(x = Sex.Bias.Category, y = Mean_DoS, color = Sex.Bias.Category)) +
  geom_point(size = 4) +  # Plot mean DoS
  geom_errorbar(aes(ymin = Mean_DoS - SE_DoS, ymax = Mean_DoS + SE_DoS), width = 0.2) +  # Error bars
  scale_color_manual(values = color_map) +  # Custom colors
  labs(title = "", x = "", y = expression(italic(DoS)))+
  
  theme_bw()+
  theme(
    axis.text.x = element_text(size = 12),
    axis.title.y = element_text(size = 14, face = "italic"),
    panel.grid = element_blank(),
    legend.position = "right"
  )+
  
  # Add bins for each category
  geom_vline(xintercept = 1.5, linetype = "solid", color = "black") +
  geom_vline(xintercept = 2.5, linetype = "solid", color = "black") +
  geom_vline(xintercept = 3.5, linetype = "solid", color = "black") +
  geom_vline(xintercept = 4.5, linetype = "solid", color = "black")


summary_tissue_data$Sex.Bias.Category <- factor(summary_tissue_data$Sex.Bias.Category, 
                                         levels = c("Highly female-biased", "Moderately female-biased", 
                                                     "Unbiased", "Moderately male-biased", "Highly male-biased"))

tissue_colors <- c("Gonads" = "purple", "Heads" = "black", "Whole Body" = "gray")

# Sex Bias/DoS -Tissue type ----

summary_tissue_data$Sex.Bias.Category <- factor(summary_tissue_data$Sex.Bias.Category, 
                                                levels = c("Highly female-biased", "Moderately female-biased", 
                                                           "Unbiased", "Moderately male-biased", "Highly male-biased"))

tissue_colors <- c("Gonads" = "purple", "Heads" = "black", "Whole Body" = "gray")


dos_tissue_plot<- ggplot(summary_tissue_data, aes(x = Sex.Bias.Category, y = Mean_DoS, color = Tissue)) +
  geom_point(position = position_dodge(width = 0.6),
            size = 3) +  # Plot mean DoS without dodge
  geom_errorbar(aes(ymin = Mean_DoS - SE_DoS, ymax = Mean_DoS + SE_DoS), width = 0.2, position = position_dodge(width = 0.6)) +  # Error bars
  scale_color_manual(values = tissue_colors) +  # Colors for tissues
  theme_bw() +
  labs(title = "", x = "", y = expression(italic(DoS)), color = "Tissue Type") +
  geom_abline(intercept = 0, slope = 0,  size = 0.5, linetype= "solid", color = "black") +
  theme_classic() +
  theme(plot.title.position = c("panel"),
        legend.title = element_blank(),
        legend.position = c("bottom"),
        legend.box.background = element_rect(),
        legend.text = element_text(size = 20, color = "black"),
        axis.text.x = element_text(size=20, margin = margin(5,0,0,0), color = "black"),
        axis.text.y = element_text(size=20, margin = margin(0,5,0,0), color = "black"), 
        axis.title.x = element_text(size=40, margin = margin(10,0,0,0), color = "black"),
        axis.title.y = element_text(size=40, margin = margin(0,10,0,0), color = "black"),
        plot.title = element_text(size=40, margin = margin(0,0,0,0), color = "black"),
        plot.margin = margin(6,6,6,6)
  )
  
print(dos_tissue_plot)

# Expression Level Whole body


SBG_no_disease_column <- SBG_no_disease%>% distinct()

SGE_expression<- SBG_whole_bodies%>%
  inner_join(Fra_master_data, by = "FlyBaseID")%>%
  select(FlyBaseID,Sex.Bias.Category,log2FoldChange,exprLevelWB)%>% 
  distinct() # Remove disease variable to avoid repeats

# Creating Bins
bins_expr <-SGE_expression %>%
  mutate(Sex.Bias.Category = case_when(
    log2FoldChange >= 2 ~ "Highly male-biased",
    log2FoldChange >= 0.5 & log2FoldChange < 2 ~ "Moderately male-biased",
    log2FoldChange >= -0.5 & log2FoldChange < 0.5 ~ "Unbiased",
    log2FoldChange >= -2 & log2FoldChange < -0.5 ~ "Moderately female-biased",
    log2FoldChange < -2 ~ "Highly female-biased"
  ))


summary_expr_data <- bins_expr %>%
  group_by(Sex.Bias.Category) %>%
  summarise(
    Mean_expr = mean(exprLevelWB, na.rm = TRUE),
    SE_expr = sd(exprLevelWB, na.rm = TRUE) / sqrt(n())  # Standard Error
  )

summary_expr_data$Sex.Bias.Category <- factor(summary_expr_data$Sex.Bias.Category, 
                                                levels = c("Highly female-biased", "Moderately female-biased", 
                                                           "Unbiased", "Moderately male-biased", "Highly male-biased"))


expr_plot<- ggplot(summary_expr_data, aes(x = Sex.Bias.Category, y = Mean_expr)) +
  geom_point(position = position_dodge(width = 0.6),
             size = 3) +  # Plot mean DoS without dodge
  geom_errorbar(aes(ymin = Mean_expr - SE_expr, ymax = Mean_expr+ SE_expr), width = 0.2, position = position_dodge(width = 0.6)) +  # Error bars
  theme_bw() +
  labs(title = "", x = "", y = "Level of Whole body expression") +
  geom_abline(intercept = 0, slope = 0,  size = 0.5, linetype= "solid", color = "black") +
  theme_classic() +
  theme(plot.title.position = c("panel"),
        legend.title = element_blank(),
        legend.position = c("bottom"),
        legend.box.background = element_rect(),
        legend.text = element_text(size = 20, color = "black"),
        axis.text.x = element_text(size=20, margin = margin(5,0,0,0), color = "black"),
        axis.text.y = element_text(size=20, margin = margin(0,5,0,0), color = "black"), 
        axis.title.x = element_text(size=40, margin = margin(10,0,0,0), color = "black"),
        axis.title.y = element_text(size=20, margin = margin(0,10,0,0), color = "black"),
        plot.title = element_text(size=40, margin = margin(0,0,0,0), color = "black"),
        plot.margin = margin(6,6,6,6)
  )

print(expr_plot)


# Expression Level Whole Body ---- 

SGE_expression$Sex.Bias.Category <- factor(SGE_expression$Sex.Bias.Category, 
                                 levels = c("Highly female-biased", "Moderately female-biased", 
                                            "Unbiased", "Moderately male-biased", "Highly male-biased"))

colour_map <- c("Highly female-biased" = "red", "Moderately female-biased" = "pink", 
               "Unbiased" = "black", "Moderately male-biased" = "lightblue", "Highly male-biased" = "blue")

ggplot(SGE_expression, aes(x = Sex.Bias.Category, y = exprLevelWB, color = Sex.Bias.Category)) +
  geom_point(size = 4, alpha = 0.7) +  # Scatter points
  geom_smooth(method = "lm", se = FALSE, linetype = "dashed") +  # Add regression line
  scale_color_manual(values = colour_map) +  # Apply custom color scheme
  labs(
    title = "Sex Bias Gene expression of whole body",
    x = "",
    y = "Expression Level (WB)",
    color = "Sex Bias Category"
  ) +
  theme_minimal()

ggplot(SGE_expression, aes(x = Sex.Bias.Category, y = exprLevelWB, fill = Sex.Bias.Category)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +  # Outlier points removed for clarity
  geom_jitter(aes(color = Sex.Bias.Category), size = 2, width = 0.2, alpha = 0.7) +  # Add jittered points
  scale_fill_manual(values = colour_map) +  # Apply custom fill colors for boxes
  scale_color_manual(values = colour_map) +  # Apply custom colors for jittered points
  labs(
    title = "Expression Level (WB) Across Sex Bias Categories",
    x = "Sex Bias Category",
    y = "Expression Level (WB)",
    fill = "Sex.Bias.Category"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels for better readability


# Tajimas D


# Go enrichment ----

write_clip(union_congruentTogether$FlyBaseID) # Target set
write_clip(data$FlyBaseID)












