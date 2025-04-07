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
library(plyr)
library(broom)
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

# Allele specific contrasts

DGE_whole_bodies <- read.csv("/Users/elliehowatson/Final_Year_Project/R/data/DGE_whole_bodies.csv")

DGE_whole_bodies <- DGE_whole_bodies%>% 
  mutate(Tissue = "Whole Body")

DGE_gonads <- read.csv("/Users/elliehowatson/Final_Year_Project/R/data/DifferentialGeneExpressionGonads.csv")

DGE_gonads <- DGE_gonads%>% 
  mutate(Tissue = "Gonads")

DGE_heads <- read.csv("/Users/elliehowatson/Final_Year_Project/R/data/DifferentialGeneExpression.csv")

DGE_heads <- DGE_heads%>% 
  mutate(Tissue = "Heads")

# Bind DGE Dataframe Together

DGE_data<- rbind(DGE_whole_bodies, DGE_gonads, DGE_heads)

# Select and rename desired columns of Differential gene expression

DGE_data <- data.frame(cbind(DGE_data$log2FoldChange, DGE_data$lfcSE, DGE_data$FlyBaseID, DGE_data$Tissue))
colnames(DGE_data) <- c("exp_SBGE_ase", "se_SBGE_ase", "FlyBaseID", "Tissue")

# Convert columns to numeric
DGE_data$exp_SBGE_ase <- as.numeric(DGE_data$exp_SBGE_ase)
DGE_data$se_SBGE_ase <- as.numeric(DGE_data$se_SBGE_ase)
str(DGE_data)

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

gene_info <- gene_info %>%
  rename("Fly.Symbol" = "external_gene_name",
         "FlyBaseID" = "flybase_gene_id")


gene_info <- merge(gene_info, all_diseases_data, by = "FlyBaseID", all = TRUE) 

# Combining flybase search genes with ortholog genes ----

gene_info_inc_human<- gene_info %>%
  mutate("Human.Symbol" = NA) %>%  # Add missing column with NA
  select("Human.Symbol" , "Fly.Symbol", "FlyBaseID", "disease")

final_data_inc_human_genes<- bind_rows(orthologs, gene_info_inc_human)%>%
  rename("FlyBaseID" = "FlyBaseID")

final_data<- bind_rows(orthologs, gene_info)%>%
  select( "Fly.Symbol", "FlyBaseID", "disease")

#________________________----
#DATA ANALYSIS----

#Chromosome locations from reference genome ----

all.genes <- read.delim(file="~/Final_Year_Project/R/data/GeneLists/all.genes.tsv", sep="\t", header=FALSE)
colnames(all.genes) = c("FlyBaseID")

Xchr <- read.delim(file="~/Final_Year_Project/R/data/GeneLists/X.chromosome.genes.tsv", sep="\t", header=TRUE)
colnames(Xchr) = c("FlyBaseID")
Xchr$Chr <- rep("X", dim(Xchr)[1])

Ychr <- read.delim(file="~/Final_Year_Project/R/data/GeneLists/Y.chromosome.genes.tsv", sep="\t", header=TRUE)
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
write_clip(Chrs_disease $FlyBaseID) # You decide which genes to use as target set
# Background set
write_clip(Chrs$FlyBaseID) 

# Background set
# These were used here: https://biit.cs.ut.ee/gprofiler/gost
# Alternatively, use: https://bioinformatics.sdstate.edu/go/

# Chromosomal distribution ----

# Count occurrences of each chromosome
chr_counts <- Chrs_disease %>%
  group_by(Chr) %>%
  summarise(count = n())


# Creating Bin Categories ---- 

# # Define cut-offs for SBGE categorization
# x1 = 1 # first cut-off
# x2 = 5 # second cut-off
# y0 = 1 # tolerance of middle bins
# 
# # Complex (5 levels) SBGE categorization
# 
# # Highly female-biased gene expression
# more.fbg.keep <- DGE_data$exp_SBGE_ase < -x2 & (DGE_data$exp_SBGE_ase + DGE_data$se_SBGE_ase) < -(x2-y0)
# more.fbg <- DGE_data[more.fbg.keep,]
# more.fbg$SBGE_comp <- rep(c("a.more.fbg"), dim(more.fbg)[1])
# 
# # Female-biased gene expression
# fbg.keep <- DGE_data$exp_SBGE_ase < -x1 & DGE_data$exp_SBGE_ase > -x2 & (DGE_data$exp_SBGE_ase + DGE_data$se_SBGE_ase) < 0 & (DGE_data$exp_SBGE_ase - DGE_data$se_SBGE_ase) > -(x2+y0)
# fbg <- DGE_data[fbg.keep,]
# fbg$SBGE_comp <- rep(c("b.fbg"), dim(fbg)[1])
# 
# # Unbiased category
# ubg.keep <- DGE_data$exp_SBGE_ase < x1 & DGE_data$exp_SBGE_ase > -x1 & (DGE_data$exp_SBGE_ase - DGE_data$se_SBGE_ase) > -(x1+y0) & (DGE_data$exp_SBGE_ase + DGE_data$se_SBGE_ase) < (x1+y0)
# ubg <- DGE_data[ubg.keep,]
# ubg$SBGE_comp <- rep(c("c.ubg"), dim(ubg)[1])
# 
# # Male-biased gene expression
# mbg.keep <- DGE_data$exp_SBGE_ase > x1 & DGE_data$exp_SBGE_ase < x2 & (DGE_data$exp_SBGE_ase - DGE_data$se_SBGE_ase) > 0 & (DGE_data$exp_SBGE_ase + DGE_data$se_SBGE_ase) < (x2+y0)
# mbg <- DGE_data[mbg.keep,]
# mbg$SBGE_comp <- rep(c("d.mbg"), dim(mbg)[1])
# 
# # Highly male-biased gene expression
# more.mbg.keep <- DGE_data$exp_SBGE_ase > x2 & (DGE_data$exp_SBGE_ase - DGE_data$se_SBGE_ase) > (x2-y0)
# more.mbg <- DGE_data[more.mbg.keep,]
# more.mbg$SBGE_comp <- rep(c("e.more.mbg"), dim(more.mbg)[1])

# # Combine the categories
# comp_bins<- rbind(more.fbg, fbg, ubg, mbg, more.mbg)
# str(comp_bins)

# Merging data frames with Fraisse et al.2023 ----

Fra_master_data <- Fra %>%
  inner_join(final_data, by = "FlyBaseID")%>%
  select(disease, everything())  # Reorder columns

# Remove exact duplicates
Fra_master_data  <- Fra_master_data  %>% distinct()

# Merging data frames with Mishra et al. 2024 ----

# DGE_master_data <- DGE_data%>%
#   inner_join(final_data, by = "FlyBaseID")


SBG_wb_data <- SBG_whole_bodies %>%
  inner_join(final_data, by = "FlyBaseID")

SBG_wb_data <- SBG_wb_data %>% distinct()
colnames(SBG_wb_data) <- c("FlyBaseID","Sex.Bias.Category","exp_SBGE_ase","Tissue","Fly.Symbol","disease")

# Removing repeats where disease listed twice on DGE AND SBGE datasets 

SBG_no_disease_column <- SBG_wb_data %>%
  select(FlyBaseID,Sex.Bias.Category,exp_SBGE_ase,Tissue)# Reorder columns

SBG_no_disease_column <- SBG_no_disease_column %>% distinct() # Remove repeats

# DGE_no_disease_column <- DGE_master_data  %>%
#   select(FlyBaseID,exp_SBGE_ase,Tissue,Fly.Symbol)# Reorder columns
# 
# DGE_no_disease_column <- DGE_no_disease_column %>% distinct() # Remove repeats

# Merged mishra Sex biased genes and Fraisse et al. 2018 ----
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

# Control genes ----

# Create a control set by excluding disease genes

disease_gene_ids <- unique(Chrs_disease$FlyBaseID)

control_genes <- Chrs%>%
  filter(!(FlyBaseID %in% disease_gene_ids))

# Merging Control genes with Mishra Singh and Agrawal ----

control_singh<- control_genes %>%
  inner_join(Singh, by = "FlyBaseID")

SBGE_control_data<- rbind(SBG_whole_bodies, SBG_gonads, SBG_heads)

control_mish_SBGE<-SBGE_control_data%>%
  inner_join( control_genes, by = "FlyBaseID")

control_mish_SBGE_wb<-SBG_whole_bodies%>%
  inner_join( control_genes, by = "FlyBaseID")

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

SBGE_control_data<- rbind(SBG_whole_bodies, SBG_gonads, SBG_heads)

control_mish_SBGE_wb<-SBG_whole_bodies%>%
  inner_join( control_genes, by = "FlyBaseID")

# Merging control Mishra genes and control Singh genes

control_Fra_mish_SBGE_wb<-control_mish_SBGE_wb%>%
  inner_join(control_Fra, by = "FlyBaseID")

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

#__________________________----
#DATA VISUALS ----
 
# Removing repeats where disease listed twice 

SBG_no_disease_column <- mish_fra_SBG  %>%
  select(FlyBaseID,Sex.Bias.Category,exp_SBGE_ase,Tissue,DoS)  # Reorder columns

SBG_no_disease_column <- SBG_no_disease_column%>% distinct()

# Change into numerical 
str(SBG_no_disease_column$exp_SBGE_ase)
SBG_no_disease_column <- SBG_no_disease_column %>%
  mutate(exp_SBGE_ase = as.numeric(exp_SBGE_ase))

#Creating Bins ----
bins <- SBG_no_disease_column %>%
  mutate(Sex.Bias.Category = case_when(
    exp_SBGE_ase >= 2 ~ "Highly male-biased",
    exp_SBGE_ase >= 0.5 & exp_SBGE_ase < 2 ~ "Moderately male-biased",
    exp_SBGE_ase >= -0.5 & exp_SBGE_ase < 0.5 ~ "Unbiased",
    exp_SBGE_ase > -2 & exp_SBGE_ase <= -0.5 ~ "Moderately female-biased",
    exp_SBGE_ase <= -2 ~ "Highly female-biased"
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

# Creating Summary table for DoS tissue data
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

SBGE_expression<- SBG_whole_bodies%>%
  inner_join(Fra_master_data, by = "FlyBaseID")%>%  
  rename("exp_SBGE_ase" = "log2FoldChange")%>%
  select(FlyBaseID,Sex.Bias.Category,exp_SBGE_ase,exprLevelWB)%>% 
  distinct() # Remove disease variable to avoid repeats

# Creating Bins
bins_expr <-SBGE_expression %>%
  mutate(Sex.Bias.Category = case_when(
    exp_SBGE_ase >= 2 ~ "Highly male-biased",
    exp_SBGE_ase >= 0.5 & exp_SBGE_ase < 2 ~ "Moderately male-biased",
    exp_SBGE_ase >= -0.5 & exp_SBGE_ase < 0.5 ~ "Unbiased",
    exp_SBGE_ase >= -2 & exp_SBGE_ase < -0.5 ~ "Moderately female-biased",
    exp_SBGE_ase < -2 ~ "Highly female-biased"
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
# Average Expression Level Whole Body ---- 

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


hist(SBGE_expression$exprLevelWB)

SBGE_expression$Sex.Bias.Category <- factor(SBGE_expression$Sex.Bias.Category, 
                                 levels = c("Highly female-biased", "Moderately female-biased", 
                                            "Unbiased", "Moderately male-biased", "Highly male-biased"))

colour_map <- c("Highly female-biased" = "red", "Moderately female-biased" = "pink", 
               "Unbiased" = "black", "Moderately male-biased" = "lightblue", "Highly male-biased" = "blue")

ggplot(SBGE_expression, aes(x = Sex.Bias.Category, y = exprLevelWB, color = Sex.Bias.Category)) +
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

ggplot(SBGE_expression, aes(x = Sex.Bias.Category, y = exprLevelWB, fill = Sex.Bias.Category)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA, width = 0.6) + # Wider box plots
  geom_jitter(aes(color = Sex.Bias.Category), size = 2, alpha = 0.7, width = 0.15) + # Less horizontal jitter
  scale_fill_manual(values = colour_map) +
  scale_color_manual(values = colour_map) +
  scale_x_discrete(expand = c(0.2, 0.2)) + # Increase spacing
  labs(
    title = "Expression Level (WB) Across Sex Bias Categories",
    x = "Sex Bias Category",
    y = "Expression Level (WB)",
    fill = "Sex.Bias.Category"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) # Rotate x-axis labels

ggplot(SBGE_expression, aes(x = Sex.Bias.Category, y = exprLevelWB, fill = Sex.Bias.Category)) +
  geom_violin(alpha = 0.5, width = 0.8, draw_quantiles = c(0.25, 0.5, 0.75)) + # Adds quartiles
  geom_jitter(aes(color = Sex.Bias.Category), size = 1.5, width = 0.15, alpha = 0.7) + # Keeps jittered points
  scale_fill_manual(values = colour_map) +
  scale_color_manual(values = colour_map) +
  labs(
    title = "Expression Level (WB) Across Sex Bias Categories",
    x = "Sex Bias Category",
    y = "Expression Level (WB)"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(SBGE_expression, aes(x = Sex.Bias.Category, y = exprLevelWB, fill = Sex.Bias.Category)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA, width = 0.6) +
  geom_jitter(aes(color = Sex.Bias.Category), size = 2, width = 0.15, alpha = 0.7) +
  scale_y_log10() + # Apply log scale
  scale_fill_manual(values = colour_map) +
  scale_color_manual(values = colour_map) +
  labs(
    title = "Expression Level (WB) Across Sex Bias Categories (Log Scale)",
    x = "Sex Bias Category",
    y = "Expression Level (WB) (log scale)"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Tajimas D ----

# Change into numerical 

SBGE_tajDN<- mish_singh_SBG%>%
  select(FlyBaseID,Sex.Bias.Category,exp_SBGE_ase,tajD.N,Tissue)%>% 
  distinct() # Remove disease variable to avoid repeats

SBGE_tajDN<- SBGE_tajDN %>%
  mutate(exp_SBGE_ase = as.numeric(exp_SBGE_ase))

SBGE_tajDS<- mish_singh_SBG%>%
  select(FlyBaseID,Sex.Bias.Category,exp_SBGE_ase,tajD.S,Tissue)%>% 
  distinct() # Remove disease variable to avoid repeats

SBGE_tajDS<- SBGE_tajDS %>%
  mutate(exp_SBGE_ase = as.numeric(exp_SBGE_ase))

# Creating Bins

SBGE_tajDN_bins<-SBGE_tajDN %>%
  mutate(Sex.Bias.Category = case_when(
    exp_SBGE_ase >= 2 ~ "Highly male-biased",
    exp_SBGE_ase >= 0.5 & exp_SBGE_ase < 2 ~ "Moderately male-biased",
    exp_SBGE_ase >= -0.5 & exp_SBGE_ase < 0.5 ~ "Unbiased",
    exp_SBGE_ase >= -2 & exp_SBGE_ase < -0.5 ~ "Moderately female-biased",
    exp_SBGE_ase < -2 ~ "Highly female-biased"
  ))

SBGE_tajDS_bins<-SBGE_tajDS %>%
  mutate(Sex.Bias.Category = case_when(
    exp_SBGE_ase >= 2 ~ "Highly male-biased",
    exp_SBGE_ase >= 0.5 & exp_SBGE_ase < 2 ~ "Moderately male-biased",
    exp_SBGE_ase >= -0.5 & exp_SBGE_ase < 0.5 ~ "Unbiased",
    exp_SBGE_ase >= -2 & exp_SBGE_ase < -0.5 ~ "Moderately female-biased",
    exp_SBGE_ase < -2 ~ "Highly female-biased"
  ))


# Compute Mean and Standard Error (SE) for tajDN and tajDS per Sex Bias Category

summary_tajDN_data <-SBGE_tajDN_bins %>%
  group_by(Sex.Bias.Category) %>%
  summarise(
    Mean_tajD.N = mean(tajD.N, na.rm = TRUE),
    SE_tajD.N = sd(tajD.N, na.rm = TRUE) / sqrt(n())  # Standard Error
  )


summary_tajDS_data <-SBGE_tajDS_bins %>%
  group_by(Sex.Bias.Category) %>%
  summarise(
    Mean_tajD.S = mean(tajD.S, na.rm = TRUE),
    SE_tajD.S = sd(tajD.S, na.rm = TRUE) / sqrt(n())  # Standard Error
  )

# Merge summary tajD data
summary_tajD_data <- merge(summary_tajDN_data,summary_tajDS_data)



# Creating Summary table for tajD tissue data

SBGE_tajD_tissue<- mish_singh_SBG%>%
  select(FlyBaseID,Sex.Bias.Category,exp_SBGE_ase,tajD.S,tajD.N,Tissue)%>% 
  distinct() # Remove disease variable to avoid repeats

summary_tissue_tajDS_data <- SBGE_tajDS_bins %>%
  group_by(Sex.Bias.Category, Tissue) %>%
  summarise( 
    Mean_tajD.S = mean(tajD.S, na.rm = TRUE),
    SE_tajD.S = sd(tajD.S, na.rm = TRUE) / sqrt(n()),  # Standard Error
    .groups = "drop"
  )

summary_tissue_tajDN_data <- SBGE_tajDN_bins %>%
  group_by(Sex.Bias.Category, Tissue) %>%
  summarise( 
    Mean_tajD.N = mean(tajD.N, na.rm = TRUE),
    SE_tajD.N = sd(tajD.N, na.rm = TRUE) / sqrt(n()),  # Standard Error
    .groups = "drop"
  )

summary_tajD_tissue_data <- merge(summary_tissue_tajDN_data,summary_tissue_tajDS_data)

# Plot tajD non-synonymous

summary_tajDN_data$Sex.Bias.Category <- factor(summary_tajDN_data$Sex.Bias.Category, 
                                             levels = c("Highly female-biased", "Moderately female-biased", 
                                                        "Unbiased", "Moderately male-biased", "Highly male-biased"))

ggplot(summary_tajDN_data, aes(x = Sex.Bias.Category, y = Mean_tajD.N, color = Sex.Bias.Category)) +
  geom_point(size = 4) +  # Plot mean DoS
  geom_errorbar(aes(ymin = Mean_tajD.N - SE_tajD.N, ymax = Mean_tajD.N + SE_tajD.N), width = 0.2) +  # Error bars
  scale_color_manual(values = color_map) +  # Custom colors
  labs(title = "", x = "", y = "Tajimas D")+
  
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

# Plot tajD ----

summary_tajDS_data$Sex.Bias.Category <- factor(summary_tajDS_data$Sex.Bias.Category, 
                                               levels = c("Highly female-biased", "Moderately female-biased", 
                                                          "Unbiased", "Moderately male-biased", "Highly male-biased"))

ggplot(summary_tajDS_data, aes(x = Sex.Bias.Category, y = Mean_tajD.S, color = Sex.Bias.Category)) +
  geom_point(size = 4) +  # Plot mean DoS
  geom_errorbar(aes(ymin = Mean_tajD.S - SE_tajD.S, ymax = Mean_tajD.S + SE_tajD.S), width = 0.2) +  # Error bars
  scale_color_manual(values = color_map) +  # Custom colors
  labs(title = "", x = "", y = "Tajimas D")+
  
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

# Plot tajD tissue

summary_tajD_tissue_data$Sex.Bias.Category <- factor(summary_tajD_tissue_data$Sex.Bias.Category, 
                                            levels = c("Highly female-biased", "Moderately female-biased", 
                                                       "Unbiased", "Moderately male-biased", "Highly male-biased"))

tajD_tissue_plot<- ggplot(summary_tajD_tissue_data, aes(x = Sex.Bias.Category, y =  Mean_tajD.S, color = Tissue)) +
  geom_point(position = position_dodge(width = 0.6),
             size = 3) +  # Plot mean DoS without dodge
  geom_errorbar(aes(ymin =  Mean_tajD.S - SE_tajD.S, ymax = Mean_tajD.S+ SE_tajD.S), width = 0.2, position = position_dodge(width = 0.6)) +  # Error bars
  scale_color_manual(values = tissue_colors) +  # Colors for tissues
  theme_bw() +
  labs(title = "", x = "", y = "Tajimas D", color = "Tissue Type") +
  geom_abline(intercept = 0, slope = 0,  size = 0.5, linetype= "solid", color = "black") +
  theme_classic() +
  theme(plot.title.position = c("panel"),
        legend.title = element_blank(),
        legend.position = c("right"),
        legend.box.background = element_rect(),
        legend.text = element_text(size = 20, color = "black"),
        axis.text.x = element_text(angle = 20, size=15, hjust = 1, margin = margin(5,0,0,0), color = "black"),
        axis.text.y = element_text(size=20, margin = margin(0,5,0,0), color = "black"), 
        axis.title.x = element_text(size=20, margin = margin(10,0,0,0), color = "black"),
        axis.title.y = element_text(size=20, margin = margin(0,10,0,0), color = "black"),
        plot.title = element_text(size=20, margin = margin(0,0,0,0), color = "black"),
        plot.margin = margin(6,6,6,6)
  )+
  
  # Add bins for each category
  geom_vline(xintercept = 1.5, linetype = "solid", color = "black") +
  geom_vline(xintercept = 2.5, linetype = "solid", color = "black") +
  geom_vline(xintercept = 3.5, linetype = "solid", color = "black") +
  geom_vline(xintercept = 4.5, linetype = "solid", color = "black")

print(tajD_tissue_plot)


# Chromosomal distribution ----

# Convert Chr to factor with an order: X, 2, then 3
chr_counts$Chr <- factor(chr_counts$Chr, levels = c("X", "2", "3"))

# Create the bar plot 
ggplot(chr_counts, aes(x = Chr, y = count, fill = Chr)) +
  geom_bar(stat = "identity") +
  labs(title = "Distribution of Chromosomes in Disease Genes",
       x = "Chromosome",
       y = "Number of Disease-Associated Genes") +
  theme_minimal() + 
  theme(panel.grid.major = element_blank(),  # Remove major gridlines
        panel.grid.minor = element_blank())  # Remove minor gridlines

# rMF ----

# rMF data 

SBGE_rMF<- mish_singh_SBG%>%
  select(FlyBaseID,Sex.Bias.Category,exp_SBGE_ase,DoS,Tissue)%>% 
  distinct() # Remove disease variable to avoid repeats

# Change into numerical 
str(SBGE_rMF$exp_SBGE_ase)
SBGE_rMF<- SBGE_rMF %>%
  mutate(exp_SBGE_ase = as.numeric(exp_SBGE_ase))

# Creating Bins
SBGE_rMF_bins<-SBGE_rMF %>%
  mutate(Sex.Bias.Category = case_when(
    exp_SBGE_ase >= 2 ~ "Highly male-biased",
    exp_SBGE_ase >= 0.5 & exp_SBGE_ase < 2 ~ "Moderately male-biased",
    exp_SBGE_ase >= -0.5 & exp_SBGE_ase < 0.5 ~ "Unbiased",
    exp_SBGE_ase >= -2 & exp_SBGE_ase < -0.5 ~ "Moderately female-biased",
    exp_SBGE_ase < -2 ~ "Highly female-biased"
  ))

# Compute Mean and Standard Error (SE) for rMF per Sex Bias Category

summary_rMF_data <-SBGE_rMF_bins %>%
  group_by(Sex.Bias.Category) %>%
  summarise(
    Mean_rmf = mean(rmf, na.rm = TRUE),
    SE_rmf = sd(rmf, na.rm = TRUE) / sqrt(n())  # Standard Error
  )


# Creating Summary table for rMF tissue data

SBGE_rMF_tissue<- mish_singh_SBG%>%
  select(FlyBaseID,Sex.Bias.Category,exp_SBGE_ase,rmf,Tissue)%>% 
  distinct() # Remove disease variable to avoid repeats

summary_tissue_rMF_data <- SBGE_rMF_bins %>%
  group_by(Sex.Bias.Category, Tissue) %>%
  summarise(
    Mean_rmf = mean(rmf, na.rm = TRUE),
    SE_rmf = sd(rmf, na.rm = TRUE) / sqrt(n()),  # Standard Error
    .groups = "drop"
  )


# rMF no tissue plot

summary_rMF_data$Sex.Bias.Category <- factor(summary_rMF_data$Sex.Bias.Category, 
                                                    levels = c("Highly female-biased", "Moderately female-biased", 
                                                               "Unbiased", "Moderately male-biased", "Highly male-biased"))

ggplot(summary_rMF_data, aes(x = Sex.Bias.Category, y = Mean_rmf, color = Sex.Bias.Category)) +
  geom_point(size = 4) +  # Plot mean DoS
  geom_errorbar(aes(ymin = Mean_rmf - SE_rmf, ymax = Mean_rmf + SE_rmf), width = 0.2) +  # Error bars
  scale_color_manual(values = color_map) +  # Custom colors
  labs(title = "", x = "", y = "rmf")+
  
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

summary_tissue_rMF_data$Sex.Bias.Category <- factor(summary_tissue_rMF_data$Sex.Bias.Category, 
                                                levels = c("Highly female-biased", "Moderately female-biased", 
                                                           "Unbiased", "Moderately male-biased", "Highly male-biased"))

# rMF tissue plot

rMF_tissue_plot<- ggplot(summary_tissue_rMF_data, aes(x = Sex.Bias.Category, y =  Mean_rmf, color = Tissue)) +
  geom_point(position = position_dodge(width = 0.6),
             size = 3) +  # Plot mean DoS without dodge
  geom_errorbar(aes(ymin =  Mean_rmf - SE_rmf, ymax = Mean_rmf + SE_rmf), width = 0.2, position = position_dodge(width = 0.6)) +  # Error bars
  scale_color_manual(values = tissue_colors) +  # Colors for tissues
  theme_bw() +
  labs(title = "", x = "", y = "rMF", color = "Tissue Type") +
  geom_abline(intercept = 0, slope = 0,  size = 0.5, linetype= "solid", color = "black") +
  theme_classic() +
  theme(plot.title.position = c("panel"),
        legend.title = element_blank(),
        legend.position = c("right"),
        legend.box.background = element_rect(),
        legend.text = element_text(size = 20, color = "black"),
        axis.text.x = element_text(angle = 45, size=15, hjust = 1, margin = margin(5,0,0,0), color = "black"),
        axis.text.y = element_text(size=20, margin = margin(0,5,0,0), color = "black"), 
        axis.title.x = element_text(size=40, margin = margin(10,0,0,0), color = "black"),
        axis.title.y = element_text(size=40, margin = margin(0,10,0,0), color = "black"),
        plot.title = element_text(size=40, margin = margin(0,0,0,0), color = "black"),
        plot.margin = margin(6,6,6,6)
  )+
  
  # Add bins for each category
  geom_vline(xintercept = 1.5, linetype = "solid", color = "black") +
  geom_vline(xintercept = 2.5, linetype = "solid", color = "black") +
  geom_vline(xintercept = 3.5, linetype = "solid", color = "black") +
  geom_vline(xintercept = 4.5, linetype = "solid", color = "black")

print(rMF_tissue_plot)

#COMPARATIVE PLOTS WITH CONTROL GENEs-----

# Tajimas D ----

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

# Compute Mean and Standard Error (SE) for tajDN and tajDS per Sex Bias Category

summary_tajDS_data_comparative <-SBGE_tajDS_bins_compare %>%
  group_by(Sex.Bias.Category,GeneType) %>%
  summarise(
    Mean_tajD.S = mean(tajD.S, na.rm = TRUE),
    SE_tajD.S = sd(tajD.S, na.rm = TRUE) / sqrt(n())  # Standard Error
  )

# summary_tissue_tajDS_data_comparative <-SBGE_tajDS_bins_compare%>%
#   group_by(Sex.Bias.Category, Tissue, GeneType) %>%
#   summarise( 
#     Mean_tajD.S = mean(tajD.S, na.rm = TRUE),
#     SE_tajD.S = sd(tajD.S, na.rm = TRUE) / sqrt(n()),  # Standard Error
#     .groups = "drop"
#   )


# # Plot using a single x-axis with grouped categories
# ggplot(summary_tajDS_data_comparative, aes(x = Sex.Bias.Category, y = Mean_tajD.S, color = GeneType)) +
#   geom_pointrange(aes(ymin = Mean_tajD.S - SE_tajD.S, ymax = Mean_tajD.S + SE_tajD.S),
#                   position = position_dodge(width = 0.5),
#                   size = 1) +  # Adjust error bar size
#   geom_point(position = position_dodge(width = 0.5), size = 3) +  # Set size of points
#   scale_color_manual(values = c("Control" = "blue", "Disease" = "red")) +  # Color points by Gene Type
#   labs(x = "Sex Bias Category", y = "Tajima's D") +
#   theme_minimal(base_size = 14) +  # Remove default gridlines
#   theme(
#     panel.grid.major = element_blank(),  # Remove major gridlines
#     panel.grid.minor = element_blank(),  # Remove minor gridlines
#     panel.border = element_rect(color = "black", fill = NA, size = 1),  # Add border around plot
#     axis.text.x = element_text(angle = 45, hjust = 1),
#     panel.background = element_blank()
#   ) +
#   geom_vline(xintercept = seq(0.5, length(unique(summary_data$Sex.Bias.Category)) + 0.5, 1),
#              linetype = "solid", color = "black", size = 0.8)  # Add bin separators


#summary_tissue_tajDS_data_comparative$Sex.Bias.Category <- factor(summary_tissue_tajDS_data_comparative$Sex.Bias.Category,
                                                     #levels = c("Highly female-biased", "Moderately female-biased",
                                                                #"Unbiased", "Moderately male-biased", "Highly male-biased"))

# 
# ggplot(summary_tissue_tajDS_data_comparative, aes(x = Sex.Bias.Category, y = Mean_tajD.S, color = GeneType)) +
#   geom_pointrange(aes(ymin = Mean_tajD.S - SE_tajD.S, ymax = Mean_tajD.S + SE_tajD.S),
#                   position = position_dodge(width = 0.5),
#                   size = 1) +  # Adjust error bar size
#   geom_point(position = position_dodge(width = 0.5), size = 3) +  # Set size of points
#   scale_color_manual(values = c("Control" = "blue", "Disease" = "red")) +  # Color points by Gene Type
#   labs(x = "Sex Bias Category", y = "Tajima's D") +
#   theme_minimal(base_size = 14) +  # Remove default gridlines
#   theme(
#     panel.grid.major = element_blank(),  # Remove major gridlines
#     panel.grid.minor = element_blank(),  # Remove minor gridlines
#     panel.border = element_rect(color = "black", fill = NA, size = 1),  # Add border around plot
#     axis.text.x = element_text(angle = 45, hjust = 1),
#     panel.background = element_blank()
#   ) +
#   geom_vline(xintercept = seq(0.5, length(unique(summary_data$Sex.Bias.Category)) + 0.5, 1),
#              linetype = "solid", color = "black", size = 0.8) +  # Add bin separators
#   facet_wrap(~Tissue, ncol = 1, )  # Facet by Tissue, one column
# 

sample_sizes_tadj <-SBGE_tajDS_bins_compare%>%
  group_by(Sex.Bias.Category, GeneType) %>%
  summarise(SampleSize = n(), .groups = "drop")  # Count observations

summary_tajDS_data_comparative <- left_join(summary_tajDS_data_comparative, sample_sizes, by = c("Sex.Bias.Category", "GeneType"))

# Then compute the "All" category stats
all_category <- summary_tajDS_data_comparative  %>%
  group_by(GeneType) %>%
  summarise(
    Mean_tajD.S = mean(Mean_tajD.S, na.rm = TRUE),  # Calculate mean across all categories
    SE_tajD.S= sqrt(sum(SE_tajD.S^2) / n()),  # Aggregate standard errors
    SampleSize = sum(SampleSize)  # Sum all sample sizes
  ) %>%
  mutate(Sex.Bias.Category = "All")


summary_tajDS_data_comparative <- bind_rows(summary_tajDS_data_comparative, all_category)

summary_tajDS_data_comparative$Sex.Bias.Category <- factor(summary_tajDS_data_comparative$Sex.Bias.Category, 
                                                     levels = c("All","Highly female-biased", "Moderately female-biased", 
                                                                "Unbiased", "Moderately male-biased", "Highly male-biased"))

# Plot using a single x-axis with grouped categories

ggplot(summary_tajDS_data_comparative, 
       aes(x = Sex.Bias.Category, y = Mean_tajD.S, color = GeneType)) + 
  geom_pointrange(aes(ymin = Mean_tajD.S - SE_tajD.S, ymax = Mean_tajD.S + SE_tajD.S),
                  position = position_dodge(width = 0.5),
                  size = 1) +  # Adjust error bar size
  geom_point(position = position_dodge(width = 0.5), size = 3) +  # Set size of points
  geom_text(aes(label = SampleSize),  
            position = position_dodge(width = 0.5),  
            color = "black", size = 3) + # Adjust label position and size
  scale_color_manual(values = c("Control" = "darkgreen", "Disease" = "purple")) +  # Match colors
  labs(x = "Sex Bias Category", y = "Tajima’s D",color = "Gene Type") +
  guides(color = guide_legend(override.aes = list(size = 1))) +
  theme_minimal(base_size = 14) +  # Remove default gridlines
  theme(
    panel.grid.major = element_blank(),  # Remove major gridlines
    panel.grid.minor = element_blank(),  # Remove minor gridlines
    panel.border = element_rect(color = "black", fill = NA, size = 1),  # Add border
    panel.background = element_blank()
  ) +
  geom_vline(xintercept = 1.5, color = "black", size = 1.2) +  # Bold divider
  geom_vline(xintercept = 2.5:5.5, color = "grey60", size = 0.5)+  # Thin grey lines
  scale_x_discrete(labels = c("All", "Highly FB", "Female-Biased", "Unbiased", "Male-Biased", "Highly MB")) 



# rMF ----

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


# Compute Mean and Standard Error (SE) for tajDN and tajDS per Sex Bias Category

summary_rMF_data_comparative <-SBGE_rMF_bins_compare%>%
  group_by(Sex.Bias.Category,GeneType) %>%
  summarise(
    Mean_rmf = mean(rmf, na.rm = TRUE),
    SE_rmf = sd(rmf, na.rm = TRUE) / sqrt(n())  # Standard Error
  )

# summary_tissue_rMF_data_comparative <-SBGE_rMF_bins_compare%>%
#   group_by(Sex.Bias.Category, Tissue, GeneType) %>%
#   summarise( 
#     Mean_rmf = mean(rmf, na.rm = TRUE),
#     SE_rmf = sd(rmf, na.rm = TRUE) / sqrt(n()),  # Standard Error
#     .groups = "drop"
#   )

rMF_sample_sizes <-SBGE_rMF_bins_compare%>%
  group_by(Sex.Bias.Category, GeneType) %>%
  summarise(SampleSize = n(), .groups = "drop")  # Count observations

summary_rMF_data_comparative <- left_join(summary_rMF_data_comparative, rMF_sample_sizes, by = c("Sex.Bias.Category", "GeneType"))

# Then compute the "All" category stats
all_category <- summary_rMF_data_comparative  %>%
  group_by(GeneType) %>%
  summarise(
    Mean_rmf = mean(Mean_rmf, na.rm = TRUE),  # Calculate mean across all categories
    SE_rmf= sqrt(sum(SE_rmf^2) / n()),  # Aggregate standard errors
    SampleSize = sum(SampleSize)  # Sum all sample sizes
  ) %>%
  mutate(Sex.Bias.Category = "All")


summary_rMF_data_comparative <- bind_rows(summary_rMF_data_comparative,all_category_rmf)

summary_rMF_data_comparative$Sex.Bias.Category <- factor(summary_rMF_data_comparative$Sex.Bias.Category, 
                                                           levels = c("All","Highly female-biased", "Moderately female-biased", 
                                                                      "Unbiased", "Moderately male-biased", "Highly male-biased"))
# rmf Plot

ggplot(summary_rMF_data_comparative, 
       aes(x = Sex.Bias.Category, y = Mean_rmf, color = GeneType)) + 
  geom_pointrange(aes(ymin = Mean_rmf - SE_rmf, ymax = Mean_rmf + SE_rmf),
                  position = position_dodge(width = 0.5),
                  size = 1) +  # Adjust error bar size
  geom_point(position = position_dodge(width = 0.5), size = 3) +  # Set size of points
  geom_text(aes(label = SampleSize),  
            position = position_dodge(width = 0.5),
            color = "black", size = 3) + # Adjust label position and size
  scale_color_manual(values = c("Control" = "darkgreen", "Disease" = "purple")) +  # Match colors
  labs(x = "Sex Bias Category", y = "rMF",color = "Gene Type") +
  guides(color = guide_legend(override.aes = list(size = 1))) +
  theme_minimal(base_size = 14) +  # Remove default gridlines
  theme(
    panel.grid.major = element_blank(),  # Remove major gridlines
    panel.grid.minor = element_blank(),  # Remove minor gridlines
    panel.border = element_rect(color = "black", fill = NA, size = 1),  # Add border
    panel.background = element_blank()
  )+
  geom_vline(xintercept = 1.5, color = "black", size = 1.2) +  # Bold divider
  geom_vline(xintercept = 2.5:5.5, color = "grey60", size = 0.5)+  # Thin grey lines
  scale_x_discrete(labels = c("All", "Highly FB", "Female-Biased", "Unbiased", "Male-Biased", "Highly MB")) 


# DoS ----
# Change into numerical 

SBGE_DoS_compare<- control_set_genes_mish_Fra_SBGE%>%
  select(FlyBaseID,Sex.Bias.Category,exp_SBGE_ase,DoS,GeneType,disease)%>% 
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

# Compute Mean and Standard Error (SE) for tajDN and tajDS per Sex Bias Category

summary_DoS_data_comparative <-SBGE_DoS_bins_compare%>%
  group_by(Sex.Bias.Category,GeneType) %>%
  summarise(
    Mean_DoS = mean(DoS, na.rm = TRUE),
    SE_DoS = sd(DoS, na.rm = TRUE) / sqrt(n())  # Standard Error
  )

DoS_sample_sizes <-SBGE_DoS_bins_compare%>%
  group_by(Sex.Bias.Category, GeneType) %>%
  summarise(SampleSize = n(), .groups = "drop")  # Count observations

summary_DoS_data_comparative <- left_join(summary_DoS_data_comparative, DoS_sample_sizes, by = c("Sex.Bias.Category", "GeneType"))

unique_DoS <- SBGE_DoS_compare %>%
  group_by(FlyBaseID, GeneType) %>%
  summarise(
    DoS = mean(DoS, na.rm = TRUE),# Or just pick one value if DoS is identical for repeats
    .groups = "drop"
  )

# Then compute the "All" category stats
DoS_all_category <- unique_DoS %>%
  group_by(GeneType) %>%
  summarise(
    Mean_DoS = mean(DoS, na.rm = TRUE),
    SE_DoS = sd(DoS, na.rm = TRUE) / sqrt(n()),  # Standard Error = SD / sqrt(n)
    SampleSize = n(),  # Unique gene count
    .groups = "drop"
  ) %>%
  mutate(Sex.Bias.Category = "All")


summary_DoS_data_comparative <- bind_rows(summary_DoS_data_comparative,DoS_all_category2)

summary_DoS_data_comparative$Sex.Bias.Category <- factor(summary_DoS_data_comparative$Sex.Bias.Category, 
                                                         levels = c("All","Highly female-biased", "Moderately female-biased", 
                                                                    "Unbiased", "Moderately male-biased", "Highly male-biased"))
# DoS Plot

ggplot(summary_DoS_data_comparative, 
       aes(x = Sex.Bias.Category, y = Mean_DoS, color = GeneType)) + 
  geom_pointrange(aes(ymin = Mean_DoS - SE_DoS, ymax = Mean_DoS + SE_DoS),
                  position = position_dodge(width = 0.5),
                  size = 1) +  # Adjust error bar size
  geom_point(position = position_dodge(width = 0.5), size = 3) +  # Set size of points
  geom_text(aes(label = SampleSize),  
            position = position_dodge(width = 0.5),  
            color = "black", size = 3) + # Adjust label position and size
  scale_color_manual(values = c("Control" = "darkgreen", "Disease" = "purple")) +  # Match colors
  labs(x = "Sex Bias Category", y = "DoS",color = "Gene Type") +
  guides(color = guide_legend(override.aes = list(size = 1))) +
  theme_minimal(base_size = 14) +  # Remove default gridlines
  theme(
    panel.grid.major = element_blank(),  # Remove major gridlines
    panel.grid.minor = element_blank(),  # Remove minor gridlines
    panel.border = element_rect(color = "black", fill = NA, size = 1),  # Add border
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.background = element_blank()
  ) +
  geom_vline(xintercept = 1.5, color = "black", size = 1.2) +  # Bold divider
  geom_vline(xintercept = 2.5:5.5, color = "grey60", size = 0.5)  # Thin grey lines







# STATISTICAL ANALYSES----

# Permutation tests----

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
TwoPerm <- function(perm_dat, x_col, groupBy, n_perm = 10000, alternative = c("two.sided","less","greater")){
  # separate data for groups 1 and 2
  trt1 <- perm_dat[perm_dat[[groupBy]] == "Control",][[x_col]]
  trt2 <- perm_dat[perm_dat[[groupBy]] == "Disease",][[x_col]]
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


# Apply the TwoPerm function to the provided dataframe
result_tajD <- TwoPerm(perm_dat = summary_tajDS_data_comparative, x_col = "Mean_tajD.S", groupBy = "GeneType", n_perm = 10000, alternative = "two.sided")

paired_tajD <-PairedTwoPerm(perm_dat = summary_tajDS_data_comparative, x_col = "Mean_tajD.S", groupBy = "GeneType", n_perm = 10000, alternative = "two.sided")

# Apply the TwoPerm_SBGE function to the provided dataframe
result_tajD_SBGE_cat <- TwoPerm_SBGE(perm_dat = summary_tajDS_data_comparative, x_col = "Mean_tajD.S", groupBy = "GeneType", SBGE_cat = "Sex.Bias.Category")

# Print the result
print(result_tajD)
print(result_tajD_SBGE_cat )

view(result_tajD_SBGE_cat )


# Apply the TwoPerm function to the provided dataframe
result_DoS <- TwoPerm(perm_dat = summary_DoS_data_comparative, x_col = "Mean_DoS", groupBy = "GeneType", n_perm = 10000, alternative = "two.sided")

# Apply the TwoPerm_SBGE function to the provided dataframe
result_DoS_SBGE_cat <- TwoPerm_SBGE(perm_dat = summary_DoS_data_comparative, x_col = "Mean_DoS", groupBy = "GeneType", SBGE_cat = "Sex.Bias.Category")

# Print the result
print(result_DoS)
print(result_DoS_SBGE_cat )

view(result_DoS_SBGE_cat )


# Apply the TwoPerm function to the provided dataframe
result_rmf <- TwoPerm(perm_dat = summary_rMF_data_comparative, x_col = "Mean_rmf", groupBy = "GeneType", n_perm = 10000, alternative = "two.sided")

# Apply the TwoPerm_SBGE function to the provided dataframe
result_rmf_SBGE_cat <- TwoPerm_SBGE(perm_dat = summary_rMF_data_comparative, x_col = "Mean_rmf", groupBy = "GeneType", SBGE_cat = "Sex.Bias.Category")

# Print the result
print(result_rmf)
print(result_rmf_SBGE_cat )

view(result_rmf_SBGE_cat )
