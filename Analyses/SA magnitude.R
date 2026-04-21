#extract genes from final_leading_snps

library(biomaRt)
library(readxl)
library(data.table)




snps <- read.delim("/Users/elliehowatson/Final_Year_Project/R/data/true_cmh_results.txt")
colnames(snps) <- c(
  "Chromosome",
  "Basepair",
  "Reference",
  paste0("Sample_", 1:28),
  "Log10_pvalue",
  "LogOdds"
)


snps$Chromosome <- as.character(snps$Chromosome)
snps$Basepair   <- as.numeric(snps$Basepair)

snps$region <- paste0(
  snps$Chromosome, ":",
  snps$Basepair, ":",
  snps$Basepair
)

head(snps$region)

ensembl <- useMart("ensembl", dataset = "dmelanogaster_gene_ensembl")

#extract genes from list
gene_hits_snps <- getBM(
  attributes = c("chromosome_name", "start_position", "end_position",
                 "ensembl_gene_id", "external_gene_name", "gene_biotype"),
  filters = "chromosomal_region",
  values = snps$region,
  mart = ensembl
)

# Convert to data.table
snps_dt <- as.data.table(snps)
genes_dt <- as.data.table(gene_hits_snps)

# Set keys for faster joining
setkey(snps_dt, Chromosome, Basepair)
setkey(genes_dt, chromosome_name, start_position, end_position)

#Find which genes are associated to each SNP
gene_snp_associations <- genes_dt[snps_dt,
                                  on = .(chromosome_name = Chromosome,
                                         start_position <= Basepair,
                                         end_position >= Basepair),
                                  nomatch = 0,
                                  allow.cartesian = TRUE
]

# Calculate mean scores for each gene
gene_scores <- gene_snp_associations[,
                                     .(mean_Log10_pvalue = mean(Log10_pvalue, na.rm = TRUE),
                                       mean_LogOdds = mean(LogOdds, na.rm = TRUE),
                                       n_snps = .N),
                                     by = .(external_gene_name, ensembl_gene_id, gene_biotype)
][!is.na(external_gene_name)][order(-mean_Log10_pvalue)]


# View results
head(gene_scores, 20)



write.csv(gene_scores, "/Users/elliehowatson/Final_Year_Project/R/data/gene_SA_magnitude.csv", row.names = FALSE)




#merge with snps data to get SA magnitude
gene_hits_snps$Chromosome <- gene_hits_snps$chromosome_name
gene_hits_snps$Basepair <- gene_hits_snps$start_position




#Merge gene names with orthlog, disease, background and candidate genes

# 1) Make a dataframe with SBGE and all orthlogs
SBGE_all_orthologs_data <- SBG_whole_bodies %>%
  inner_join(all_orthologs_list, by = "FlyBaseID")

SBGE_all_orthologs_data$GeneType<- "Orthologs"

# 2) Remove orthologues that are also in the candidate disease set
SBGE_all_orthologs_data  <- SBGE_all_orthologs_data [ !(SBGE_all_orthologs_data $FlyBaseID %in% SBGE_disease_data$FlyBaseID), ]

# 3) Bind together
SBGE_disease_data <-SBGE_disease_data%>%
  rename(
    log2FoldChange = exp_SBGE_ase
  )
SBGE_disease_data<- SBGE_disease_data[, c("FlyBaseID", "Sex.Bias.Category", "log2FoldChange", "GeneType")]
candidate_ortholog_SBGE<- rbind(SBGE_disease_data, SBGE_all_orthologs_data)
candidate_ortholog_SBGE<- candidate_ortholog_SBGE[!duplicated(candidate_ortholog_SBGE), ]

#Do the same for control genes

SBGE_control_data  <- control_mish_SBGE [ !(control_mish_SBGE $FlyBaseID %in% SBGE_disease_data$FlyBaseID), ]

#Bind together control SBGE and candidate

SBGE_control_data <-SBGE_control_data%>%
  rename(
    log2FoldChange = exp_SBGE_ase
  )
SBGE_control_data<- SBGE_control_data[, c("FlyBaseID", "Sex.Bias.Category", "log2FoldChange", "GeneType")]

candidate_control_SBGE<- rbind(SBGE_disease_data, SBGE_control_data)

#Do the same for all disease genes

# 1) Make a dataframe with SBGE and all orthlogs
SBGE_all_disease_data <- SBG_whole_bodies %>%
  inner_join(all_disease_genes_list, by = "FlyBaseID")

SBGE_all_disease_data$GeneType<- "Disease"

# 2) Remove disease genes that are also in the candidate disease set
SBGE_all_disease_data <- SBGE_all_disease_data [ !(SBGE_all_disease_data $FlyBaseID %in% SBGE_disease_data$FlyBaseID), ]

# 3) Bind together

candidate_all_disease_SBGE<- rbind(SBGE_disease_data, SBGE_all_disease_data )
candidate_all_disease_SBGE<- candidate_all_disease_SBGE[!duplicated(candidate_all_disease_SBGE), ]


# Merge each gene set with the snps dataset

gene_hits_snps <-gene_hits_snps%>%
  rename(
    FlyBaseID= ensembl_gene_id
  )

#ortholog
ortholog_snps <- left_join(gene_hits_snps, candidate_ortholog_SBGE, by = "FlyBaseID")

#all disease
all_disease_snps <- left_join(gene_hits_snps, candidate_all_disease_SBGE, by = "FlyBaseID")

#background
background_snps <- left_join(gene_hits_snps, candidate_control_SBGE, by = "FlyBaseID")



# Ortholog SA magnitude Summary table ----

ortholog_snps_compare<- ortholog_snps %>%
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


write.csv(gene_hits_snps, file = "/Users/elliehowatson/Final_Year_Project/R/data/gene_hits_snps.csv", row.names = FALSE) #save as csv
