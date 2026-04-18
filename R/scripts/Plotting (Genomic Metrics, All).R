#PLOTTING
#_____________________----

library(ggplot2)
library(patchwork)


#Upload Summary Plotting Data----


tajD_plot_data_clean <-read.csv("/Users/elliehowatson/Final_Year_Project/R/data/Summary Tables/tajD Plot Data.csv")
rmf_plot_data_clean  <-read.csv("/Users/elliehowatson/Final_Year_Project/R/data/Summary Tables/rMF Plot Data.csv")
rlnd_plot_data_clean  <-read.csv("/Users/elliehowatson/Final_Year_Project/R/data/Summary Tables/rlnd Plot Data.csv")
tissuespec_plot_data_clean <-read.csv("/Users/elliehowatson/Final_Year_Project/R/data/Summary Tables/Tissue Specificity Plot Data.csv")
DoS_plot_data_clean<-read.csv("/Users/elliehowatson/Final_Year_Project/R/data/Summary Tables/DoS Plot Data.csv")


#PLOTS----
#ALL 4 POINTS ON ONE PLOT----
#  ---- rlnd ----

rlnd_plot_data_clean $GeneType <- factor(rlnd_plot_data_clean $GeneType, 
                                         levels = c("SSDOs","All Disease Orthologues", "All Orthologues", 
                                                    "Background"))

Fig1A<- ggplot(rlnd_plot_data_clean , 
               aes(x = GeneType, y = q50, color = GeneType)) + 
  geom_pointrange(aes(ymin = q05, ymax = q95), 
                  position = position_dodge(width = 0.5), 
                  size = 1.2) +
  geom_point(aes(color = GeneType),
             position = position_dodge(width = 0.5),
             size = 2) +  
  geom_errorbar(aes(ymin = q05, ymax = q95), 
                position = position_dodge(width = 0.5),
                width = 0.4, 
                size = 1) +
  geom_text(aes(y = q95 + 0.006, label = SampleSize, group = GeneType), 
            position = position_dodge(width = 0.5),
            color = "black", size = 6) +
  scale_color_manual(values = c("All Disease Orthologues" = "darkgreen", "SSDOs" = "purple", 
                                "All Orthologues" = "orange", "Background" = "blue")) +
  labs(x = "Gene Type", 
       y = expression(tilde(pi)[italic(r)*N])) +
  theme_classic(base_size = 16) +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 14),
    axis.text.y = element_text(size = 14),
    axis.title = element_text(size = 16, face = "bold"),
    legend.position = "none",
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 1.2),
    panel.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(20, 20, 20, 20)
  )
print(Fig1A)


# ---- DoS ----

DoS_plot_data_clean $GeneType <- factor(DoS_plot_data_clean $GeneType, 
                                        levels = c("SSDOs","All Disease Orthologues", "All Orthologues", 
                                                   "Background"))

Fig1B<- ggplot(DoS_plot_data_clean , 
               aes(x = GeneType, y = q50, color = GeneType)) + 
  geom_pointrange(aes(ymin = q05, ymax = q95), 
                  position = position_dodge(width = 0.5), 
                  size = 1.2) +
  geom_point(aes(color = GeneType),
             position = position_dodge(width = 0.5),
             size = 2) +  
  geom_errorbar(aes(ymin = q05, ymax = q95), 
                position = position_dodge(width = 0.5),
                width = 0.4, 
                size = 1) +
  geom_text(aes(y = q95 + 0.006, label = SampleSize, group = GeneType), 
            position = position_dodge(width = 0.5),
            color = "black", size = 6) +
  scale_color_manual(values = c("All Disease Orthologues" = "darkgreen", "SSDOs" = "purple", 
                                "All Orthologues" = "orange", "Background" = "blue")) +
  labs(x = "Gene Type", 
       y = expression(italic(DoS))) +
  theme_classic(base_size = 16) +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 14),
    axis.text.y = element_text(size = 14),
    axis.title = element_text(size = 16, face = "bold"),
    legend.position = "none",
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 1.2),
    panel.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(20, 20, 20, 20)
  )

print(Fig1B)
# --- rMF ----

rmf_plot_data_clean $GeneType <- factor(rmf_plot_data_clean $GeneType, 
                                        levels = c("SSDOs","All Disease Orthologues", "All Orthologues", 
                                                   "Background"))
#rmf plot----

Fig2A<- ggplot(rmf_plot_data_clean , 
               aes(x = GeneType, y = q50, color = GeneType)) + 
  geom_pointrange(aes(ymin = q05, ymax = q95), 
                  position = position_dodge(width = 0.5), 
                  size = 1.2) +
  geom_point(aes(color = GeneType),
             position = position_dodge(width = 0.5),
             size = 2) +  
  geom_errorbar(aes(ymin = q05, ymax = q95), 
                position = position_dodge(width = 0.5),
                width = 0.4, 
                size = 1) +
  geom_text(aes(y = q95 + 0.006, label = SampleSize, group = GeneType), 
            position = position_dodge(width = 0.5),
            color = "black", size = 6) +
  scale_color_manual(values = c("All Disease Orthologues" = "darkgreen", "SSDOs" = "purple", 
                                "All Orthologues" = "orange", "Background" = "blue")) +
  labs(x = "Gene Type", 
       y = expression(italic(r)[MF])) +
  theme_classic(base_size = 16) +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 14),
    axis.text.y = element_text(size = 14),
    axis.title = element_text(size = 16, face = "bold"),
    legend.position = "none",
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 1.2),
    panel.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(20, 20, 20, 20)
  )
print(Fig2A)

# ---- Tajimas D ----


tajD_plot_data_clean $GeneType <- factor(tajD_plot_data_clean $GeneType, 
                                         levels = c("SSDOs","All Disease Orthologues", "All Orthologues", 
                                                    "Background"))


Fig2B<- ggplot(tajD_plot_data_clean , 
                       aes(x = GeneType, y = q50, color = GeneType)) + 
  geom_pointrange(aes(ymin = q05, ymax = q95), 
                  position = position_dodge(width = 0.5), 
                  size = 1.2) +
  geom_point(aes(color = GeneType),
             position = position_dodge(width = 0.5),
             size = 2) +  
  geom_errorbar(aes(ymin = q05, ymax = q95), 
                position = position_dodge(width = 0.5),
                width = 0.4, 
                size = 1) +
  geom_text(aes(y = q95 + 0.006, label = SampleSize, group = GeneType), 
            position = position_dodge(width = 0.5),
            color = "black", size = 6) +
  scale_color_manual(values = c("All Disease Orthologues" = "darkgreen", "SSDOs" = "purple", 
                                "All Orthologues" = "orange", "Background" = "blue")) +
  labs(x = "Gene Type", 
       y = expression("Tajima's " * italic(D))) +
  theme_classic(base_size = 16) +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 14),
    axis.text.y = element_text(size = 14),
    axis.title = element_text(size = 16, face = "bold"),
    legend.position = "none",
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 1.2),
    panel.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(20, 20, 20, 20)
  )
print(Fig2B)



# TissueSpecificty ----

tissuespec_plot_data_clean $GeneType <- factor(tissuespec_plot_data_clean $GeneType, 
                                               levels = c("SSDOs","All Disease Orthologues", "All Orthologues", 
                                                          "Background"))
#tissuespec plot----

Fig2C<- ggplot(tissuespec_plot_data_clean , 
                             aes(x = GeneType, y = q50, color = GeneType)) + 
  geom_pointrange(aes(ymin = q05, ymax = q95), 
                  position = position_dodge(width = 0.5), 
                  size = 1.2) +
  geom_point(aes(color = GeneType),
             position = position_dodge(width = 0.5),
             size = 2) +  
  geom_errorbar(aes(ymin = q05, ymax = q95), 
                position = position_dodge(width = 0.5),
                width = 0.4, 
                size = 1) +
  geom_text(aes(y = q95 + 0.006, label = SampleSize, group = GeneType), 
            position = position_dodge(width = 0.5),
            color = "black", size = 6) +
  scale_color_manual(values = c("All Disease Orthologues" = "darkgreen", "SSDOs" = "purple", 
                                "All Orthologues" = "orange", "Background" = "blue")) +
  labs(x = "Gene Type", 
       y = "Tissue Specificity") +
  theme_classic(base_size = 16) +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 14),
    axis.text.y = element_text(size = 14),
    axis.title = element_text(size = 16, face = "bold"),
    legend.position = "none",
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 1.2),
    panel.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(20, 20, 20, 20)
  )
print(Fig2C)


