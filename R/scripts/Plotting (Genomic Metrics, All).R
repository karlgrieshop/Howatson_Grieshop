#PLOTTING
#_____________________----

library(ggplot2)
library(patchwork)

pub_theme <- function(base_size = 14) {
  theme_classic(base_size = base_size) +
    theme(
      axis.text        = element_text(size = base_size, colour = "black"),
      axis.title       = element_text(size = base_size + 2, face = "bold"),
      legend.text      = element_text(size = base_size),
      legend.title     = element_text(size = base_size + 1, face = "bold"),
      legend.position  = "top",
      panel.border     = element_rect(colour = "black", fill = NA, linewidth = 0.8),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      plot.margin      = margin(12, 12, 12, 12)
    )
}


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
                  linewidth = 1.2) +
  geom_text(aes(y = q95 + 0.006, label = SampleSize, group = GeneType), 
            position = position_dodge(width = 0.5),
            color = "black", size = 6) +
  scale_color_manual(values = c("SSDOs" = "#D55E00", "All Disease Orthologues" = "#0072B2",
                                "All Orthologues" = "#009E73", "Background" = "#999999")) +
  labs(x = "Gene Type", 
       y = expression(tilde(pi)[italic(r)*N])) +
  pub_theme() +
  theme(legend.position = "none")
print(Fig1A)


# ---- DoS ----

DoS_plot_data_clean $GeneType <- factor(DoS_plot_data_clean $GeneType, 
                                        levels = c("SSDOs","All Disease Orthologues", "All Orthologues", 
                                                   "Background"))

Fig1B<- ggplot(DoS_plot_data_clean , 
               aes(x = GeneType, y = q50, color = GeneType)) + 
  geom_pointrange(aes(ymin = q05, ymax = q95), 
                  position = position_dodge(width = 0.5), 
                  linewidth = 1.2) +
  geom_text(aes(y = q95 + 0.006, label = SampleSize, group = GeneType), 
            position = position_dodge(width = 0.5),
            color = "black", size = 6) +
  scale_color_manual(values = c("SSDOs" = "#D55E00", "All Disease Orthologues" = "#0072B2",
                                "All Orthologues" = "#009E73", "Background" = "#999999")) +
  labs(x = "Gene Type", 
       y = expression(italic(DoS))) +
  pub_theme() +
  theme(legend.position = "none")

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
                  linewidth = 1.2) +
  geom_text(aes(y = q95 + 0.006, label = SampleSize, group = GeneType), 
            position = position_dodge(width = 0.5),
            color = "black", size = 6) +
  scale_color_manual(values = c("SSDOs" = "#D55E00", "All Disease Orthologues" = "#0072B2",
                                "All Orthologues" = "#009E73", "Background" = "#999999")) +
  labs(x = "Gene Type", 
       y = expression(italic(r)[MF])) +
  pub_theme() +
  theme(legend.position = "none")
print(Fig2A)

# ---- Tajimas D ----


tajD_plot_data_clean $GeneType <- factor(tajD_plot_data_clean $GeneType, 
                                         levels = c("SSDOs","All Disease Orthologues", "All Orthologues", 
                                                    "Background"))


Fig2B<- ggplot(tajD_plot_data_clean , 
                       aes(x = GeneType, y = q50, color = GeneType)) + 
  geom_pointrange(aes(ymin = q05, ymax = q95), 
                  position = position_dodge(width = 0.5), 
                  linewidth = 1.2) +
  geom_text(aes(y = q95 + 0.006, label = SampleSize, group = GeneType), 
            position = position_dodge(width = 0.5),
            color = "black", size = 6) +
  scale_color_manual(values = c("SSDOs" = "#D55E00", "All Disease Orthologues" = "#0072B2",
                                "All Orthologues" = "#009E73", "Background" = "#999999")) +
  labs(x = "Gene Type", 
       y = expression("Tajima's " * italic(D))) +
  pub_theme() +
  theme(legend.position = "none")
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
                  linewidth = 1.2) +
  geom_text(aes(y = q95 + 0.006, label = SampleSize, group = GeneType), 
            position = position_dodge(width = 0.5),
            color = "black", size = 6) +
  scale_color_manual(values = c("SSDOs" = "#D55E00", "All Disease Orthologues" = "#0072B2",
                                "All Orthologues" = "#009E73", "Background" = "#999999")) +
  labs(x = "Gene Type", 
       y = "Tissue Specificity") +
  pub_theme() +
  theme(legend.position = "none")
print(Fig2C)


