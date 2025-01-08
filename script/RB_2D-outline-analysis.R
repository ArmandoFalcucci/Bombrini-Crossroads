# 0. Loading necessary packages ----
library(Momocs)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(tidyr)
library(viridis)
library(cowplot)
library(vegan)
library(pairwiseAdonis)
library(MetBrewer)
library(dispRity)

zapotec_palette <- c(
  "Local" = "#D14A34",        # Bold red
  "Circum-local" = "#C09B29", # Rich blue
  "Distant" = "#3D6F9E",      # Earthy orange
  "Very distant" = "#F1A64C"# Golden yellow
)

A1_A2_palette <- c(
  "A1" = "#6F9E58",   # Calming green for A1
  "A2" = "#8A3E28"    # Deep red for A2
)

Local_NonLocal_palette <- c(
  "Local" = "#C75C4A",   # Warm terracotta red for Local (Adobe color)
  "Non-local" = "#A58D65" # Muted sandy beige for Non-local (Desert earth tone)
)


# 1. Reading file names ----
lf <- list.files("data/outlines", pattern = "\\.txt$", full.names = TRUE)

# 2. Modifying the dataset ----
Dataset_2DGM <- Dataset %>%
  filter(Blank == "Bladelet", 
         Preservation == "Complete", 
         Layer %in% c("A1", "A2"), 
         raw.material.source != "not.available" & raw.material.source != "Undetermined", 
         Class == "Blank")

# Exclude specific IDs
ids_to_exclude <- c("RB_231", "RB_446", "RB_451", "RB_504", "RB_825") 
Dataset_2DGM <- Dataset_2DGM %>%
  filter(!ID %in% ids_to_exclude)

# 3. Extracting coordinates and attaching IDs ----
Coordinates <- import_txt(lf)
Coordinates_2DGM <- Dataset_2DGM %>%
  pull("ID") %>%
  as.character()

Coordinates_2DGM_out <- Coordinates[Coordinates_2DGM]
GM_bladelets <- Out(Coordinates_2DGM_out, fac = select(Dataset_2DGM, c("Layer", "raw.material.source", "raw.material.provenience", "raw.material.lumped")))

# 4. GMM procedures to centre, scale, and rotate the coordinates ----
GM_bladelets_centered <- Momocs::coo_centre(GM_bladelets)
GM_bladelets_centered_scaled <- Momocs::coo_scale(GM_bladelets_centered)
GM_bladelets_centered_scaled <- Momocs::coo_slidedirection(GM_bladelets_centered_scaled)
GM_bladelets_centered_scaled_rotated <- Momocs::coo_rotatecenter(GM_bladelets_centered_scaled, theta = -pi/2)

# 5. Calibrating harmonic power ----
harmonic_power <- calibrate_harmonicpower_efourier(GM_bladelets_centered_scaled_rotated)

# 6. EFA analysis ----
GM_bladelets_centered_scaled_rotated.EFA <- efourier(GM_bladelets_centered_scaled_rotated, nb.h = 26)

# 7. PCA (Principal Component Analysis) ----
GM_bladelets_centered_scaled_rotated.PCA <- PCA(GM_bladelets_centered_scaled_rotated.EFA)

# 8. Screeplot ----
GM_bladelets_screeplot <- Momocs::scree_plot(GM_bladelets_centered_scaled_rotated.PCA, nax = 1:8) +
  cowplot::theme_minimal_grid() +
  theme(plot.title = element_blank(),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 12),
        text = element_text(size = 14)) +
  labs(x = "Components", y = "Proportion")

# Save screeplot
ggsave(GM_bladelets_screeplot,
       filename = "output/figures/Figure_S7.png",
       width = 4, height = 4, dpi = 300, units = "in", device = "png")

# 9. Shape variation plot ----
GM_bladelets_shape_variation <- Momocs::PCcontrib(GM_bladelets_centered_scaled_rotated.PCA, nax = 1:3, sd.r = c(-2, -1, 0, 1, 2))

# Flipped shape variation plot
GM_bladelets_shape_variation_flipped <- GM_bladelets_shape_variation$gg +
  scale_y_continuous(labels = NULL) +
  scale_x_reverse(labels = NULL, sec.axis = sec_axis(~., name = "PC")) +
  coord_flip() +
  labs(title = "Mean + SD") +
  theme(plot.title = element_text(hjust = 0.5, size = 14)) +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        text = element_text(size = 14))

# Save shape variation plot
ggsave(GM_bladelets_shape_variation_flipped,
       filename = "output/figures/Figure_S8.png",
       width = 5, height = 5, dpi = 300, units = "in", device = "png")

#10 Merging datasets to produce figures and run other analyses ####
GM_PCScores <- as.data.frame(GM_bladelets_centered_scaled_rotated.PCA$x)

GM_PCScores$ID <- Dataset_2DGM$ID

Dataset_2DGM <- left_join(Dataset_2DGM, GM_PCScores[ ,c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10", "ID")], by = "ID") 

# 11. Correlation between measurements and principal components ----

# Calculate Spearman correlation for Length and PC1
cor_length_pc1 <- cor.test(Dataset_2DGM$Length, Dataset_2DGM$PC1, method = "spearman", exact = FALSE)

# Calculate Spearman correlation for Width and PC1
cor_width_pc1 <- cor.test(Dataset_2DGM$Width, Dataset_2DGM$PC1, method = "spearman", exact = FALSE)

# Calculate Spearman correlation for Thickness and PC1
cor_thickness_pc1 <- cor.test(Dataset_2DGM$Thickness, Dataset_2DGM$PC1, method = "spearman", exact = FALSE)

# Calculate Spearman correlation for Length and PC2
cor_length_pc2 <- cor.test(Dataset_2DGM$Length, Dataset_2DGM$PC2, method = "spearman", exact = FALSE)

# Calculate Spearman correlation for Width and PC2
cor_width_pc2 <- cor.test(Dataset_2DGM$Width, Dataset_2DGM$PC2, method = "spearman", exact = FALSE)

# Calculate Spearman correlation for Thickness and PC2
cor_thickness_pc2 <- cor.test(Dataset_2DGM$Thickness, Dataset_2DGM$PC2, method = "spearman", exact = FALSE)

# Calculate Spearman correlation for Length and PC3
cor_length_pc3 <- cor.test(Dataset_2DGM$Length, Dataset_2DGM$PC3, method = "spearman", exact = FALSE)

# Calculate Spearman correlation for Width and PC3
cor_width_pc3 <- cor.test(Dataset_2DGM$Width, Dataset_2DGM$PC3, method = "spearman", exact = FALSE)

# Calculate Spearman correlation for Thickness and PC3
cor_thickness_pc3 <- cor.test(Dataset_2DGM$Thickness, Dataset_2DGM$PC3, method = "spearman", exact = FALSE)


# Create a data frame to store correlation results
correlation_results <- data.frame(
  Variable1 = c("Length", "Width", "Thickness", "Length", "Width", "Thickness", "Length", "Width", "Thickness"),
  Variable2 = rep(c("PC1", "PC2", "PC3"), each = 3),
  Spearman_Correlation = c(cor_length_pc1$estimate, cor_width_pc1$estimate, cor_thickness_pc1$estimate,
                           cor_length_pc2$estimate, cor_width_pc2$estimate, cor_thickness_pc2$estimate,
                           cor_length_pc3$estimate, cor_width_pc3$estimate, cor_thickness_pc3$estimate),
  p_value = c(cor_length_pc1$p.value, cor_width_pc1$p.value, cor_thickness_pc1$p.value,
              cor_length_pc2$p.value, cor_width_pc2$p.value, cor_thickness_pc2$p.value,
              cor_length_pc3$p.value, cor_width_pc3$p.value, cor_thickness_pc3$p.value)
)

# Optional: Round p-values for better readability
correlation_results$p_value <- round(correlation_results$p_value, 3)

# Print the table
print(correlation_results)

# Clean the column names
# correlation_results <- clean_names(correlation_results)

colnames(correlation_results) <- c("Measurement", "Principal Component", "Spearman", "p-value")

# Save the correlation results as an Rdata file
save(correlation_results, file = "data/correlation_results_table.Rdata")

# 12. Plotting results of PCA ----

Means_PC1_PC3.layer <- Dataset_2DGM %>%                                               #subset dataset by Artifact Class and calculate means for PCs for each Class
  group_by(Layer) %>% 
  dplyr::summarise(PC1 = mean(PC1),
                   PC2 = mean(PC2),
                   PC3 = mean(PC3))

Means_PC1_PC3.provenience <- Dataset_2DGM %>%                                      #subset dataset by Artifact Class and calculate means for PCs for each Class
  group_by(raw.material.provenience) %>% 
  dplyr::summarise(PC1 = mean(PC1),
                   PC2 = mean(PC2),
                   PC3 = mean(PC3))

Means_PC1_PC3.source <- Dataset_2DGM %>% 
  mutate(raw.material.lumped = fct_relevel(raw.material.lumped, "Local", "Circum.local", "Distant", "Very.distant")) %>%
  mutate(raw.material.lumped = recode(raw.material.lumped, Circum.local = "Circum-local", Very.distant = "Very distant")) %>%
  #subset dataset by Artifact Class and calculate means for PCs for each Class
  group_by(raw.material.lumped) %>% 
  dplyr::summarise(PC1 = mean(PC1),
                   PC2 = mean(PC2),
                   PC3 = mean(PC3))



PC1toPC2.layer <- ggplot(data = Dataset_2DGM, aes(x = PC1, y = PC2, color = Layer)) +
  geom_point(size = 1.5, alpha = 0.5) +
  labs(y= "PC2 (14% of total variance)", x = "PC1 (62% of total variance)") +
  geom_hline(yintercept = 0, 
             linetype = "dashed", 
             color = "black",
             size = 1, 
             alpha = 0.5) +
  geom_vline(xintercept = 0, 
             linetype = "dashed", 
             color = "black",
             size = 1, 
             alpha = 0.5)  +
  geom_point(data = Means_PC1_PC3.layer,
             aes(x = PC1, y = PC2),
             size = 4.5, shape = 16) +
  guides(colour = guide_legend(override.aes = list(size=5))) +
  theme_pubclean() +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size=14),
        axis.title.x = element_text(size = 14, hjust = -0.2),
        axis.title.y = element_text(size = 14),
        axis.text = element_text(size=12)) +
  scale_color_manual(labels=c("A1 (n = 171)", "A2 (n = 198)"), values=A1_A2_palette)

PC1toPC2.provenience <- Dataset_2DGM %>%
  mutate(raw.material.lumped = fct_relevel(raw.material.lumped, "Local", "Circum.local", "Distant", "Very.distant")) %>%
  mutate(raw.material.lumped = recode(raw.material.lumped, Circum.local = "Circum-local", Very.distant = "Very distant")) %>%
  ggplot(aes(x = PC1, y = PC2, color = raw.material.lumped)) +
  geom_point(size = 1.5, alpha = 0.5) +
  labs(y= "PC2 (14% of total variance)", x = "PC1 (62% of total variance)") +
  geom_hline(yintercept = 0, 
             linetype = "dashed", 
             color = "black",
             size = 1, 
             alpha = 0.5) +
  geom_vline(xintercept = 0, 
             linetype = "dashed", 
             color = "black",
             size = 1, 
             alpha = 0.5)  +
  geom_point(data = Means_PC1_PC3.source,
             aes(x = PC1, y = PC2),
             size = 4.5, shape = 16) +
  guides(colour = guide_legend(override.aes = list(size=5), ncol = 2)) +
  theme_pubclean() +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size=14),
        axis.title.x = element_text(size = 14, hjust = -0.2),
        axis.title.y = element_text(size = 14),
        axis.text = element_text(size=12)) +
  scale_color_manual(labels=c("Local (n = 121)", "Circum-local (n = 44)", "Distant (n = 95)", "Very distant (n = 109)"), 
                     values=zapotec_palette)


# Combine the plots into a single figure
PC1toPC2.plots.combined <- ggarrange(PC1toPC2.layer, PC1toPC2.provenience,
                                     labels = c("a", "b"),
                                     nrow = 1,
                                     widths = c(1, 1),
                                     legend = "bottom",
                                     align = "h") + 
  theme(
    text = element_text(size = 18),
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 18),
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 18),
    plot.title = element_text(size = 20)
  )


ggsave("output/figures/Figure_7.tiff", PC1toPC2.plots.combined, width = 11.7, height = 6.3, units = "in")


# 13. PERMANOVA ----
min_n_PCs.2DGM <- Momocs::scree_min(GM_bladelets_centered_scaled_rotated.PCA, prop = 0.95) 

# Create the Y matrix of variables under comparison:
Y.PERMANOVA <- Dataset_2DGM[, c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9")]

# Perform one-way PERMANOVA
GM_PERMANOVA.layer <- vegan::adonis2(Y.PERMANOVA ~ Dataset_2DGM$Layer, method = "euclidean", permutations = 10000)
GM_PERMANOVA.rawmatsource <- vegan::adonis2(Y.PERMANOVA ~ Dataset_2DGM$raw.material.lumped, method = "euclidean", permutations = 10000)
GM_PERMANOVA.rawmatprov <- vegan::adonis2(Y.PERMANOVA ~ Dataset_2DGM$raw.material.provenience, method = "euclidean", permutations = 10000)

# Pairwise differences:
Pairwise_PERMANOVA.rawmatsource <- pairwiseAdonis::pairwise.adonis(Y.PERMANOVA, Dataset_2DGM$raw.material.lumped, sim.method = "euclidean", p.adjust.m = "bonferroni", perm = 10000)

# 14. Disparity Test ----
Dataset_2DGM$Layer_RawMat <- paste(Dataset_2DGM$Layer, Dataset_2DGM$raw.material.lumped, sep = ", ")
Dataset_2DGM$Layer_RawMat <- dplyr::recode(Dataset_2DGM$Layer_RawMat, "A1, Local" = "A1, Local", "A2, Local" = "A2, Local", "A1, Circum.local" = "A1, Non-local", "A2, Circum.local" = "A2, Non-local", "A1, Distant" = "A1, Non-local", "A2, Distant" = "A2, Non-local", "A1, Very.distant" = "A1, Non-local", "A2, Very.distant" = "A2, Non-local") %>%
  fct_relevel("A1, Local", "A2, Local", "A1, Non-local", "A2, Non-local")

# Ensure Dataset_2DGM is a data frame or matrix
Dataset_2DGM <- as.data.frame(Dataset_2DGM)

# Subset the columns correctly (make sure PC1 to PC9 are the correct column names)
PC.data.subset.GM <- Dataset_2DGM[, c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9")]

rownames_DATASETS <- as.factor(Dataset_2DGM$Layer_RawMat)

# Disparity test 
data_subsets <- dispRity::custom.subsets(PC.data.subset.GM, group = rownames_DATASETS)
data_boot <- boot.matrix(data_subsets, bootstraps = 10000)
data_disp <- dispRity(data_boot, metric = c(sum, variances))

# Wilcox.test
pairwise_results.disparity.shape <- test.dispRity(data_disp, test = wilcox.test, comparisons = "pairwise", correction = "bonferroni")

# Format and display results
disparity.shape.formatted_results <- as.data.frame(pairwise_results.disparity.shape)
disparity.shape.formatted_results <- disparity.shape.formatted_results %>%
  mutate(p.value = format.pval(p.value, digits = 3))

# Disparity Plot
data_names <- names(data_disp$disparity)
disparity_df_list <- list()
for(i in data_names){
  disparity_df_list[[i]] <- data.frame(Context = paste0(i, 
                                                        "\n(n=",nrow(data_disp$subsets[[i]]$elements),")"),
                                       disparity = as.vector(data_disp$disparity[[i]][[2]]),
                                       nelements = nrow(data_disp$subsets[[i]]$elements),
                                       TS = i)
}
disparity_df_discrete_GM <- do.call(rbind.data.frame, disparity_df_list)

# Plot with the updated categories (Local vs Non-local)
disparity_df_discrete_GM.plot <- 
  disparity_df_discrete_GM %>%
  ggplot(aes(x = Context, y = disparity)) +
  geom_violin(aes(fill = TS)) + 
  geom_boxplot(notch = TRUE, width = 0.1, fill = "white", color = "black") +
  theme_bw() +
  ggtitle(NULL) +
  xlab("") + 
  ylab("Disparity") +
  theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
        axis.text = element_text(size = 16),
        axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.5, size = 16),
        axis.title.y = element_text(vjust = 0),
        axis.title = element_text(size = 18)) +
  scale_fill_manual(values = c("A1, Local" = "#C75C4A",  # Second Zapotec color for Local
                               "A2, Local" = "#C75C4A",  # Second Zapotec color for Local
                               "A1, Non-local" = "#A58D65",  # Third Zapotec color for Non-local
                               "A2, Non-local" = "#A58D65"   # Third Zapotec color for Non-local
  )) +
  guides(color = FALSE, fill = FALSE)


# Display the plot
disparity_df_discrete_GM.plot

# Save the disparity test plot
ggsave("output/figures/Figure_8.tiff", disparity_df_discrete_GM.plot, width = 6, height = 4.5, units = "in")
