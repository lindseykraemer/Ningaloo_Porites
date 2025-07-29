##################################### Diff Abundance (DHW) ##############################

# VST and differential abundance testing using actual formula

dds2 <- phyloseq_to_deseq2(new_phylo, ~ DHW_8)
dds2


counts2 <-counts(dds2)
#See which count values are > 0
counts(dds2) > 0


dim(dds2)
#remove ASVs with no reads
keep2 <-rowSums(counts(dds2)) > 0
dds2 <- dds2[keep2, ]
dim(dds2)
#still the same , all ASVs had reads


gm_mean <- function(x, na.rm=TRUE){exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans2 <- apply(counts(dds2), 1, gm_mean) #counts function in deseq2 package
dds2 <- estimateSizeFactors(dds2, geoMeans = geoMeans2)


#`DESeq()` does the rest of the testing for differential expression (abundance), in this case with default testing framework, but you can actually use alternatives. 
#*Note:* The default multiple-inference correction is Benjamini-Hochberg, and occurs within the `DESeq` function.

dds_l2 <- DESeq(dds2, fitType = "local")

plotDispEsts(dds_l2)
title("Local Fit")


#perform the variance-stabilizing transformation
vst_dds2 <- varianceStabilizingTransformation(dds_l2, blind = TRUE)

#to inspect this transformation
head(assay(vst_dds2), 3)

#Plot PCA to see what's going on before transforming negative values
#DHW > 8
plotPCA(vst_dds2, intgroup = "DHW_8") + geom_vline(xintercept = 0.0, color = "gray", size = 0.5) + geom_hline(yintercept = 0.0, color = "gray", size = 0.5) + facet_grid(~DHW_8) + labs(color = NULL) + theme(legend.position = "none") + scale_color_manual(values = dhw_colors)


#Put vst data back into phyloseq

#create copy of phyloseq object
new_phylo2 <- new_phylo


#The `assay` function is used to extract the matrix of transformed values.

vst_assay2 <- assay(vst_dds2)
#this will be variance stabilised OTU table to merge with PHYLOSEQ


#QAQC

#inspect new count abundance 
abund_sums2 <-data.frame(sum = colSums(vst_assay2, na.rm = TRUE),
                         sample = colnames(vst_assay2),
                         type = "DESeq2")

ggplot(abund_sums2) +
  geom_histogram(aes(x = sum), binwidth = 20) +
  xlab("Total abundance within sample")

#Q-Q plot
qqnorm(abund_sums2$sum, main = "Q-Q Plot of Total Abundance")
qqline(abund_sums2$sum, col = "red")


#Set to zero all values less than zero.
vst_assay2[vst_assay2 < 0.0] <- 0


#merge DESeq dataset with phyloseq object

#to stay consistent with set up of original phyloseq object we will need to rearrange vst_assay so ASVs are columns
vst_assay2 <-(t(vst_assay2))
dim(vst_assay2) #36 990 same dimensions as original new_phylo
otu_table(new_phylo2) <- otu_table(vst_assay2, taxa_are_rows = FALSE)

sum(taxa_sums(new_phylo2) == 0) #0
sum(taxa_sums(new_phylo2) == 1) #0


#Now to see which symbiont taxa are driving these statistically significant differences. 
#GLM in DESeq2 (differential abundance)

#Want to know: 
  
#* Significant Taxa: Identify which taxa have significant differences in abundance between heat stress categories.
#* Log Fold Change: Understand the direction and magnitude of change for each significant taxon.
#* help understand which taxa are driving the observed differences in community composition

#*Steps to differential abundance analysis:*
  
#1. normalization
#2. model fitting
#3. hypothesis testing

#DESeq automatically creates names for coefficients of model and the can be extracted using resultsNames function

resultsNames(dds_l2)
#only shows comparisons between 1-3 and others, missing some info


## Pairwise Comparison

#Step 1: Identify Site levels

levels(sample_data(new_phylo2)$DHW_8) #NULL


#DHW_8 was not treated as a factor, making it a factor now.

sample_data(new_phylo2)$DHW_8 <- as.factor(sample_data(new_phylo2)$DHW_8)
#confirm the levels of the Site variable
dhw_levels <- levels(sample_data(new_phylo2)$DHW_8)

#Step 2: Perform Pairwise Comparisons and Add Taxa Names:
  #Use the `contrast` argument in the results function to perform pairwise comparisons between all site combinations. Gather the results including taxa names from the ITS2_type column.

#perform pairwise comparisons and add taxa names
comparisons <- combn(dhw_levels, 2, simplify = FALSE)

result_list <- lapply(comparisons, function(x) {
  contrast <- c("DHW_8", x[1], x[2])
  res <- results(dds_l2, contrast = contrast)
  res <- res[order(res$padj, na.last = NA), ]  #order by adjusted p-value
  alpha = 0.05
  sig_res <- res[which(res$padj < alpha), ]  #filter significant results
  sig_res <- cbind(as(sig_res, "data.frame"), as(tax_table(new_phylo2)[rownames(sig_res), ], "matrix"))
  sig_res$Comparison <- paste(x[1], "vs", x[2])
  sig_res$Taxa <- sig_res$ITS2_type
  return(sig_res)
})

names(result_list) <- sapply(comparisons, function(x) paste(x, collapse = "_vs_"))


#example to visualize results for a specific comparison
res_l_vs_m <- result_list[["1-3_vs_3-6"]]
head(res_l_vs_m)


#Rather than running and visualising the above example (res_Bundegi_vs_CoralBay) for every pairwise comparison individually, I want to make this efficient and easy to see in one plot.

#Step 3: Combine Results
#combine all the significant results into one data frame for visualization.

combined_results <- do.call(rbind, result_list)

write.csv(combined_results, file = "sig_diff_abund_dhw.csv")

#Then I want to arrange the pairwise plots by descending order in no. of significantly different taxa (ITS2_type). 

#Step 4: Count Significant Taxa
#count the number of significant taxa for each comparison and sort the comparisons accordingly.

comparison_counts <- combined_results %>%
  dplyr::group_by(Comparison) %>%
  dplyr::summarize(Count = n_distinct(Taxa)) %>%
  dplyr::arrange(desc(Count))

write.csv(comparison_counts, file = "comparison_counts_dhw.csv")

sorted_comparisons <- comparison_counts$Comparison
combined_results$Comparison <- factor(combined_results$Comparison, levels = sorted_comparisons)

#Step 5: Visualize Using ggplot2 with Facets
#facet by the pairwise comparisons in the sorted order.

theme_set(theme_bw())

#order Taxa (ITS2_type) by the maximum log2FoldChange for better visualization
combined_results <- combined_results %>%
  group_by(Taxa) %>%
  mutate(MaxLog2FC = max(log2FoldChange, na.rm = TRUE)) %>%
  ungroup() %>%
  arrange(desc(MaxLog2FC)) %>%
  mutate(Taxa = factor(Taxa, levels = unique(Taxa)))

##### TRY 1 ######
#title = "Differential Abundance of Taxa Across Pairwise Comparisons",
p1 <- ggplot(combined_results, aes(x = log2FoldChange, y = Taxa, color = Taxa)) +
  geom_vline(xintercept = 0.0, color = "gray", size = 0.5) +
  geom_point(size = 4) +
  facet_wrap(~ Comparison, scales = "free_y", ncol = 3) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0, vjust = 0.5),
    strip.text = element_text(face = "bold", size = 10),
    legend.position = "none",
    strip.background = element_rect(fill = "white"),
    plot.title = element_text(face = "bold", size = 20, hjust = -0.1) 
  ) +
  labs(
    x = "Log2 Fold Change",
    y = expression(italic("Cladocopium")~" taxa")) + 
  #color = expression(italic("Cladocopium") ~ "ASVs")) +
  scale_color_manual(values = ASV_colors)

#save the plot as a PNG file
ggsave("differential_abundance_plot_dhw_nocolor.png", plot = p1, width = 6, height = 6, units = "in")




######################### TRY 2 ##########
#I want positive side of the x-axis to correspond to the first DHW level in the comparison, and the negative side to correspond to the second DHW level
#Define color mapping for each comparison
dhw_comp_colors <- list(
  "3-6 vs 6-9" = c(neg = "#404040", pos = "#707070"),
  "1-3 vs 3-6" = c(neg = "#707070", pos = "#E0E0E0"),
  "1-3 vs 6-9" = c(neg = "#404040", pos = "#E0E0E0")
)

#Add a new column for the colors
combined_results <- combined_results %>%
  mutate(
    neg_color = sapply(Comparison, function(x) dhw_comp_colors[[x]]["neg"]),
    pos_color = sapply(Comparison, function(x) dhw_comp_colors[[x]]["pos"])
  )

#Generate the plot
p <- ggplot(combined_results, aes(x = log2FoldChange, y = Taxa, color = Taxa)) +
  geom_rect(
    aes(
      xmin = -Inf, xmax = 0, ymin = -Inf, ymax = Inf, fill = neg_color
    ),
    alpha = 0.8, inherit.aes = FALSE
  ) +
  geom_rect(
    aes(
      xmin = 0, xmax = Inf, ymin = -Inf, ymax = Inf, fill = pos_color
    ),
    alpha = 0.8, inherit.aes = FALSE
  ) +
  geom_vline(xintercept = 0.0, color = "gray", size = 0.5) +
  geom_point(size = 4) +
  facet_wrap(~ Comparison, scales = "free_y", ncol = 3) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5),
    strip.text = element_text(face = "bold", size = 10),
    legend.position = "none",
    strip.background = element_rect(fill = "white"),
    plot.title = element_text(face = "bold", size = 20, hjust = -0.1) 
  ) +
  labs(
    #title = "(c)",
    x = "Log2 Fold Change",
    y = expression("Taxa (" * italic("Cladocopium") * " ASVs)")
  ) +
  scale_color_manual(values = ASV_colors) +
  scale_fill_identity() 

# Save the updated plot
ggsave("differential_abundance_plot_dhw.png", plot = p, width = 8, height = 6, units = "in")

####### Additional map?#######



### Heat map
#Show a matrix of significant taxa and their fold changes across comparisons.

library(pheatmap)


heatmap_data <- combined_results %>%
  group_by(Taxa, Comparison) %>%
  summarize(log2FoldChange = mean(log2FoldChange, na.rm = TRUE)) %>%  # Use mean to resolve duplicates
  pivot_wider(
    names_from = Comparison,
    values_from = log2FoldChange,
    values_fill = 0
  ) %>%
  column_to_rownames("Taxa")

pheatmap(
  as.matrix(heatmap_data),
  color = colorRampPalette(c("blue", "white", "red"))(100),
  cluster_rows = TRUE,  # Cluster rows (taxa)
  cluster_cols = TRUE,  # Cluster columns (comparisons)
  main = "Differential Abundance Heatmap",
  display_numbers = FALSE,  # Optionally display values in cells
  fontsize_row = 8,  # Adjust row font size
  fontsize_col = 10  # Adjust column font size
)

#interesting to look at when trying to interpret the above plot "p"

############################### test - showing events category without sites ###########

#reorder columns manually
comparison_order <- c("1-3 vs 3-6", "3-6 vs 6-9", "1-3 vs 6-9")

heatmap_data2 <- combined_results %>%
  group_by(Taxa, Comparison) %>%
  summarize(log2FoldChange = mean(log2FoldChange, na.rm = TRUE)) %>%
  pivot_wider(
    names_from = Comparison,
    values_from = log2FoldChange,
    values_fill = 0
  ) %>%
  column_to_rownames("Taxa")

#make sure columns are in the right order
heatmap_data2 <- heatmap_data2[, comparison_order]

#create heatmap
pheatmap(
  as.matrix(heatmap_data2),
  color = colorRampPalette(c("blue", "white", "red"))(100),
  cluster_rows = TRUE,        # Cluster ITS2 types
  cluster_cols = FALSE,       # DO NOT cluster comparisons
  gaps_col = c(1,2),          # Add visual gaps between comparisons
  main = "Differential Abundance Heatmap (Ordered)",
  display_numbers = FALSE,
  fontsize_row = 8,
  fontsize_col = 10
)
