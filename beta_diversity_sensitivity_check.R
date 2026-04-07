## Sensitivity check: weighted UniFrac from raw counts vs relative abundances

library(phyloseq)
library(vegan)
library(dplyr)
library(ggplot2)
library(cowplot)


#raw-count workflow
set.seed(711)
vPorWUF_raw <- UniFrac(new_phylo, weighted = TRUE, normalized = TRUE, parallel = FALSE)

set.seed(711)
nmds_raw <- ordinate(new_phylo, method = "NMDS", distance = vPorWUF_raw)

ad_main_raw <- adonis2(
  vPorWUF_raw ~ Possible_Species + DHW_8 + clim_temp_range_c,
  data = meta1,
  permutations = 999,
  by = "margin"
)


#relative-abundance workflow
phy_rel <- transform_sample_counts(new_phylo, function(x) x / sum(x))

set.seed(711)
vPorWUF_rel <- UniFrac(phy_rel, weighted = TRUE, normalized = TRUE, parallel = FALSE)

set.seed(711)
nmds_rel <- ordinate(phy_rel, method = "NMDS", distance = vPorWUF_rel)

#rebuild metadata to match the same samples/order
meta_rel <- data.frame(sample_data(phy_rel))
rownames(meta_rel) <- sample_names(phy_rel)
meta_rel$Possible_Species <- factor(meta_rel$Possible_Species)
meta_rel$Site <- factor(meta_rel$Site)
meta_rel$DHW_8 <- factor(meta_rel$DHW_8, levels = c("1-3", "3-6", "6-9"))
meta_rel <- meta_rel[labels(vPorWUF_rel), , drop = FALSE]

ad_main_rel <- adonis2(
  vPorWUF_rel ~ Possible_Species + DHW_8 + clim_temp_range_c,
  data = meta_rel,
  permutations = 999,
  by = "margin"
)


#PERMDISP comparison
beta_host_raw <- betadisper(vPorWUF_raw, meta1$Possible_Species, bias.adjust = TRUE)
beta_dhw_raw  <- betadisper(vPorWUF_raw,  meta1$DHW_8,            bias.adjust = TRUE)

beta_host_rel <- betadisper(vPorWUF_rel, meta_rel$Possible_Species, bias.adjust = TRUE)
beta_dhw_rel  <- betadisper(vPorWUF_rel,  meta_rel$DHW_8,            bias.adjust = TRUE)

permdisp_compare <- bind_rows(
  data.frame(
    Workflow = "Raw counts",
    Factor = "Host species",
    F = permutest(beta_host_raw, permutations = 999)$tab[1, "F"],
    p_value = permutest(beta_host_raw, permutations = 999)$tab[1, "Pr(>F)"]
  ),
  data.frame(
    Workflow = "Raw counts",
    Factor = "DHW",
    F = permutest(beta_dhw_raw, permutations = 999)$tab[1, "F"],
    p_value = permutest(beta_dhw_raw, permutations = 999)$tab[1, "Pr(>F)"]
  ),
  data.frame(
    Workflow = "Relative abundance",
    Factor = "Host species",
    F = permutest(beta_host_rel, permutations = 999)$tab[1, "F"],
    p_value = permutest(beta_host_rel, permutations = 999)$tab[1, "Pr(>F)"]
  ),
  data.frame(
    Workflow = "Relative abundance",
    Factor = "DHW",
    F = permutest(beta_dhw_rel, permutations = 999)$tab[1, "F"],
    p_value = permutest(beta_dhw_rel, permutations = 999)$tab[1, "Pr(>F)"]
  )
)

print(permdisp_compare)


#clean PERMANOVA extraction
extract_main_terms <- function(adonis_obj, workflow_name) {
  out <- as.data.frame(adonis_obj)
  out$Term <- rownames(out)
  rownames(out) <- NULL
  out <- out[, c("Term", "Df", "SumOfSqs", "R2", "F", "Pr(>F)")]
  names(out) <- c("Term", "Df", "SumOfSqs", "R_squared", "pseudo_F", "p_value")
  out$Workflow <- workflow_name
  out
}

perm_compare <- bind_rows(
  extract_main_terms(ad_main_raw, "Raw counts"),
  extract_main_terms(ad_main_rel, "Relative abundance")
)

print(perm_compare)
write.csv(perm_compare, "weighted_unifrac_raw_vs_relative_PERMANOVA.csv", row.names = FALSE)


#compare NMDS stress values
stress_compare <- data.frame(
  Workflow = c("Raw counts", "Relative abundance"),
  Stress = c(nmds_raw$stress, nmds_rel$stress)
)

print(stress_compare)
write.csv(stress_compare, "weighted_unifrac_raw_vs_relative_stress.csv", row.names = FALSE)


#plot side-by-side NMDS
scores_raw <- as.data.frame(scores(nmds_raw, display = "sites"))
scores_raw$SampleID <- rownames(scores_raw)
scores_raw <- left_join(scores_raw, meta1 %>% tibble::rownames_to_column("SampleID"), by = "SampleID")
scores_raw$Workflow <- "Raw counts"

scores_rel <- as.data.frame(scores(nmds_rel, display = "sites"))
scores_rel$SampleID <- rownames(scores_rel)
scores_rel <- left_join(scores_rel, meta_rel %>% tibble::rownames_to_column("SampleID"), by = "SampleID")
scores_rel$Workflow <- "Relative abundance"

scores_both <- bind_rows(scores_raw, scores_rel)

#shared axis limits
xlims <- range(scores_both$NMDS1, na.rm = TRUE)
ylims <- range(scores_both$NMDS2, na.rm = TRUE)

p_compare <- ggplot(scores_both, aes(NMDS1, NMDS2)) +
  geom_point(aes(shape = Possible_Species, fill = clim_temp_range_c),
             size = 3.5, color = "black", stroke = 0.8) +
  scale_shape_manual(values = c("P. lobata" = 21, "P. lutea" = 24)) +
  scale_fill_gradient(low = "#dbe9f6", high = "#08519c",
                      name = "Temperature range (°C)") +
  facet_wrap(~ Workflow, nrow = 1) +
  coord_equal(xlim = xlims, ylim = ylims) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(face = "bold"),
    legend.position = "bottom"
  ) +
  labs(x = "NMDS1", y = "NMDS2")

print(p_compare)
ggsave("weighted_unifrac_raw_vs_relative_NMDS.png",
       p_compare, width = 10, height = 5.5, units = "in", dpi = 300)