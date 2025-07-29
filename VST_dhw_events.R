
# using variance normalised data

# load packages
library(phyloseq); packageVersion("phyloseq") #'1.46.0'
library(ggplot2); packageVersion("ggplot2") #'3.5.1'
library(tidyverse); packageVersion("tidyverse") #‘2.0.0’

# load phyloseq object to subset from
vst_phylo.percent <- readRDS("vst_phylo.rds")

############################ Subset by No. Accum. Heat Stress Events #########################

new_phylo_l <- subset_samples(vst_phylo.percent,DHW_8=="1-3")
new_phylo_l <- prune_taxa(taxa_sums(new_phylo_l) > 0, new_phylo_l)

new_phylo_m <- subset_samples(vst_phylo.percent,DHW_8=="3-6")
new_phylo_m <- prune_taxa(taxa_sums(new_phylo_m) > 0, new_phylo_m)

new_phylo_h <- subset_samples(vst_phylo.percent,DHW_8=="6-9")
new_phylo_h <- prune_taxa(taxa_sums(new_phylo_h) > 0, new_phylo_h)


## Determining Prevalence and Relative Abundance of each ITS2 Type across Samples that they are found in  

#extract ASV table
otu_table_data_l <- otu_table(new_phylo_l)
otu_table_data_l <- t(otu_table_data_l)

#convert to data frames
otu_df_l <- as.data.frame(otu_table_data_l)

#ensure the row names (OTU IDs) match between the ASV and tax tables
otu_df_l$OTU <- rownames(otu_df_l)

#compute prevalence of each feature, store as df
prevdf_l = apply(X = otu_table(new_phylo_l),
                 MARGIN = ifelse(taxa_are_rows(new_phylo_l), yes = 1, no = 2),
                 FUN = function(x){sum(x > 0)})

#add taxonomy and total read counts to this df
prevdf_l = data.frame(Prevalence = prevdf_l,
                      TotalAbundance = taxa_sums(new_phylo_l),
                      tax_table(new_phylo_l)) %>%
  rownames_to_column() %>%
  dplyr::rename(OTU='rowname')

prev_otu_l <- prevdf_l %>%
  left_join( otu_df_l, by = "OTU") %>%
  tidyr::pivot_longer(cols = c(5:20), names_to = "Sample", values_to = "ASV_count_per_sample") %>%
  filter(ASV_count_per_sample != 0) %>%
  group_by(Sample) %>%
  mutate(Total_reads_per_sample = sum(ASV_count_per_sample)) %>%
  mutate(ITS2_diversity_per_sample = n_distinct(ITS2_type)) 

#find abundance of each ITS2 type that ASVs map to within each sample
prev_otu2_l <- prev_otu_l %>%
  group_by(Sample, ITS2_type) %>%
  mutate(Total_ITS2_count_per_sample = sum(ASV_count_per_sample)) %>%
  mutate(abundance_per_sample = Total_ITS2_count_per_sample / Total_reads_per_sample) %>%
  summarize(abundance_per_sample = mean(abundance_per_sample), .groups = "keep") %>%
  ungroup()

#find relative abundance of each ITS2 type across samples that it's found in (plus SE)
prev_otu3_l <- prev_otu2_l %>%
  group_by(ITS2_type) %>%
  summarize(mean_abund = mean_se(abundance_per_sample)) %>%
  mutate(se_abund = mean_abund$y - mean_abund$ymin) %>%
  arrange(desc(mean_abund))

#extract prevalence data for each ITS2 type (not just specific ASVs) and combine into df with relative abundance
prev_otu4_l <- prev_otu_l %>%
  group_by(ITS2_type) %>%
  mutate(Prevalence = n_distinct(Sample)) %>%
  summarize(Prevalence = max(Prevalence)) %>%
  left_join(prev_otu3_l, by = "ITS2_type") %>%
  arrange(desc(Prevalence)) 

write.csv(prev_otu4_l, file = "vst_prevalence_abundance_low.csv")



# Med

new_phylo_m 

## Determining Prevalence and Relative Abundance of each ITS2 Type across Samples that they are found in  

#extract ASV table
otu_table_data_m <- otu_table(new_phylo_m)
otu_table_data_m <- t(otu_table_data_m)

#convert to data frames
otu_df_m <- as.data.frame(otu_table_data_m)

#ensure the row names (OTU IDs) match between the ASV and tax tables
otu_df_m$OTU <- rownames(otu_df_m)

#compute prevalence of each feature, store as df
prevdf_m = apply(X = otu_table(new_phylo_m),
                 MARGIN = ifelse(taxa_are_rows(new_phylo_m), yes = 1, no = 2),
                 FUN = function(x){sum(x > 0)})

#add taxonomy and total read counts to this df
prevdf_m = data.frame(Prevalence = prevdf_m,
                      TotalAbundance = taxa_sums(new_phylo_m),
                      tax_table(new_phylo_m)) %>%
  rownames_to_column() %>%
  dplyr::rename(OTU='rowname')

prev_otu_m <- prevdf_m %>%
  left_join( otu_df_m, by = "OTU") %>%
  tidyr::pivot_longer(cols = c(5:20), names_to = "Sample", values_to = "ASV_count_per_sample") %>%
  filter(ASV_count_per_sample != 0) %>%
  group_by(Sample) %>%
  mutate(Total_reads_per_sample = sum(ASV_count_per_sample)) %>%
  mutate(ITS2_diversity_per_sample = n_distinct(ITS2_type)) 

#find abundance of each ITS2 type that ASVs map to within each sample
prev_otu2_m <- prev_otu_m %>%
  group_by(Sample, ITS2_type) %>%
  mutate(Total_ITS2_count_per_sample = sum(ASV_count_per_sample)) %>%
  mutate(abundance_per_sample = Total_ITS2_count_per_sample / Total_reads_per_sample) %>%
  summarize(abundance_per_sample = mean(abundance_per_sample), .groups = "keep") %>%
  ungroup()

#find relative abundance of each ITS2 type across samples that it's found in (plus SE)
prev_otu3_m <- prev_otu2_m %>%
  group_by(ITS2_type) %>%
  summarize(mean_abund = mean_se(abundance_per_sample)) %>%
  mutate(se_abund = mean_abund$y - mean_abund$ymin) %>%
  arrange(desc(mean_abund))

#extract prevalence data for each ITS2 type (not just specific ASVs) and combine into df with relative abundance
prev_otu4_m <- prev_otu_m %>%
  group_by(ITS2_type) %>%
  mutate(Prevalence = n_distinct(Sample)) %>%
  summarize(Prevalence = max(Prevalence)) %>%
  left_join(prev_otu3_m, by = "ITS2_type") %>%
  arrange(desc(Prevalence)) 

write.csv(prev_otu4_m, file = "vst_prevalence_abundance_med.csv")

# High

new_phylo_h 

## Determining Prevalence and Relative Abundance of each ITS2 Type across Samples that they are found in  

#extract ASV table
otu_table_data_h <- otu_table(new_phylo_h)
otu_table_data_h <- t(otu_table_data_h)

#convert to data frames
otu_df_h <- as.data.frame(otu_table_data_h)

#ensure the row names (OTU IDs) match between the ASV and tax tables
otu_df_h$OTU <- rownames(otu_df_h)

#compute prevalence of each feature, store as df
prevdf_h = apply(X = otu_table(new_phylo_h),
                 MARGIN = ifelse(taxa_are_rows(new_phylo_h), yes = 1, no = 2),
                 FUN = function(x){sum(x > 0)})

#add taxonomy and total read counts to this df
prevdf_h = data.frame(Prevalence = prevdf_h,
                      TotalAbundance = taxa_sums(new_phylo_h),
                      tax_table(new_phylo_h)) %>%
  rownames_to_column() %>%
  dplyr::rename(OTU='rowname')

prev_otu_h <- prevdf_h %>%
  left_join( otu_df_h, by = "OTU") %>%
  tidyr::pivot_longer(cols = c(5:9), names_to = "Sample", values_to = "ASV_count_per_sample") %>%
  filter(ASV_count_per_sample != 0) %>%
  group_by(Sample) %>%
  mutate(Total_reads_per_sample = sum(ASV_count_per_sample)) %>%
  mutate(ITS2_diversity_per_sample = n_distinct(ITS2_type)) 

#find abundance of each ITS2 type that ASVs map to within each sample
prev_otu2_h <- prev_otu_h %>%
  group_by(Sample, ITS2_type) %>%
  mutate(Total_ITS2_count_per_sample = sum(ASV_count_per_sample)) %>%
  mutate(abundance_per_sample = Total_ITS2_count_per_sample / Total_reads_per_sample) %>%
  summarize(abundance_per_sample = mean(abundance_per_sample), .groups = "keep") %>%
  ungroup()

#find relative abundance of each ITS2 type across samples that it's found in (plus SE)
prev_otu3_h <- prev_otu2_h %>%
  group_by(ITS2_type) %>%
  summarize(mean_abund = mean_se(abundance_per_sample)) %>%
  mutate(se_abund = mean_abund$y - mean_abund$ymin) %>%
  arrange(desc(mean_abund))

#extract prevalence data for each ITS2 type (not just specific ASVs) and combine into df with relative abundance
prev_otu4_h <- prev_otu_h %>%
  group_by(ITS2_type) %>%
  mutate(Prevalence = n_distinct(Sample)) %>%
  summarize(Prevalence = max(Prevalence)) %>%
  left_join(prev_otu3_h, by = "ITS2_type") %>%
  arrange(desc(Prevalence)) 

write.csv(prev_otu4_h, file = "vst_prevalence_abundance_high.csv")

############################ ASVs by group #########################

# Which ASVs per group (DHW for now)

l.ASVs <- rownames(data.frame(tax_table(new_phylo_l)))
m.ASVs <- rownames(data.frame(tax_table(new_phylo_m)))
h.ASVs <- rownames(data.frame(tax_table(new_phylo_h)))


##################### Venn Diagram of ASVs by DHW #################

## caption: Venn diagram shows the number of Symbiodiniaceae ASVs present in each sample type, as well as the amount of overlap between and among sample types.


# Load necessary packages
library(VennDiagram)
library(eulerr)
library(gridExtra)

dhw_colors <- c("1-3" = "#E0E0E0", 
                "3-6" = "#707070", 
                "6-9" = "#404040")


ASV_colors <- c(
  "C15" = "#E69F00",         
  "C116" = "#D73027",        
  "C15.6" = "#F28E2B",       
  "C15h" = "#FFBF00",        
  "C15.9" = "#E68310",       
  "C56" = "#66A61E",         
  "C55" = "#1B9E77",         
  "C15.8" = "#FFA500",       
  "C15.2.2" = "#E49B0F",     
  "C15L" = "#F5B041",        
  "C3" = "#0569AD",          
  "C3w" = "#6A9EC2",         
  "C15.2.1" = "#F29E38",     
  "C15a" = "#FFB732",        
  "C15.7" = "#FF8C00",       
  "C15d" = "#FFB84D",        
  "C60" = "#228B22",         
  "C15i" = "#FFC857",        
  "C15f" = "#FFAF3F",        
  "C15m" = "#E89928",        
  "C91b" = "#78C679",        
  "C1" = "#377EB8",          
  "C15.1" = "#EAA800",       
  "C15c" = "#FFB74D",        
  "C15g" = "#FDB462",        
  "C3z" = "#5086B4",         
  "C54a" = "#006D2C",        
  "C91_C15" = "#FFDD71",     
  "C91a" = "#31A354",        
  "Cspe" = "#41AB5D",        
  "Cspf" = "#2E8B57",        
  "C15j" = "#FFAC3E",        
  "C48" = "#5A9E6F",         
  "C50"= "#4DBD33",          
  "Cspc" = "#3E9141",        
  "D5" = "#006400"           
)





# Determine the intersects (i.e. which ASVs are shared) between different heat stress for Venn plotting
lm <- intersect(l.ASVs,m.ASVs)
lmh <- intersect(lm,h.ASVs)
lm.only <- length(lm)-length(lmh)
lmh.length <- length(lmh)
lh <- intersect(l.ASVs,h.ASVs)
lh.only <- length(lh)-length(lmh)
mh <- intersect(m.ASVs,h.ASVs)
mh.only <- length(mh)-length(lmh)
low <- length(l.ASVs) - lm.only - lh.only - lmh.length
med <- length(m.ASVs) - mh.only - lm.only - lmh.length
high <- length(h.ASVs) - mh.only - lh.only - lmh.length

#what ITS2_types do these shared ASVs map to?

lmh_ITS_types <- unique(data.frame(tax_table(vst_phylo.percent))[rownames(data.frame(tax_table(vst_phylo.percent))) %in% lmh, "ITS2_type"])
lm_only_ITS_types <- unique(data.frame(tax_table(vst_phylo.percent))[rownames(data.frame(tax_table(vst_phylo.percent))) %in% lm, "ITS2_type"])
lh_only_ITS_types <- unique(data.frame(tax_table(vst_phylo.percent))[rownames(data.frame(tax_table(vst_phylo.percent))) %in% lh, "ITS2_type"])
mh_only_ITS_types <- unique(data.frame(tax_table(vst_phylo.percent))[rownames(data.frame(tax_table(vst_phylo.percent))) %in% mh, "ITS2_type"])

#what ITS2_types do unique ASVs map to?
# Unique ASVs for each group
l_u_ASVs <- setdiff(l.ASVs, union(m.ASVs, h.ASVs)) 
m_u_ASVs <- setdiff(m.ASVs, union(l.ASVs, h.ASVs)) 
h_u_ASVs <- setdiff(h.ASVs, union(l.ASVs, m.ASVs)) 

# Extract ITS2 types for unique ASVs
l_only_ITS <- unique(data.frame(tax_table(vst_phylo.percent))[rownames(data.frame(tax_table(vst_phylo.percent))) %in% l_u_ASVs, "ITS2_type"])
m_only_ITS <- unique(data.frame(tax_table(vst_phylo.percent))[rownames(data.frame(tax_table(vst_phylo.percent))) %in% m_u_ASVs, "ITS2_type"])
h_only_ITS <- unique(data.frame(tax_table(vst_phylo.percent))[rownames(data.frame(tax_table(vst_phylo.percent))) %in% h_u_ASVs, "ITS2_type"])



#Summarise
vst_summary_ITS_DHW <- list(
  "Shared Among All (LMH)" = lmh_ITS_types,
  "Shared Low & Medium Only" = lm_only_ITS_types,
  "Shared Low & High Only" = lh_only_ITS_types,
  "Shared Medium & High Only" = mh_only_ITS_types,
  "Low Only" = l_only_ITS,
  "Medium Only" = m_only_ITS,
  "High Only" = h_only_ITS
)



## Make figure - jpg
#jpeg(filename="DHW_VD_ASV.jpg", 
#     width = 4, height = 4, units="in",res = 300) # Open jpg
# Create Venn object
VD_DHW <- euler(c("1-3 Events" = (low), "3-6 Events" = (med), "6-9 Events" = (high), "1-3 Events&3-6 Events" = (lm.only), "1-3 Events&6-9 Events" = (lh.only), "3-6 Events&6-9 Events" = (mh.only), "1-3 Events&3-6 Events&6-9 Events" = (lmh.length)))
# Plot Venn object
VD_DHW_plot <- plot(VD_DHW, quantities = TRUE, font=1, cex=1, alpha=0.8,
                    fill=dhw_colors,col=dhw_colors,border=dhw_colors,lwd=c(2,2,2)) 
#dev.off() # Close jpg

# Save the updated plot
ggsave("VD_DHW.png", plot = VD_DHW_plot, width = 6, height = 6, units = "in")


##################################### Pie charts for frequency of ITS2_types that are unique or shared ##############################

#1. Prepare data

#lm
lm_ITS <- data.frame(
  ITS2_type = data.frame(tax_table(vst_phylo.percent))[rownames(data.frame(tax_table(vst_phylo.percent))) %in% lm, "ITS2_type"]
)

#lh
lh_ITS <- data.frame(
  ITS2_type = data.frame(tax_table(vst_phylo.percent))[rownames(data.frame(tax_table(vst_phylo.percent))) %in% lh, "ITS2_type"]
)

#mh
mh_ITS <- data.frame(
  ITS2_type = data.frame(tax_table(vst_phylo.percent))[rownames(data.frame(tax_table(vst_phylo.percent))) %in% mh, "ITS2_type"]
)

#lmh
lmh_ITS <- data.frame(
  ITS2_type = data.frame(tax_table(vst_phylo.percent))[rownames(data.frame(tax_table(vst_phylo.percent))) %in% lmh, "ITS2_type"]
)

#l only
l_ITS <- data.frame(
  ITS2_type = data.frame(tax_table(vst_phylo.percent))[rownames(data.frame(tax_table(vst_phylo.percent))) %in% l_u_ASVs, "ITS2_type"]
)

#m only
m_ITS <- data.frame(
  ITS2_type = data.frame(tax_table(vst_phylo.percent))[rownames(data.frame(tax_table(vst_phylo.percent))) %in% m_u_ASVs, "ITS2_type"]
)

#h only 
h_ITS <- data.frame(
  ITS2_type = data.frame(tax_table(vst_phylo.percent))[rownames(data.frame(tax_table(vst_phylo.percent))) %in% h_u_ASVs, "ITS2_type"]
)

#l all
l_ITS_all <- data.frame(ITS2_type = data.frame(tax_table(new_phylo_l))$ITS2_type)

#m all 
m_ITS_all <- data.frame(ITS2_type = data.frame(tax_table(new_phylo_m))$ITS2_type)

#h all
h_ITS_all <- data.frame(ITS2_type = data.frame(tax_table(new_phylo_h))$ITS2_type)



# Count occurrences of each ITS2 type
#lm
lm_ITS_summary <- as.data.frame(table(lm_ITS$ITS2_type))
colnames(lm_ITS_summary) <- c("ITS2_type", "Frequency")

lm_ITS_summary <- lm_ITS_summary %>%
  mutate(
    Fraction = Frequency / sum(Frequency))

#lh
lh_ITS_summary <- as.data.frame(table(lh_ITS$ITS2_type))
colnames(lh_ITS_summary) <- c("ITS2_type", "Frequency")

lh_ITS_summary <- lh_ITS_summary %>%
  mutate(
    Fraction = Frequency / sum(Frequency))

#mh
mh_ITS_summary <- as.data.frame(table(mh_ITS$ITS2_type))
colnames(mh_ITS_summary) <- c("ITS2_type", "Frequency")

mh_ITS_summary <- mh_ITS_summary %>%
  mutate(
    Fraction = Frequency / sum(Frequency))

#lmh
lmh_ITS_summary <- as.data.frame(table(lmh_ITS$ITS2_type))
colnames(lmh_ITS_summary) <- c("ITS2_type", "Frequency")

lmh_ITS_summary <- lmh_ITS_summary %>%
  mutate(
    Fraction = Frequency / sum(Frequency))

#l only
l_ITS_summary <- as.data.frame(table(l_ITS$ITS2_type))
colnames(l_ITS_summary) <- c("ITS2_type", "Frequency")

l_ITS_summary <- l_ITS_summary %>%
  mutate(
    Fraction = Frequency / sum(Frequency))

#m only
m_ITS_summary <- as.data.frame(table(m_ITS$ITS2_type))
colnames(m_ITS_summary) <- c("ITS2_type", "Frequency")

m_ITS_summary <- m_ITS_summary %>%
  mutate(
    Fraction = Frequency / sum(Frequency))

#h only
h_ITS_summary <- as.data.frame(table(h_ITS$ITS2_type))
colnames(h_ITS_summary) <- c("ITS2_type", "Frequency")

h_ITS_summary <- h_ITS_summary %>%
  mutate(
    Fraction = Frequency / sum(Frequency))


#l all
l_ITS_summary_all <- as.data.frame(table(l_ITS_all$ITS2_type))
colnames(l_ITS_summary_all) <- c("ITS2_type", "Frequency")

l_ITS_summary_all <- l_ITS_summary_all %>% 
  mutate(
    Fraction = Frequency / sum(Frequency))

#m all 
m_ITS_summary_all <- as.data.frame(table(m_ITS_all$ITS2_type))
colnames(m_ITS_summary_all) <- c("ITS2_type", "Frequency")

m_ITS_summary_all <- m_ITS_summary_all %>% 
  mutate(
    Fraction = Frequency / sum(Frequency))


#h all
h_ITS_summary_all <- as.data.frame(table(h_ITS_all$ITS2_type))
colnames(h_ITS_summary_all) <- c("ITS2_type", "Frequency")

h_ITS_summary_all <- h_ITS_summary_all %>% 
  mutate(
    Fraction = Frequency / sum(Frequency))



#2. Create Pie chart

pie_LM <- ggplot(lm_ITS_summary, aes(x = "", y = Frequency, fill = ITS2_type)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar(theta = "y") +
  labs(title = "LM Segment",
       fill = expression("Taxa (" * italic("Cladocopium") * " ASVs)")) +
  theme_void() +
  scale_fill_manual(values = ASV_colors) +
  theme(legend.position = "right")


pie_LH <- ggplot(lh_ITS_summary, aes(x = "", y = Frequency, fill = ITS2_type)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar(theta = "y") +
  labs(title = "LH Segment",
       fill = expression("Taxa (" * italic("Cladocopium") * " ASVs)")) +
  theme_void() +
  scale_fill_manual(values = ASV_colors) +
  theme(legend.position = "right")


pie_MH <- ggplot(mh_ITS_summary, aes(x = "", y = Frequency, fill = ITS2_type)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar(theta = "y") +
  labs(title = "MH Segment",
       fill = expression("Taxa (" * italic("Cladocopium") * " ASVs)")) +
  theme_void() +
  scale_fill_manual(values = ASV_colors) +
  theme(legend.position = "right")

pie_LMH <- ggplot(lmh_ITS_summary, aes(x = "", y = Frequency, fill = ITS2_type)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar(theta = "y") +
  labs(title = "(b)",
       fill = expression("Shared " * italic("Cladocopium") * " ASVs")) +
  theme_void() +
  scale_fill_manual(values = ASV_colors) +
  theme(legend.position = "right")

#without labels
ggplot(lmh_ITS_summary, aes(x = "", y = Frequency, fill = ITS2_type)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar(theta = "y") +
  theme_void() +
  scale_fill_manual(values = ASV_colors) +
  theme(
    legend.position = "none", # Remove legend
    plot.title = element_blank() # Remove title
  )


pie_L <- ggplot(l_ITS_summary, aes(x = "", y = Frequency, fill = ITS2_type)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar(theta = "y") +
  labs(title = "Low (1-3 Events)",
       fill = expression("Taxa (" * italic("Cladocopium") * " ASVs)")) +
  theme_void() +
  scale_fill_manual(values = ASV_colors) +
  theme(legend.position = "right")

#without labels
pie_L <- ggplot(l_ITS_summary, aes(x = "", y = Frequency, fill = ITS2_type)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar(theta = "y") +
  theme_void() +
  scale_fill_manual(values = ASV_colors) +
  theme(
    legend.position = "none", # Remove legend
    plot.title = element_blank() # Remove title
  )


pie_L_all <- ggplot(l_ITS_summary_all, aes(x = "", y = Frequency, fill = ITS2_type)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar(theta = "y") +
  labs(title = "Low (1-3 Events)",
       fill = expression("Taxa (" * italic("Cladocopium") * " ASVs)")) +
  theme_void() +
  scale_fill_manual(values = ASV_colors) +
  theme(legend.position = "right")


pie_M <- ggplot(m_ITS_summary, aes(x = "", y = Frequency, fill = ITS2_type)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar(theta = "y") +
  labs(title = "Medium (3-6 Events)",
       fill = expression("Taxa (" * italic("Cladocopium") * " ASVs)")) +
  theme_void() +
  scale_fill_manual(values = ASV_colors) +
  theme(legend.position = "right")

#without labels
pie_M <- ggplot(m_ITS_summary, aes(x = "", y = Frequency, fill = ITS2_type)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar(theta = "y") +
  theme_void() +
  scale_fill_manual(values = ASV_colors) +
  theme(
    legend.position = "none", # Remove legend
    plot.title = element_blank() # Remove title
  )


pie_M_all <- ggplot(m_ITS_summary_all, aes(x = "", y = Frequency, fill = ITS2_type)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar(theta = "y") +
  labs(title = "Medium (3-6 Events)",
       fill = expression("Taxa (" * italic("Cladocopium") * " ASVs)")) +
  theme_void() +
  scale_fill_manual(values = ASV_colors) +
  theme(legend.position = "right")


pie_H <- ggplot(h_ITS_summary, aes(x = "", y = Frequency, fill = ITS2_type)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar(theta = "y") +
  labs(title = "High (6-9 Events)",
       fill = expression("Taxa (" * italic("Cladocopium") * " ASVs)")) +
  theme_void() +
  scale_fill_manual(values = ASV_colors) +
  theme(legend.position = "right")

#without labels
pie_H <- ggplot(h_ITS_summary, aes(x = "", y = Frequency, fill = ITS2_type)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar(theta = "y") +
  theme_void() +
  scale_fill_manual(values = ASV_colors) +
  theme(
    legend.position = "none", # Remove legend
    plot.title = element_blank() # Remove title
  )


pie_H_all <- ggplot(h_ITS_summary_all, aes(x = "", y = Frequency, fill = ITS2_type)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar(theta = "y") +
  labs(title = "High (6-9 Events)",
       fill = expression("Taxa (" * italic("Cladocopium") * " ASVs)")) +
  theme_void() +
  scale_fill_manual(values = ASV_colors) +
  theme(legend.position = "right")


################ Combine Venn Diagram and Pie Charts ###############


# Combine Venn diagram with pie charts
final_VD_plot <- ggdraw() +
  draw_plot(VD_DHW_plot, x = 0, y = 0, width = 1, height = 1) + # Venn diagram
  draw_plot(pie_L, x = 0.08, y = 0.08, width = 0.25, height = 0.25) + # Bottom-left for Low
  draw_plot(pie_M, x = 0.75, y = 0.08, width = 0.25, height = 0.25) + # Bottom-right for Medium
  draw_plot(pie_H, x = 0.55, y = 0.76, width = 0.25, height = 0.25) + # Top-right for High
  #draw_plot(pie_LMH, x = 0.75, y = 0.75, width = 0.25, height = 0.25) + # Far upper-right for Shared ASVs
  draw_label(
    "(a)",                     # Title text
    x = 0.05,                  # Horizontal position
    y = 1,                     # Vertical position (top of the plot)
    hjust = 0,                 # Left-align the label
    vjust = 1,                 # Align to top
    fontface = "bold",         # Bold text
    size = 18                  # Font size
  )

# Save the final plot
ggsave("vst_VD_DHW_with_pies.png", plot = final_VD_plot, width = 8, height = 8, units = "in")



## Subset Plot of Shared ASVs

pie_LMH <- ggplot(lmh_ITS_summary, aes(x = "", y = Frequency, fill = ITS2_type)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar(theta = "y") +
  theme_void() +
  scale_fill_manual(
    values = ASV_colors, 
    #name = expression(italic("Cladocopium") ~ "ASVs") # Match legend title
  ) +
  theme(
    legend.position = "right",                # Match legend position
    legend.key.size = unit(0.5, "cm"),        # Match legend key size
    legend.text = element_text(size = 8),    # Match legend text size
    legend.title = element_text(size = 8),   # Match legend title size
    plot.title = element_text(face = "bold") # Make title bold
  ) +
  labs(
    title = "(b)",        # Add title
    fill = expression("Shared " * italic("Cladocopium") * " ASVs") 
  )


# Add a border around the entire plot using ggdraw
pie_LMH_with_border <- ggdraw() +
  draw_plot(pie_LMH, x = 0.02, y = 0, width = 0.95, height = 1) + # Add the pie chart
  theme(
    plot.background = element_rect(color = "black", size = 1) # Add black border
  )


# Save the final plot
ggsave("vst_VD_DHW_subset.png", plot = pie_LMH_with_border , width = 5, height = 5, units = "in")
