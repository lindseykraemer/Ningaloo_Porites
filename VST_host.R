#load packages
library(phyloseq); packageVersion("phyloseq") #'1.46.0'
library(ggplot2); packageVersion("ggplot2") #'3.5.1'
library(tidyverse); packageVersion("tidyverse") #‘2.0.0’

#load phyloseq object to subset from vst normalised count data
vst_phylo.percent <- readRDS("vst_phylo.rds")

#load host colours
species_colors <- c(
  "P. lobata" = "#009E73",           
  "P. lutea" =  "#D55E00"           
)



#load colour scheme for ASVs
#(Cladocopium taxa: C15 variants (Orange / Yellow Tones), Thermotolerant Types (Red Tones), Thermosensitive Types (Blue Tones), Unknown Thermal Tolerance (Green Tones))
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


########################## Subset by species ################################


#subset vst_phylo.percent to each species
new_phylo_lob <- subset_samples(vst_phylo.percent,Possible_Species=="P. lobata")
new_phylo_lob <- prune_taxa(taxa_sums(new_phylo_lob) > 0, new_phylo_lob)


new_phylo_lut <- subset_samples(vst_phylo.percent,Possible_Species=="P. lutea")
new_phylo_lut <- prune_taxa(taxa_sums(new_phylo_lut) > 0, new_phylo_lut)


####### Finding out proportion and frequency of taxa per host spp ######

#calculate total number of sequences
total_seqs <- sum(taxa_sums(vst_phylo.percent))
lob_seqs <- sum(taxa_sums(new_phylo_lob))
lut_seqs <- sum(taxa_sums(new_phylo_lut))

############################ ASVs by group #########################

#number of ASVs per group

n.lob.ASVs <- ntaxa(new_phylo_lob)
n.lut.ASVs <- ntaxa(new_phylo_lut)

#which ASVs per group 

lob.ASVs <- rownames(data.frame(tax_table(new_phylo_lob)))
lut.ASVs <- rownames(data.frame(tax_table(new_phylo_lut)))


##################### ITS2_type by group #######################


unique(data.frame(tax_table(new_phylo_lob))$ITS2_type)
unique(data.frame(tax_table(new_phylo_lut))$ITS2_type)


##################### Venn Diagram of ASVs by Species #################

#load packages
library(VennDiagram)
library(eulerr)
library(gridExtra)

#extract ASVs belonging to each species (below vectors from above^^)
lut.ASVs <- rownames(data.frame(tax_table(new_phylo_lut)))
lob.ASVs <- rownames(data.frame(tax_table(new_phylo_lob)))

#find shared and unique ASVs
lutlob <- intersect(lut.ASVs,lob.ASVs)
lutlob_length <- length(lutlob)
lut_only <- length(lut.ASVs) - lutlob_length
lob_only <- length(lob.ASVs) - lutlob_length
lut_u_ASVs <- setdiff(lut.ASVs, lob.ASVs)
lob_u_ASVs <- setdiff(lob.ASVs, lut.ASVs)

#what ITS2_types do these shared and unique ASVs map to?
lutlob_ITS_types <- unique(data.frame(tax_table(vst_phylo.percent))[rownames(data.frame(tax_table(vst_phylo.percent))) %in% lutlob, "ITS2_type"])
lut_only_ITS_types <- unique(data.frame(tax_table(vst_phylo.percent))[rownames(data.frame(tax_table(vst_phylo.percent))) %in% lut_u_ASVs, "ITS2_type"])
lob_only_ITS_types <- unique(data.frame(tax_table(vst_phylo.percent))[rownames(data.frame(tax_table(vst_phylo.percent))) %in% lob_u_ASVs, "ITS2_type"])

#organise for pie charts and VD
lutlob_ITS <- data.frame(
  ITS2_type = data.frame(tax_table(vst_phylo.percent))[rownames(data.frame(tax_table(vst_phylo.percent))) %in% lutlob, "ITS2_type"]
)
lut_only_ITS <- data.frame(
  ITS2_type = data.frame(tax_table(vst_phylo.percent))[rownames(data.frame(tax_table(vst_phylo.percent))) %in% lut_u_ASVs, "ITS2_type"]
)
lob_only_ITS <- data.frame(
  ITS2_type = data.frame(tax_table(vst_phylo.percent))[rownames(data.frame(tax_table(vst_phylo.percent))) %in% lob_u_ASVs, "ITS2_type"]
)

#lut all
lut_ITS_all <- data.frame(ITS2_type = data.frame(tax_table(new_phylo_lut))$ITS2_type)

#lob all 
lob_ITS_all <- data.frame(ITS2_type = data.frame(tax_table(new_phylo_lob))$ITS2_type)


## Summarise

#shared ITS2 types
lutlob_ITS_summary <- as.data.frame(table(lutlob_ITS$ITS2_type))
colnames(lutlob_ITS_summary) <- c("ITS2_type", "Frequency")
lutlob_ITS_summary <- lutlob_ITS_summary %>%
  mutate(Fraction = Frequency / sum(Frequency))

#P. lutea unique ITS2 types
lut_ITS_summary <- as.data.frame(table(lut_only_ITS$ITS2_type))
colnames(lut_ITS_summary) <- c("ITS2_type", "Frequency")
lut_ITS_summary <- lut_ITS_summary %>%
  mutate(Fraction = Frequency / sum(Frequency))

#P. lobata unique ITS2 types
lob_ITS_summary <- as.data.frame(table(lob_only_ITS$ITS2_type))
colnames(lob_ITS_summary) <- c("ITS2_type", "Frequency")
lob_ITS_summary <- lob_ITS_summary %>%
  mutate(Fraction = Frequency / sum(Frequency))

#P. lutea all ITS2 types
lut_ITS_summary_all <- as.data.frame(table(lut_ITS_all$ITS2_type))
colnames(lut_ITS_summary_all) <- c("ITS2_type", "Frequency")

lut_ITS_summary_all <- lut_ITS_summary_all %>% 
  mutate(
    Fraction = Frequency / sum(Frequency))

#P. lobata all ITS2 types
lob_ITS_summary_all <- as.data.frame(table(lob_ITS_all$ITS2_type))
colnames(lob_ITS_summary_all) <- c("ITS2_type", "Frequency")

lob_ITS_summary_all <- lob_ITS_summary_all %>% 
  mutate(
    Fraction = Frequency / sum(Frequency))


## Pie charts

#shared ITS2 Pie Chart
pie_lutlob <- ggplot(lutlob_ITS_summary, aes(x = "", y = Frequency, fill = ITS2_type)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar(theta = "y") +
  labs(title = "Both host species",
       fill = expression("Taxa (" * italic("Cladocopium") * " ASVs)")) +
  theme_void() +
  scale_fill_manual(values = ASV_colors) +
  theme(legend.position = "right")

#without labels
pie_lutlob <- ggplot(lutlob_ITS_summary, aes(x = "", y = Frequency, fill = ITS2_type)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar(theta = "y") +
  theme_void() +
  scale_fill_manual(values = ASV_colors) +
  theme(
    legend.position = "none", #remove legend
    plot.title = element_blank() #remove title
  )

#P. lutea unique ITS2 Pie Chart
pie_lut <- ggplot(lut_ITS_summary, aes(x = "", y = Frequency, fill = ITS2_type)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar(theta = "y") +
  labs(title = "P. lutea",
       fill = expression("Taxa (" * italic("Cladocopium") * " ASVs)")) +
  theme_void() +
  scale_fill_manual(values = ASV_colors) +
  theme(legend.position = "right")

#without labels
pie_lut <- ggplot(lut_ITS_summary, aes(x = "", y = Frequency, fill = ITS2_type)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar(theta = "y") +
  theme_void() +
  scale_fill_manual(values = ASV_colors) +
  theme(
    legend.position = "none", #remove legend
    plot.title = element_blank() #remove title
  )

#P. lobata unique ITS2 Pie Chart
pie_lob <- ggplot(lob_ITS_summary, aes(x = "", y = Frequency, fill = ITS2_type)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar(theta = "y") +
  labs(title = "P. lobata",
       fill = expression("Taxa (" * italic("Cladocopium") * " ASVs)")) +
  theme_void() +
  scale_fill_manual(values = ASV_colors) +
  theme(legend.position = "right")

#without labels
pie_lob <- ggplot(lob_ITS_summary, aes(x = "", y = Frequency, fill = ITS2_type)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar(theta = "y") +
  theme_void() +
  scale_fill_manual(values = ASV_colors) +
  theme(
    legend.position = "none", #remove legend
    plot.title = element_blank() #remove title
  )

#P. lutea all ITS2 Pie Chart
pie_lut_all <- ggplot(lut_ITS_summary_all, aes(x = "", y = Frequency, fill = ITS2_type)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar(theta = "y") +
  labs(title = "P. lutea all",
       fill = expression("Taxa (" * italic("Cladocopium") * " ASVs)")) +
  theme_void() +
  scale_fill_manual(values = ASV_colors) +
  theme(legend.position = "right")

#P. lobata all ITS2 Pie Chart
pie_lob_all <- ggplot(lob_ITS_summary_all, aes(x = "", y = Frequency, fill = ITS2_type)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar(theta = "y") +
  labs(title = "P. lobata all",
       fill = expression("Taxa (" * italic("Cladocopium") * " ASVs)")) +
  theme_void() +
  scale_fill_manual(values = ASV_colors) +
  theme(legend.position = "right")


## Create one big happy VD plot with pie charts

library(cowplot)

#create Euler/Venn Diagram for shared ASVs
VD_Spp <- euler(c("P. lutea" = (lut_only), "P. lobata" = (lob_only), "P. lutea&P. lobata" = (lutlob_length))) 

VD_Spp_plot <- plot(VD_Spp, quantities = TRUE, font=1, cex=1, alpha=0.5,
                    fill=species_colors,col=species_colors,border=species_colors,lwd=c(2,2,2))

#combine into a single figure
final_VDspp_plot <- ggdraw() +
  draw_plot(VD_Spp_plot, x = 0, y = 0, width = 1, height = 1) +
  draw_plot(pie_lob, x = 0.08, y = 0.08, width = 0.25, height = 0.25) +
  draw_plot(pie_lut, x = 0.75, y = 0.08, width = 0.25, height = 0.25) +
  draw_plot(pie_lutlob, x = 0.50, y = 0.55, width = 0.25, height = 0.25) #+
#draw_label("(a)", x = 0.05, y = 1, hjust = 0, vjust = 1, fontface = "bold", size = 18)

#save the final figure
ggsave("VD_Spp_with_pies.png", plot = final_VDspp_plot, width = 8, height = 8, units = "in")


########### Prevalence and abundance per host species ##############

#create empty lists to store results for each site
prev_otu_list <- list()
otu_per_sample_list <- list()
prev_otu2_list <- list()
prev_otu4_list <- list()


#define site names and corresponding phyloseq objects
species_list <- list(
  "P.lobata" = new_phylo_lob,
  "P.lutea" = new_phylo_lut
)

for (species in names(species_list)) {
  
  phylo_obj <- species_list[[species]]  #get phyloseq object for the site
  
  #extract ASV table
  otu_table_data <- t(otu_table(phylo_obj))
  otu_df <- as.data.frame(otu_table_data)
  otu_df$OTU <- rownames(otu_df)
  
  #compute prevalence of each feature
  prevdf <- apply(X = otu_table(phylo_obj), 
                  MARGIN = ifelse(taxa_are_rows(phylo_obj), yes = 1, no = 2),
                  FUN = function(x) sum(x > 0))
  
  #add taxonomy and total read counts
  prevdf <- data.frame(Prevalence = prevdf,
                       TotalAbundance = taxa_sums(phylo_obj),
                       tax_table(phylo_obj)) %>%
    rownames_to_column() %>%
    dplyr::rename(OTU = 'rowname')
  
  #merge prevalence data with OTU table
  prev_otu <- prevdf %>%
    left_join(otu_df, by = "OTU") %>%
    tidyr::pivot_longer(cols = 5:ncol(.), names_to = "Sample", values_to = "ASV_count_per_sample") %>%
    filter(ASV_count_per_sample != 0) %>%
    group_by(Sample) %>%
    mutate(Total_reads_per_sample = sum(ASV_count_per_sample),
           ITS2_diversity_per_sample = n_distinct(ITS2_type))
  
  #store prev_otu in the list
  prev_otu_list[[species]] <- prev_otu
  
  #compute total OTUs per sample 
  otu_per_sample <- prev_otu %>%
    group_by(Sample) %>%
    summarise(Total_OTUs = n_distinct(OTU))
  
  #store OTU count per sample in the list
  otu_per_sample_list[[species]] <- otu_per_sample
  
  
  #compute abundance of each ITS2 type within each sample
  prev_otu2 <- prev_otu %>%
    group_by(Sample, ITS2_type) %>%
    mutate(Total_ITS2_count_per_sample = sum(ASV_count_per_sample),
           abundance_per_sample = Total_ITS2_count_per_sample / Total_reads_per_sample) %>%
    summarize(abundance_per_sample = mean(abundance_per_sample), .groups = "keep") %>%
    ungroup()
  
  #store prev_otu2 in the list
  prev_otu2_list[[species]] <- prev_otu2
  
  
  #compute relative abundance across samples
  prev_otu3 <- prev_otu2 %>%
    group_by(ITS2_type) %>%
    summarize(mean_abund = mean_se(abundance_per_sample)) %>%
    mutate(se_abund = mean_abund$y - mean_abund$ymin) %>%
    arrange(desc(mean_abund))
  
  #extract prevalence of ITS2 types and combine with relative abundance
  prev_otu4 <- prev_otu %>%
    group_by(ITS2_type) %>%
    mutate(Prevalence = n_distinct(Sample)) %>%
    summarize(Prevalence = max(Prevalence)) %>%
    left_join(prev_otu3, by = "ITS2_type") %>%
    arrange(desc(Prevalence))
  
  #store prev_otu4 in the list
  prev_otu4_list[[species]] <- prev_otu4
  
  #save results
  file_name <- paste0("vst_prevalence_abundance_", species, ".csv")
  write.csv(prev_otu4, file = file_name)
  
  print(paste("Processed species:", species, " -> Saved as:", file_name))
}


##example: To view OTU counts for MUI
otu_per_sample_list[["P.lobata"]]

#or to see all
otu_per_sample_list

#to access data for individual sites from the lists above (MUI prev_otu2 as example)
prev_otu2_list[["P.lobata"]]

prev_otu2_list[["P.lutea"]]

prev_samp_lut <- prev_otu2_list[["P.lutea"]]


#prevalence of Specific ASVs
prev_otu_lut <- prev_otu_list[["P.lutea"]]

prev_otu_lot <- prev_otu_list[["P.lobata"]]

#count the number of unique ASVs per sample
otu_per_sample_lut <- prev_otu_lut %>%
  group_by(Sample) %>%
  summarise(Total_OTUs = n_distinct(OTU))

