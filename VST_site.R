# load packages
library(phyloseq); packageVersion("phyloseq") #'1.46.0'
library(ggplot2); packageVersion("ggplot2") #'3.5.1'
library(tidyverse); packageVersion("tidyverse") #‘2.0.0’


#using vst normalised count data

vst_phylo.percent <- readRDS("vst_phylo.rds")

############################# Subset by Site ###########################
new_phylo_mui <- subset_samples(vst_phylo.percent,Site=="MUI")
new_phylo_mui <- prune_taxa(taxa_sums(new_phylo_mui) > 0, new_phylo_mui)

new_phylo_vlf <- subset_samples(vst_phylo.percent,Site=="VLF")
new_phylo_vlf <- prune_taxa(taxa_sums(new_phylo_vlf) > 0, new_phylo_vlf)

new_phylo_bun <- subset_samples(vst_phylo.percent,Site=="BUN")
new_phylo_bun <- prune_taxa(taxa_sums(new_phylo_bun) > 0, new_phylo_bun)

new_phylo_tan <- subset_samples(vst_phylo.percent,Site=="TAN")
new_phylo_tan <- prune_taxa(taxa_sums(new_phylo_tan) > 0, new_phylo_tan)

new_phylo_lktr <- subset_samples(vst_phylo.percent,Site=="LK/TR")
new_phylo_lktr <- prune_taxa(taxa_sums(new_phylo_lktr) > 0, new_phylo_lktr)

new_phylo_osp <- subset_samples(vst_phylo.percent,Site=="OSP")
new_phylo_osp <- prune_taxa(taxa_sums(new_phylo_osp) > 0, new_phylo_osp)

new_phylo_cor <- subset_samples(vst_phylo.percent,Site=="COR")
new_phylo_cor <- prune_taxa(taxa_sums(new_phylo_cor) > 0, new_phylo_cor)



#create empty lists to store results for each site
prev_otu_list <- list()
otu_per_sample_list <- list()
prev_otu2_list <- list()
prev_otu4_list <- list()


#define site names and corresponding phyloseq objects
site_list <- list(
  "MUI" = new_phylo_mui,
  "VLF" = new_phylo_vlf,
  "BUN" = new_phylo_bun,
  "TAN" = new_phylo_tan,
  "LKTR" = new_phylo_lktr,
  "OSP" = new_phylo_osp,
  "COR" = new_phylo_cor
)

for (site in names(site_list)) {
  
  phylo_obj <- site_list[[site]]  #get phyloseq object for the site
  
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
  prev_otu_list[[site]] <- prev_otu
  
  #compute total ASVs per sample 
  otu_per_sample <- prev_otu %>%
    group_by(Sample) %>%
    summarise(Total_OTUs = n_distinct(OTU))
  
  #store ASV count per sample in the list
  otu_per_sample_list[[site]] <- otu_per_sample
  
  
  #compute abundance of each ITS2 type within each sample
  prev_otu2 <- prev_otu %>%
    group_by(Sample, ITS2_type) %>%
    mutate(Total_ITS2_count_per_sample = sum(ASV_count_per_sample),
           abundance_per_sample = Total_ITS2_count_per_sample / Total_reads_per_sample) %>%
    summarize(abundance_per_sample = mean(abundance_per_sample), .groups = "keep") %>%
    ungroup()
  
  #store prev_otu2 in the list
  prev_otu2_list[[site]] <- prev_otu2
  
  
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
  prev_otu4_list[[site]] <- prev_otu4
  
  #save results
  file_name <- paste0("vst_prevalence_abundance_", site, ".csv")
  write.csv(prev_otu4, file = file_name)
  
  print(paste("Processed site:", site, " -> Saved as:", file_name))
}

#example: To view ASV counts for MUI
otu_per_sample_list[["MUI"]]

#or to see all
otu_per_sample_list

#to access data for individual sites from the lists above (MUI prev_otu2 as example)
prev_otu2_list[["MUI"]]

prev_samp_mui <- prev_otu2_list[["MUI"]]
prev_samp_vlf <- prev_otu2_list[["VLF"]]
prev_samp_bun <- prev_otu2_list[["BUN"]]
prev_samp_tan <- prev_otu2_list[["TAN"]]
prev_samp_lktr <- prev_otu2_list[["LKTR"]]
prev_samp_osp <- prev_otu2_list[["OSP"]]
prev_samp_cor <- prev_otu2_list[["COR"]]


#prevalence of specific ASVs
prev_otu_m <- prev_otu_list[["MUI"]]

#count the number of unique ASVs per sample
otu_per_sample_m <- prev_otu_m %>%
  group_by(Sample) %>%
  summarise(Total_OTUs = n_distinct(OTU))


####### Finding out proportion and frequency of taxa per site ######

#loop through each site and compute the proportion and frequency of ITS2 types
for (site in names(site_list)) {
  
  #access the phyloseq object for the current site
  phylo_obj <- site_list[[site]]
  
  #extract ITS2 type data for the current site
  ITS2_data <- data.frame(ITS2_type = data.frame(tax_table(phylo_obj))$ITS2_type)
  
  #create a summary of ITS2 types (frequency and proportion)
  ITS_summary <- as.data.frame(table(ITS2_data$ITS2_type))
  colnames(ITS_summary) <- c("ITS2_type", "Frequency")
  
  #add the fraction (proportion) of each ITS2 type
  ITS_summary <- ITS_summary %>% mutate(Fraction = Frequency / sum(Frequency))
  
  #store the result in a list
  assign(paste0(site, "_ITS_summary_all"), ITS_summary)
  
  #print a message for each site processed
  print(paste("Processed", site, "ITS2 type summary"))
}

