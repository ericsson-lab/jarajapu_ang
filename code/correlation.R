library(tidyverse)
library(readxl)
library(scales)

# Pull in metadata, table, taxa_table
source("./code/load_data.R")

# Load values to correlate
output_metrics <- read_excel("data/metabolites.xlsx") %>%  
  column_to_rownames(var = "outcome") %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "outcome") %>% 
  as_tibble() %>% 
  drop_na(`FITC-dextran`) %>% 
  column_to_rownames(var = "outcome")
colnames(output_metrics) <- c("Zonulin", "FITC")

metadata <- metadata %>% arrange(desc(age))

drop_samples <- c("JY42", "YP010")

# Generate table at specified taxonomic level
filter_table <- function(taxa_level = "Famly") {
  table %>% 
    pivot_longer(-featureid, names_to = "sampleid", values_to = "count") %>% 
    filter(!sampleid %in% drop_samples) %>%
    left_join(., taxonomy, by = "featureid") %>% 
    filter(level == taxa_level) %>% 
    group_by(sampleid,taxon) %>% 
    summarise(count = sum(count)) %>% 
    tidyr::pivot_wider(names_from = "taxon", values_from = "count", values_fn = sum) %>% 
    column_to_rownames(var = "sampleid") %>% 
    select_if(colSums(.) != 0) # Removes columns with 0 counts
}

# Generate family table
family_table <- filter_table(taxa_level = "Family")

# Filter down table, default is taxa must be > 0.1%
# Pulls taxa names
filter_table_names <- function(input_table, filter_level = 0.001) {
  input_table %>% 
    rownames_to_column(var = "sampleid") %>% 
    pivot_longer(-sampleid, names_to = "taxon", values_to = "count") %>% 
    group_by(sampleid) %>% 
    mutate(rel_abund = count/sum(count)) %>% 
    group_by(taxon) %>% 
    summarize(mean_rel_abund = mean(rel_abund)) %>% 
    arrange(desc(mean_rel_abund)) %>% 
    filter(mean_rel_abund > filter_level) %>%
    pull(taxon)
}

family_table_filtered_names <- filter_table_names(family_table, 0) 

# select genera with more than 0.01% rel_abund across all samples
family_0.1p <- family_table %>% select(all_of(family_table_filtered_names))

# spearman correlation of 16S at genera level and metabolites
res <- cor(output_metrics, family_0.1p, method = "spearman")

## Attempting to perform significance testing and plotting only sig correlations
## Following Associations between the gut microbiome and metabolome in early life" by Nguyen et al. 2021+:
## https://github.com/qpmnguyen/infant_metabolome_microbiome/blob/master/R/analysis_spearman_corr.R

# pairwise testing correlation
# res <- cor(genus_0.01p, metabolite, method = "spearman")
p_mat <- t(apply(family_0.1p, 2, function(x){
  p_vals <- c()
  print(length(x))
  for (i in 1:ncol(output_metrics)){
    p_vals[i] <- cor.test(x = x, y = output_metrics[,i], method = "spearman")$p.value
  }
  return(p_vals)
}))
# p-value adjustment using BH
adj_mat <- matrix(p.adjust(as.vector(p_mat), method = "BH"),
                  ncol = ncol(p_mat), nrow = nrow(p_mat), byrow = F)

output <- list(
  cor_mat = res,
  p_mat = p_mat,
  bh_p_mat = adj_mat)


colnames(output$p_mat) <-  rownames(output$cor_mat)
rownames(output$p_mat) <- colnames(output$cor_mat)

colnames(output$bh_p_mat) <-  rownames(output$cor_mat)
rownames(output$bh_p_mat) <- colnames(output$cor_mat)

cor_save <- output$cor_mat %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "metabolite") %>% 
  pivot_longer(-metabolite, names_to = "taxon", values_to = "rho") %>% 
  mutate(id = paste0(taxon, ":", metabolite)) %>% 
  select(id, rho)

p_save <- output$p_mat %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "taxon") %>% 
  pivot_longer(-taxon, names_to = "metabolite", values_to = "pvalue") %>% 
  mutate(id = paste0(taxon, ":", metabolite)) %>% 
  select(id, pvalue)

bh_p_save <- output$bh_p_mat %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "genus") %>% 
  pivot_longer(-genus, names_to = "metabolite", values_to = "pvalue_bh") %>% 
  mutate(id = paste0(genus, ":", metabolite)) %>% 
  select(id, pvalue_bh)

comb <- inner_join(cor_save, bh_p_save) %>% 
  inner_join(p_save) %>% 
  arrange(id) %>% 
  separate(col = id, into = c("taxon", "output"), sep = ":") %>% 
  mutate(sig = case_when(pvalue_bh < 0.05 ~ "*",
                         pvalue_bh < 0.01 ~ "**",
                         pvalue_bh < 0.001 ~ "***"))

sig_list<- comb %>%
  drop_na(sig) %>% 
  pull(taxon) %>% 
  unique()


  