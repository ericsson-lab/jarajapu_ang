library(tidyverse)
library(readxl)

metadata <- read_tsv("data/metadata.txt")

raw_table <- read_excel("data/210928_results.xlsx", 
                        sheet = "annotated_table",
                        skip = 1) %>% 
  rename(featureid = `#OTU ID`)
  

table_taxa <- read_excel("data/210928_results.xlsx", 
                         sheet = "annotated_table",
                         skip = 1) %>% 
  rename(featureid = `#OTU ID`) %>% 
  select(featureid, metadata$sampleid, Taxon)

table <- table_taxa %>% 
  select(-Taxon)

rarefy_depth <- table %>% 
  pivot_longer(-featureid) %>% 
  group_by(name) %>% 
  summarise(count = sum(value)) %>% 
  arrange(count) %>% 
  pull(count[1]) %>% 
  magrittr::extract2(1) -1

# rarefied_table <- table %>%
#   column_to_rownames(var = "featureid") %>%
#   t() %>%
#   rrarefy(., rarefy_depth) %>%
#   t() %>%
#   as.data.frame() %>%
#   rownames_to_column(var = "featureid") %>%
#   write_tsv("data/rarefied_table_91951.tsv")

table <- read_tsv("data/rarefied_table_91951.tsv")

taxonomy <- table_taxa %>% 
  select(featureid, Taxon) %>% 
  separate(Taxon, 
           into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus"),
           sep = "; ") %>% 
  pivot_longer(-featureid,
               names_to = "level",
               values_to = "taxon") %>% 
  group_by(featureid) %>% 
  mutate(taxon = str_sub(taxon, 4)) %>% 
  fill(taxon)



