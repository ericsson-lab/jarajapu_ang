library(vegan)
library(microbiome)

source("./code/load_data.R")

alpha.stats <- table %>% 
  pivot_longer(-featureid, 
               names_to = "sampleid",
               values_to = "count") %>% 
  group_by(sampleid) %>% 
  summarize(richness = richness(count, index = "observed"),
            shannon = vegan::diversity(count, 
                                       index = "shannon",
                                       base = exp(1)),
            chao1 = richness(count, index = c("chao1")) %>% as.double(),
            simpson = vegan::diversity(count, 
                                       index = "simpson")) %>% 
  inner_join(., metadata, by = "sampleid")




