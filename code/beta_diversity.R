library(vegan)
library(ape)

source("code/load_data.R")

permutations = 9999

generate_dist <- function(distance = "bray"){
  
  table_t <- table %>% 
    column_to_rownames(var = "featureid") %>% 
    t()
  
  table.transformed <- table_t^1/4
  
  dist <- vegan::vegdist(table.transformed, method= distance)
}

dist <- generate_dist(distance = "bray")

generate_pcoa <- function(dist = "dist") {
  pcoa <- pcoa(dist, correction = "cailliez") 
  
  p_var <- (pcoa$values$Eigenvalues/pcoa$trace)*100
  
  pcoa_vectors <- pcoa$vectors %>% 
    as_tibble(rownames = "sampleid") %>% 
    select(sampleid, Axis.1,Axis.2)
  colnames(pcoa_vectors) <- c("sampleid", "PCo1", "PCo2")
  
  pcoa.metadata <- inner_join(pcoa_vectors, metadata, by = "sampleid") 

  adonis.output <- adonis(dist ~ treatment,
                          permutations = permutations,
                          data = metadata)

  f_value <- round(adonis.output$aov.tab$F.Model[1],2)

  p_value <- ifelse(adonis.output$aov.tab$`Pr(>F)`[1] < 0.001,
                    "< 0.001",
                    adonis.output$aov.tab$`Pr(>F)`[1])

  output <- list(pcoa.metadata, p_var, f_value, p_value)
  return(output)
}



