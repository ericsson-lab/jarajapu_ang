library(ggpubr)
library(rstatix)
library(patchwork)
library(glue)
library(EcolUtils)
library(gt)

source("./code/beta_diversity.R")
source("./code/alpha_stats.R")
source("./code/correlation.R")

labels = c("Young", "Old", "Old\nAng(1-7)")

alpha.stats$group <- factor(alpha.stats$group, levels = c("young - control", "old - control", "old - Ang-(1-7)"))

chao1_stat<- alpha.stats %>% 
  tukey_hsd(chao1 ~ group) %>% 
  add_xy_position()
shannon_stat <- alpha.stats %>% 
  tukey_hsd(shannon~group) %>% 
  add_xy_position()

alpha_plot <- function(alpha_stat, stat, y_lab){
  
  alpha.stats %>% 
    ggplot(aes(x = group, 
               y = {{alpha_stat}})) +
    geom_boxplot(width = 0.4, size = 1) +
    geom_point(aes(fill = group),
               shape = 21, size = 2.5, color = "black")+
    scale_fill_manual(values = c("white", "black", "gray"),
                      breaks = c("young - control", "old - control", "old - Ang-(1-7)")) +
    scale_x_discrete(labels = labels) +

    theme_classic(base_size = 12) +
    labs(y = y_lab) +
    theme(
      axis.text = element_text(face = "bold", color = "black"),
      axis.title = element_text(face = "bold", color = "black"),
      axis.title.x = element_blank(),
      
      legend.position = "none",
      
    ) +
    
    stat_pvalue_manual(data = stat, hide.ns = T,
                       label.size = 8, tip.length = 0,
                       bracket.size = 2, linetype = 1)
  
}

chao1_plot <- alpha_plot(alpha_stat = chao1, stat = chao1_stat, y_lab = "Chao1 Index") +    
  scale_y_continuous(limits = c(0, 800), expand = c(0,0)) 
shannon_plot <- alpha_plot(alpha_stat = shannon, stat = shannon_stat, y_lab = "Shannon Index") +
  scale_y_continuous(limits = c(2, 6), expand = c(0,0)) 


bc_pcoa <- generate_pcoa(dist)

bc_plot <- bc_pcoa[1] %>% 
  as.data.frame() %>% 
  ggplot(aes(x = PCo1, y = PCo2, group = group)) +
  geom_point(aes(fill = group), size = 3, shape = 21, color = "black") +
  stat_ellipse() +
  
  scale_fill_manual(values = c("white", "black", "gray"),
                    breaks = c("young - control", "old - control", "old - Ang-(1-7)"),
                    labels = labels) +
  
  labs(x = paste0("PCo1 - ", round(bc_pcoa[[2]][1], 2), "%"),
       y = paste0("PCo2 - ", round(bc_pcoa[[2]][2], 2), "%"),
       caption = glue("F: {bc_pcoa[3]}, p = {bc_pcoa[4]}")) +
  
  scale_y_continuous(limits = c(-.55, 0.75)) +
  scale_x_continuous(limits = c(-0.75, 0.85)) +
  theme_classic(base_size = 12) +
  theme(
    axis.text = element_text(face = "bold", color = "black"),
    axis.title = element_text(face = "bold", color = "black"),
    
    legend.title = element_blank(),
    legend.text = element_text(face = "bold", color = "black", hjust = 0.5),
    legend.position = c(0.85, 0.8),
    
    plot.caption = element_text(face = "bold", size = 12)
  ) 


cor_plot <- comb %>% 
  filter(taxon %in% sig_list) %>% 
  ggplot(aes(x = output, fill= rho, label = sig)) +
  geom_tile(aes(y = taxon), color = NA) +
  geom_text(aes(label = sig, y = taxon), color = "white", size = 8, nudge_y = -0.1) +
  scale_fill_gradient2(low = muted("blue"), high = muted("red"), limits = c(-1, 1) ) +
  scale_y_discrete(expand = c(0,0)) +
  scale_x_discrete(expand = c(0,0)) +
  theme_minimal() +
  theme(
    axis.text = element_text(face = "bold", color = "black"),
    # axis.text.x = element_text(angle = 25, hjust = 1),
    axis.title = element_blank(),
    legend.text = element_text(face = "bold", color = "black"),
    legend.title = element_text(face = "bold.italic", color = "black", hjust = 0.5)
  ) +
  
  guides(fill=guide_colorbar(ticks.colour = NA, title = "Rho"))
cor_plot
#   diversity <- (chao1_plot + shannon_plot) / bc_plot
# 
# (diversity | plot_spacer() | cor_plot) + plot_layout(width = c(2,0.1,0.5)) &
#   plot_annotation(tag_levels = "A") + patchwork::align_patches()
# 
# 
# 
#   plot_annotation(tag_levels = "A") &
#   theme(
#     plot.tag = element_text(face = "bold")
#   )

  
library(cowplot)
  
alpha <- plot_grid(chao1_plot, shannon_plot, labels = c("A", "B"), hjust = 0.75, vjust = 1, label_size = 16)
beta <- plot_grid(alpha, bc_plot, nrow = 2, labels = c(NA, "C"), hjust = 0.75, vjust = 1, label_size = 16)
cor <- plot_grid(beta, cor_plot, rel_widths = c(0.5, 0.5), labels = c(NA, "D"), hjust = 0.75, vjust = 1, label_size = 16) + 
  theme(plot.margin = unit(c(0.4, 0.4, 0.4, 0.4), "cm"))
ggsave(plot =  cor, 
       "plots/microbiome_plot.png",
       width = 7.5,
       height = 6,
       units = c("in"),
       bg = "white")

 dist

adonis.pair(dist, factor(metadata$group)) %>% 
  gt() %>% 
  gtsave("plots/pairwise_adonis.png")



