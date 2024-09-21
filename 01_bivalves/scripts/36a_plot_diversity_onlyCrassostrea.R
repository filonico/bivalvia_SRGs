library(tidyr)
library(dplyr)
library(ggplot2)
library(ggtree)


######################
#     READ INPUT     #
######################

# amino acid diversity file
diversity_file <- "13a_distribution_divergence_onlyCrassostrea/distance_median_values.tsv"


############################
#     OUTPUT FILENAMES     #
############################

distance_median_values_quant <- "13a_distribution_divergence_onlyCrassostrea/distance_median_values_quant.tsv"

diversity_panel <- "13a_distribution_divergence_onlyCrassostrea/02_plot_diversity/supp_fig_S15.png"


##################################
#     PLOT DISTRIBUTION DATA     #
##################################

# load diversity data
diversity_data <- read.table(diversity_file, header = TRUE, sep = "\t")

# compute density
density.median <- density(diversity_data$median)
density.median.df <- data.frame(x = density.median$x, y = density.median$y)

# compute quantiles
quantile = c(0.95, 0.99)
density.median.quant <- quantile(diversity_data$median, prob = quantile)

density.median.df$quant <- factor(findInterval(density.median.df$x, density.median.quant))
diversity_data$quant <- factor(findInterval(diversity_data$median, density.median.quant))

# plot density
number_of_genes <- nrow(diversity_data)

plot_density <- ggplot(density.median.df, aes(x = x, y = y)) +
  
  geom_ribbon(aes(ymin = 0, ymax = y, fill = quant)) +
  geom_line(aes(color = as.factor(quant)), linewidth = 1.3, lineend = "round", show.legend = FALSE) +
  
  coord_flip() +
  
  scale_colour_manual(values = c("#0e65d1", "#b353df", "#cf1c47"),
                      labels = c("< 95%", paste0(quantile*100, "%"))) +
  
  scale_fill_manual(values = c("#7eaee8", "#d79ef1", "#ee809a"),
                    labels = c("< 95%", paste0(quantile*100, "%"))) +  
  
  guides(fill = guide_legend(title = "Quantiles")) +
  
  scale_x_continuous(limits = c(0, 3), breaks = seq(0, 3, 1)) +
  
  scale_y_reverse(limits = c(1,0)) +
  
  xlab("Amino acid divergence") +
  ylab("Density") +
  
  ggtitle(paste0("Distribution out of\n", number_of_genes, " genes")) +
  
  theme_minimal() +
  
  theme(panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_blank(),
        # plot.borde = element_blank(),
        # panel.background = element_rect(fill = 'transparent'),
        # plot.background = element_rect(fill = 'transparent', color = NA),
        axis.line = element_line(linewidth = 0.5),
        axis.text = element_text(size = 8),
        legend.background = element_rect(fill = 'transparent', color = NA))

plot_density

##########################################################
#     CREATE THE DATAFRAME TO PLOT DIVERSITY OF SRGs     #
##########################################################

# rename genes
diversity_data_SRGs <- diversity_data %>%
  filter(stringr::str_detect(gene, stringr::regex("dmrt|sox|fox", ignore_case = TRUE))) %>%
  dplyr::mutate(gene = stringr::str_extract(gene, "[^/]+$")) %>%
  dplyr::mutate(gene = stringr::str_extract(gene, "^[^_]+_[^_]+_[^_]+"))

diversity_data_SRGs

diversity_data_SRGs$quant <- as.factor(diversity_data_SRGs$quant)
diversity_data_SRGs$median <- as.numeric(diversity_data_SRGs$median)


#######################################
#     PLOT DIVERSITY DATA OF SRGs     #
#######################################

# plot median diversity (remove one Fox-B1/2, which is a duplicate)
plot_points <- diversity_data_SRGs %>%
  # filter(stringr::str_detect(gene_names, "Sox")) %>%
  
  ggplot(aes(y = median, x = gene, colour = quant, fill = quant)) +
  
  geom_segment(aes(y = -Inf, x = gene, yend = median, xend = gene),
               col = "#706d73", alpha = 0.3, linewidth = 0.8) +
  
  geom_point(size = 4, shape = 21, stroke = 1.5) +
  
  geom_text(aes(y = median, x = gene, label = gene),
            col = "#706d73", size = 3, fontface = "italic",
            angle = 90, hjust = 0, nudge_y = 0.1) +
  
  coord_cartesian(clip = "off") +
  
  xlab("Genes") +
  
  scale_y_continuous(limits = c(0, 3), breaks = seq(0, 3, 1)) +
  # scale_x_discrete(breaks = gene_order[-c(5,12)]) +
  
  # scale_colour_manual(values = c("#fcba9b", "#f57b83", "#e91674"), na.value = "transparent",
  #                     labels = c("<90%", paste0(quantile*100, "%")),
  #                     guide = "none") +
  # 
  # scale_fill_manual(values = c("#ffdac9", "#f4999f", "#ec4690"),
  #                   labels = c("< 90%", paste0(quantile*100, "%")),
  #                   guide = "none") +
  
  scale_colour_manual(values = c("#0e65d1", "#b353df", "#cf1c47"),
                      labels = c("< 95%", paste0(quantile*100, "%")),
                      guide = "none") +
  
  scale_fill_manual(values = c("#7eaee8", "#d79ef1", "#ee809a"),
                    labels = c("< 95%", paste0(quantile*100, "%")),
                    guide = "none") +
  
  scale_size_continuous(
                        # breaks = c(min(joined_data_SRGs$species, na.rm = TRUE), max(joined_data_SRGs$species, na.rm = TRUE)),
                        name = "Quantiles Number of species") +
  
  # guides(colour = guide_legend(title = "Quantiles")) +
  
  theme_minimal() +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.line.y = element_blank(),
        axis.text.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.y = element_blank(),
        axis.line = element_line(size = 0.5))

plot_points


######################
#     PLOT PANEL     #
######################


final_panel <- ggpubr::ggarrange(plot_density, plot_points,
                                 nrow = 1, ncol = 2,
                                 align = "hv",
                                 common.legend = TRUE,
                                 legend = "top",
                                 widths = c(0.3,1),
                                 labels = c("A", "B"))
                 
final_panel


########################
#     SAVE OUTPUTS     #
########################

write.table(diversity_data, distance_median_values_quant, quote = FALSE, sep = "\t", row.names = FALSE)

ggsave(diversity_panel,
       plot = final_panel, device = "png",
       dpi = 300, height = 7, width = 9, units = ("in"), bg = 'white')
