library(tidyr)
library(dplyr)
library(ggplot2)
library(ggtree)


######################
#     READ INPUT     #
######################

# gene trees
dmrtTree_data <- read.tree(text = "(Dsx,(Dmrt-99B,(Dmrt-93B,Dmrt-11E)));")
soxTree_data <- read.tree(text = "(Sox-15,(Sox-100B,(Sox-14,(Dichaete,(Sox-21B,(Sox-N,Sox-21A))))));")
foxTree_data <- read.tree(text = "(Fox-3F,((Ches-1,Jumeau),((Hcm-1,(Fox-P,Fox-O)),(Fd-102C,(Biniou,(Fd-2/Fox-L1,((Fd-19B,(Slp-1,Slp-2)),(Fd-3/Fd-59A,(Crocodile,(Fkh,(Fd-5/Fd-96Cb,Fd-4/Fd-96Ca)))))))))));")

# amino acid diversity file
diversity_file <- "13_distribution_divergence/distance_median_values.tsv"

# gene names conversion tables
file_paths <- list.files(path = "./10_SRG_decomposition/", pattern = "*conversion.tsv", full.names = TRUE)

# occurrence table
occurrence_data <- read.table("06_possvm_orthology/01_plot_occurrence/occurrence_matrix_withOutgroups.tsv",
                              header = TRUE, sep = "\t", check.names=FALSE)


############################
#     OUTPUT FILENAMES     #
############################

distance_median_values_quant <- "13_distribution_divergence/distance_median_values_quant.tsv"

diversity_panel <- "13_distribution_divergence/02_plot_diversity/diversity_panel.pdf"


##########################
#     PLOT GENE TREE     #
##########################

# load gene trees
dmrtTree_plot <- ggtree(dmrtTree_data) +
  # geom_tiplab() +
  coord_flip() +
  ggpubr::theme_transparent() +
  theme(plot.margin = unit(c(0, 0, 0, 0), "points"))

dmrtTree_plot

soxTree_plot <- ggtree(soxTree_data) +
  # geom_tiplab() +
  coord_flip() +
  ggpubr::theme_transparent() +
  theme(plot.margin = unit(c(0, 0, 0, 0), "points"))

soxTree_plot

foxTree_plot <- ggtree(foxTree_data) +
  # geom_tiplab() +
  coord_flip() +
  ggpubr::theme_transparent() +
  theme(plot.margin = unit(c(0, 0, 0, 0), "points"))

foxTree_plot

panel_geneTrees <- ggpubr::ggarrange(dmrtTree_plot, soxTree_plot, foxTree_plot,
                                     ncol = 3, nrow = 1,
                                     align = "hv",
                                     widths = c(0.5, 0.7, 1.6))

panel_geneTrees


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
  geom_line(aes(color = as.factor(quant)), linewidth = 1.5, lineend = "round", show.legend = FALSE) +
  
  coord_flip() +
  
  scale_colour_manual(values = c("#0e65d1", "#b353df", "#cf1c47"),
                      labels = c("< 95%", paste0(quantile*100, "%"))) +
  
  scale_fill_manual(values = c("#7eaee8", "#d79ef1", "#ee809a"),
                    labels = c("< 95%", paste0(quantile*100, "%"))) +
  
  guides(fill = guide_legend(title = "Quantiles")) +
  
  scale_x_continuous(limits = c(-0.11, 3), breaks = seq(0, 3, 0.5)) +
  
  scale_y_reverse(limits = c(2.1,0)) +
  
  xlab("Amino acid divergence") +
  ylab("Density") +
  
  ggtitle(paste0("Distribution out of\n", format(round(as.numeric(number_of_genes), " genes")) +
  
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

data_list <- lapply(file_paths, function(file) {
  read.table(file, header = FALSE, sep = "\t")
})

gene_conversion <- do.call(rbind, data_list)

# rename conversion tables
gene_conversion <- gene_conversion %>%
  mutate(V1 = stringr::str_extract(V1, "[^.]+"),
         V2 = stringr::str_extract(V2, "[^/]+"),
         V2 = stringr::str_replace_all(V2, c("dmrt" = "Dmrt", "dsx" = "Dsx", "OG10_Dmrt-A2" = "OG10_Dmrt-93B", "OG1_Dmrt-A2" = "OG1_Dmrt-99B", "Dmrt-2" = "Dmrt-11E",
                                             "hyphally-cell-wall-protein-3" = "Ches-1","fox-I1c" = "Hcm-1", "slp" = "Slp", "fd" = "Fd", "Fd-4" = "Fd-4/Fd-96Ca", "Fd-5" = "Fd-5/Fd-96Cb","fkh" = "Fkh", "croc" = "Croc", "Fd-3" = "Fd-3/Fd-59A", "Fd-2" = "Fd-2/Fox-L1", "biniou" = "Biniou", "fox" = "Fox", "Fox-N2" = "Jumeau",
                                             "nan" = "Sox-21B", "NA_GATA-zinc-finger-10" = "Sox-15", "GATA-zinc-finger-10" = "Sox-N", "dicha" = "Dicha", "sox" = "Sox", "Sox-4-8-10" = "Sox-100B", "Sox-1a-14-21a_Sox-14" = "Sox-21A"))) %>%
  mutate(V2 = stringr::str_extract(V2, "[^_]+$"))

gene_conversion

# rename genes
diversity_data_SRGs <- diversity_data %>%
  filter(stringr::str_detect(gene, stringr::regex("dmrt|sox|fox", ignore_case = TRUE))) %>%
  dplyr::mutate(gene = stringr::str_extract(gene, "[^/]+$")) %>%
  dplyr::mutate(gene = stringr::str_extract(gene, "^[^_]+_[^_]+_[^_]+"))

diversity_data_SRGs

# define species to keep
species_toKeep <- c("Agam", "Dbus", "Dalb", "Dgri", "Dari", "Dhyd", "Dwil", "Dmir", "Dpse", "Dbip", "Dana", "Dkik", "Dser", "Dele", "Dsuz", "Dere", "Dmel", "Dsec")

# get counts of genes
gene_counts <- occurrence_data %>%
  filter(species %in% species_toKeep) %>%
  pivot_longer(cols = names(occurrence_data[,-1]), names_to = "gene") %>%
  aggregate(value ~ gene, FUN = sum)

# create a joined dataset
joined_data <- left_join(diversity_data_SRGs, gene_conversion, by = join_by(gene == V1)) %>%
  drop_na() %>%
  left_join(setNames(gene_counts, c("gene", "sum")), by = join_by(V2 == gene)) %>%
  rename(gene_names = V2, disco_genes = gene)

joined_data

joined_data[nrow(joined_data) + 1,] <- list(NA, NA, NA, "dmrt5", NA)
joined_data[nrow(joined_data) + 1,] <- list(NA, NA, NA, "foxZ", NA)

joined_data$quant <- as.factor(joined_data$quant)


#######################################
#     PLOT DIVERSITY DATA OF SRGs     #
#######################################

# set gene order to plot
# sox-1 and sox-5 are not present because didn't make through the alignment filter step
gene_order <- as.factor(c("Dsx", "Dmrt-99B", "Dmrt-93B", "Dmrt-11E", "dmrt5",
                          "Sox-15", "Sox-100B", "Sox-14", "Dichaete", "Sox-21B", "Sox-N", "Sox-21A",  
                          "foxZ", "Fox-3F", "Ches-1", "Jumeau", "Hcm-1", "Fox-P", "Fox-O", "Fd-102C", "Biniou", "Fd-2/Fox-L1", "Fd-19B", "Slp-1", "Slp-2", "Fd-3/Fd-59A", "Crocodile", "Fkh", "Fd-5/Fd-96Cb", "Fd-4/Fd-96Ca"))

gene_order_toplot <- factor(joined_data$gene_names,
                            levels = gene_order)

# plot median diversity (remove one Fox-B1/2, which is a duplicate)
plot_points <- joined_data %>%
  # filter(stringr::str_detect(gene_names, "Sox")) %>%
  
  ggplot(aes(y = median, x = gene_order_toplot, colour = quant, fill = quant)) +
  
  geom_segment(aes(y = -Inf, x = gene_order_toplot, yend = median, xend = gene_order_toplot),
               col = "#706d73", alpha = 0.3, linewidth = 0.8) +
  
  geom_point(aes(size = sum),
             shape = 21, stroke = 1.5) +
  
  geom_text(aes(y = median, x = gene_order_toplot), label = gene_order_toplot,
            col = "#706d73", size = 3, fontface = "italic",
            angle = 90, hjust = 0, nudge_y = 0.15) +
  
  coord_cartesian(clip = "off") +
  
  xlab("") +
  
  scale_y_continuous(limits = c(-0.11, 3), breaks = seq(0, 3, 0.5)) +
  scale_x_discrete(breaks = gene_order[-c(8,29)]) +
  
  scale_colour_manual(values = c("#0e65d1", "#b353df", "#cf1c47"),
                      labels = c("< 95%", paste0(quantile*100, "%")),
                      guide = "none") +
  
  scale_fill_manual(values = c("#7eaee8", "#d79ef1", "#ee809a"),
                    labels = c("< 95%", paste0(quantile*100, "%")),
                    guide = "none") +
  
  scale_size_continuous(range = c(5, 10),
                        breaks = c(11, 19),
                        name = "Quantiles Number of genes") +
  
  # guides(colour = guide_legend(title = "Quantiles")) +
  
  theme_minimal() +
  theme(axis.title.y = element_blank(),
	axis.line.y = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.y = element_blank(),
        axis.line = element_line(linewidth = 0.5))

plot_points


######################
#     PLOT PANEL     #
######################

empty_plot <- ggplot() + theme_void()

final_panel <- ggpubr::ggarrange(ggpubr::ggarrange(plot_density, plot_points,
                                                   nrow = 1, ncol = 2,
                                                   align = "hv",
                                                   # common.legend = TRUE,
                                                   legend = "top",
                                                   widths = c(0.3,1),
                                                   labels = c("A", "B")),
                                 
                                 ggpubr::ggarrange(empty_plot, panel_geneTrees, empty_plot,
                                                   nrow = 1, ncol = 3,
                                                   align = "hv",
                                                   widths = c(0.3,1,0.1)),
                                 ncol = 1,
                                 nrow = 2,
                                 align = "hv",
                                 heights = c(2, 0.3))

final_panel


########################
#     SAVE OUTPUTS     #
########################

write.table(diversity_data, distance_median_values_quant, quote = FALSE, sep = "\t", row.names = FALSE)

ggsave(diversity_panel,
       plot = final_panel, device = "pdf",
       dpi = 300, height = 9, width = 12, units = ("in"), bg = 'transparent')

