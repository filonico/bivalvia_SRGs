#!/usr/bin/env Rscript

library(tidyr)
library(dplyr)
library(ggplot2)
library(ggtree)


######################
#     READ INPUT     #
######################

# gene trees
dmrtTree_data <- read.tree(text = "(Dmrt-1L,(Dmrt-3,(Dmrt-2,Dmrt-4/5)));")
soxTree_data <- read.tree(text = "((Sox-OG0/NA,Sox-OG1/NA),(Sox-H,(Sox-D,(Sox-B1/2,(Sox-C,(Sox-F,Sox-E))))));") %>%
  ape::drop.tip(c("soxOG1", "soxOG0"))
foxTree_data <- read.tree(text = "(Fox-OG2/NA,((Fox-O,Fox-P),(Fox-J2/3,(((Fox-OG13/NA,Fox-N2/3),(Fox-OG16/NA,Fox-N1/4)),((Fox-K,Fox-J1),((Fox-OG15/NA,(Fox-Q2,Fox-OG28/NA)),(Fox-G,(Fox-L2,((Fox-L1,Fox-C),((Fox-F,Fox-H),((Fox-E,Fox-D),(Fox-OG39/NA,(Fox-B,Fox-A)))))))))))));") %>%
  ape::drop.tip("foxK")

# amino acid diversity file
diversity_file <- "13_distribution_divergence/distance_median_values.tsv"

# gene names conversion tables
file_paths <- list.files(path = "./10_SRG_decomposition/", pattern = "*conversion.tsv", full.names = TRUE)

# occurrence table
occurrence_data <- read.table("06_possvm_orthology/01_plot_occurrence/occurrence_matrix_withOutgroups.tsv",
                              header = TRUE, sep = "\t", check.names = FALSE)


############################
#     OUTPUT FILENAMES     #
############################

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
                                     widths = c(0.45, 0.6, 2))

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
  geom_line(aes(color = as.factor(quant)), linewidth = 1.3, lineend = "round", show.legend = FALSE) +
  
  coord_flip() +
  
  scale_colour_manual(values = c("#0e65d1", "#b353df", "#cf1c47"),
                      labels = c("< 95%", paste0(quantile*100, "%"))) +
  
  scale_fill_manual(values = c("#7eaee8", "#d79ef1", "#ee809a"),
                    labels = c("< 95%", paste0(quantile*100, "%"))) +  
  
  guides(fill = guide_legend(title = "Quantiles")) +
  
  scale_x_continuous(limits = c(0, 5), breaks = seq(0, 5, 1)) +
  
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

data_list <- lapply(file_paths, function(file) {
  read.table(file, header = FALSE, sep = "\t")
})

gene_conversion <- do.call(rbind, data_list)

# rename conversion tables
gene_conversion <- gene_conversion %>%
  mutate(V1 = stringr::str_extract(V1, "[^.]+"),
         V2 = stringr::str_extract(sub("_NA", "", sub("_NA", "ogNA", V2)), "[^/]+$")) %>%
  mutate(V2 = stringr::str_extract(V2, "[^_]+$"),
         V2 = stringr::str_replace_all(V2, c("dmrt" = "Dmrt", "sox" = "Sox", "fox" = "Fox",
                                             "mab-23" = "Dmrt-1L", "2-3" = "2/3", "1-4" = "1/4", "Dmrt-5" = "Dmrt-4/5", "Sox-B1" = "Sox-B1/2",
                                             "Fox-J3" = "Fox-J2/3", "Fox-M" = "OG2/NA",
                                             "og" = "/")),
         V2 = ifelse(stringr::str_detect(V1, stringr::regex("sox", ignore_case = TRUE)),
                     stringr::str_replace(V2, "OG", "Sox-OG"),
                     ifelse(stringr::str_detect(V1, stringr::regex("Fox", ignore_case = TRUE)),
                            stringr::str_replace(V2, "OG", "Fox-OG"),
                            V2)))
gene_conversion

# rename genes
diversity_data <- diversity_data %>%
  filter(stringr::str_detect(gene, stringr::regex("dmrt|sox|fox", ignore_case = TRUE))) %>%
  dplyr::mutate(gene = stringr::str_extract(gene, "[^/]+$")) %>%
  dplyr::mutate(gene = stringr::str_extract(gene, "^[^_]+_[^_]+_[^_]+"))

diversity_data

# define species to keep
species_toKeep <- c("Sbro","Airc","Apec","Amar","Cflu","Cvir","Cpli","Csin","Dpol","Hbia","Lorb","Mchi","Cgig","Mmar","Mner","Mmer","Pyes","Mmod","Mare","Mgal","Oedu","Pgen","Pmax","Pvir","Ppur","Poku","Pmar","Pcor","Pstr","Rdec","Rphi","Sglo","Scon","Sgra","Tgra","Tsqu")

# get counts of genes
gene_counts <- occurrence_data %>%
  filter(species %in% species_toKeep) %>%
  pivot_longer(cols = names(occurrence_data[,-1]), names_to = "gene") %>%
  aggregate(value ~ gene, FUN = sum)

# create a joined dataset
joined_data <- left_join(diversity_data, gene_conversion, by = join_by(gene == V1)) %>%
  drop_na() %>%
  left_join(setNames(gene_counts, c("gene", "sum")), by = join_by(V2 == gene)) %>%
  rename(gene_names = V2, disco_genes = gene)

joined_data

joined_data[nrow(joined_data) + 1,] <- list(NA, NA, NA, "dmrt5", NA)
joined_data[nrow(joined_data) + 1,] <- list(NA, NA, NA, "foxZ", NA)

# remove one Sox-B1/2 gene
joined_data <- joined_data[-6,]

joined_data$quant <- as.factor(joined_data$quant)


#######################################
#     PLOT DIVERSITY DATA OF SRGs     #
#######################################

# set gene order to plot
gene_order <- as.factor(c("Dmrt-1L", "Dmrt-3", "Dmrt-2", "Dmrt-4/5", "dmrt5",
                          "Sox-H", "Sox-D", "Sox-B1/2", "Sox-C", "Sox-F", "Sox-E",
                          "foxZ", "Fox-OG2/NA", "Fox-O", "Fox-P", "Fox-J2/3", "Fox-OG13/NA", "Fox-N2/3", "Fox-OG16/NA", "Fox-N1/4", "Fox-J1", "Fox-OG15/NA", "Fox-Q2", "Fox-OG28/NA", "Fox-G", "Fox-L2", "Fox-L1", "Fox-C", "Fox-F", "Fox-E", "Fox-D", "Fox-OG39/NA", "Fox-B", "Fox-A"))

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
            angle = 90, hjust = 0, nudge_y = 0.25) +
  
  coord_cartesian(clip = "off") +
  
  xlab("") +
  
  scale_y_continuous(limits = c(0, 5), breaks = seq(0, 5, 1)) +
  scale_x_discrete(breaks = gene_order[-c(5,12)]) +
  
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
  
  scale_size_continuous(range = c(5, 17),
                        breaks = c(17, 74),
                        name = "Quantiles Number of genes") +
  
  # guides(colour = guide_legend(title = "Quantiles")) +
  
  theme_minimal() +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_text(size = 8),
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

ggsave(diversity_panel,
       plot = final_panel, device = "pdf",
       dpi = 300, height = 9, width = 12, units = ("in"), bg = 'transparent')

