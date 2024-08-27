library(tidyr)
library(ggplot2)
library(ggtree)


######################
#     READ INPUT     #
######################

# occurrence matrix
occurrence_data <- "06_possvm_orthology/01_plot_occurrence/possvm_orthology_all_withOUT.tsv"

# trees
speciesTree_data <- read.tree(file = "00_input/species_tree_ALL.nwk")
dmrtTree_data <- read.tree(text = "(Dsx,(Dmrt-99B,(Dmrt-93B,Dmrt-11E)));")
soxTree_data <- read.tree(text = "(Sox-15,(Sox-100B,(Sox-14,(Dichaete,(Sox-21B,(Sox-N,Sox-21A))))));")
foxTree_data <- read.tree(text = "(Fox-3F,((Ches-1,Jumeau),((Hcm-1,(Fox-P,Fox-O)),(Fd-102C,(Biniou,(Fd-2/Fox-L1,((Fd-19B,(Slp-1,Slp-2)),(Fd-3/Fd-59A,(Crocodile,(Fkh,(Fd-5/Fd-96Cb,Fd-4/Fd-96Ca)))))))))));")


############################
#     OUTPUT FILENAMES     #
############################

occurrence_matrix_withOutgroups <- "06_possvm_orthology/01_plot_occurrence/occurrence_matrix_withOutgroups.tsv"
SRGs_occurrence_panel <- "06_possvm_orthology/01_plot_occurrence/SRGs_occurrence_panel.png"


#####################
#   FUNCTIONS   #
#####################

compute_occurrenceMatrix <- function(possvm_orthology_file, species_order) {
 
 # read in the possvm orthology file
 possvm_orthology <- read.table(possvm_orthology_file, header = FALSE, sep = "\t")
 
 # rename columns
 names(possvm_orthology) <- c("species", "orthogroup", "reference_gene")
 
 # get a count matrix
 gene_occurrence <- as.data.frame.matrix(table(possvm_orthology[c("species","reference_gene")]))
 
 # filter out orthogroups with less than 50% of species
 gene_occurrence_filter <- gene_occurrence[, colSums(gene_occurrence == 0)/nrow(gene_occurrence) <= 0.55, drop = FALSE] %>%
  tibble::rownames_to_column(var = "species") %>%
  dplyr::arrange(factor(species, levels = species_order))
 
 # # rename headers
 names(gene_occurrence_filter) <- sub("biniou", "Biniou", colnames(gene_occurrence_filter))
 names(gene_occurrence_filter) <- sub("crocodile", "Crocodile", colnames(gene_occurrence_filter))
 names(gene_occurrence_filter) <- sub("dichaete", "Dichaete", colnames(gene_occurrence_filter))
 names(gene_occurrence_filter) <- sub("dmrt-11E/dmrt-2", "Dmrt-11E", colnames(gene_occurrence_filter))
 names(gene_occurrence_filter) <- sub("dmrt-93B/dmrt-A2", "Dmrt-93B", colnames(gene_occurrence_filter))
 names(gene_occurrence_filter) <- sub("dmrt-99B/dmrt-A2/dmrt-A2-like", "Dmrt-99B", colnames(gene_occurrence_filter))
 names(gene_occurrence_filter) <- sub("dsx", "Dsx", colnames(gene_occurrence_filter))
 names(gene_occurrence_filter) <- sub("fd-102C/fox-C1-A/nan", "Fd-102C", colnames(gene_occurrence_filter))
 names(gene_occurrence_filter) <- sub("fd-19B/slp-1", "Fd-19B", colnames(gene_occurrence_filter))
 names(gene_occurrence_filter) <- sub("fd-2/fox-L1", "Fd-2/Fox-L1", colnames(gene_occurrence_filter))
 names(gene_occurrence_filter) <- sub("fd-3/fd-59A", "Fd-3/Fd-59A", colnames(gene_occurrence_filter))
 names(gene_occurrence_filter) <- sub("fd-4/fd-96Ca", "Fd-4/Fd-96Ca", colnames(gene_occurrence_filter))
 names(gene_occurrence_filter) <- sub("fd-5/fd-96Cb", "Fd-5/Fd-96Cb", colnames(gene_occurrence_filter))
 names(gene_occurrence_filter) <- sub("fkh", "Fkh", colnames(gene_occurrence_filter))
 names(gene_occurrence_filter) <- sub("fox-3F/homeobox-5/hormone-receptor-4", "Fox-3F", colnames(gene_occurrence_filter))
 names(gene_occurrence_filter) <- sub("fox-I1c/hcm-1/hnf-3b/nan", "Hcm-1", colnames(gene_occurrence_filter))
 names(gene_occurrence_filter) <- sub("fox-N2/jumeau/nan", "Jumeau", colnames(gene_occurrence_filter))
 names(gene_occurrence_filter) <- sub("fox-3F/homeobox-5/hormone-receptor-4", "Fox-3F", colnames(gene_occurrence_filter))
 names(gene_occurrence_filter) <- sub("fox-O", "Fox-O", colnames(gene_occurrence_filter))
 names(gene_occurrence_filter) <- sub("fox-P/fox-P1/fox-P4", "Fox-P", colnames(gene_occurrence_filter))
 names(gene_occurrence_filter) <- sub("GATA-zinc-finger-10/sox-1/sox-N", "Sox-N", colnames(gene_occurrence_filter))
 names(gene_occurrence_filter) <- sub("hyphally-cell-wall-protein-3/serine-rich-repeat-protein", "Ches-1", colnames(gene_occurrence_filter))
 names(gene_occurrence_filter) <- sub("nan/sox-1/sox-21b", "Sox-21B", colnames(gene_occurrence_filter))
 names(gene_occurrence_filter) <- sub("slp", "Slp", colnames(gene_occurrence_filter))
 names(gene_occurrence_filter) <- sub("sox-10/sox-4/sox-8/sox100B", "Sox-100B", colnames(gene_occurrence_filter))
 names(gene_occurrence_filter) <- sub("sox-14/sox-1a/sox-21a", "Sox-21A", colnames(gene_occurrence_filter))
 names(gene_occurrence_filter) <- sub("sox", "Sox", colnames(gene_occurrence_filter))
 
 return(gene_occurrence_filter)
 
}


plotOccurrence <- function(occurrence_matrix, species_order, speciesNames_toPlot, gene_order) {
 
 # create a tibble to plot
 gene_occurrence_filter_longer <- occurrence_matrix %>%
  # tibble::rownames_to_column(var = "species") %>%
  tidyr::pivot_longer(!species, names_to = "gene", values_to = "count")
 
 species_order_toplot <- factor(gene_occurrence_filter_longer$species,
                                levels = species_order)
 
 gene_order_toplot <- factor(gene_occurrence_filter_longer$gene,
                             levels = gene_order)
 
 plot <- ggplot(gene_occurrence_filter_longer,
                aes(x = gene_order_toplot, y = species_order_toplot,
                    fill = stage(count, after_scale = alpha(fill, 0.7)), colour = count)) +

  # annotate drosophila
  annotate(xmin = -Inf, xmax = Inf,
           ymin = -Inf, ymax = 1.5,
           geom = "rect", fill = "#324542", alpha = 0.1) +
  
  
  geom_point(aes(size = ifelse(count==0, NA, count)), shape = 21, stroke = 0.9) +
  
  guides(colour = guide_colourbar(barwidth = 2, direction = "horizontal"),
         fill = guide_colourbar(barwidth = 2, direction = "horizontal")) +

  coord_cartesian(clip = "off") +
  
  geom_text(gene_occurrence_filter_longer,
            mapping = aes(x = gene_order_toplot, y = species_order_toplot, label = ifelse(count > 1, count, NA)),
            size = 2,
            colour = "white",
            fontface = 2,
            show.legend = FALSE) +
  
  scale_y_discrete(position = "right", labels = speciesNames_toPlot) +
  scale_x_discrete(position = "top",
                   breaks = gene_order[-c(5,13)]) +
  
  scale_size_continuous(range = c(4,4),
             # breaks = c(1,seq(4,16,4)),
             guide = "none") +
  
  scale_fill_gradient(low = "#639585",
                      high = "#a59955",
                      name = stringr::str_wrap("Number of genes", width = 8),
                      breaks = c(1,2),
                      limits = c(1,2)) +
  scale_colour_gradient(low = "#639585",
                        high = "#a59955",
                        name = stringr::str_wrap("Number of genes", width = 8),
                        breaks = c(1,2),
                        limits = c(1,2)) +
  theme_minimal() +
  theme(panel.border = element_blank(),
     axis.text = element_text(color = "black"),
     axis.title.x = element_blank(),
     axis.title.y = element_blank(),
     axis.text.x = element_text(face = "italic", size = 8, angle = 45, hjust = 0),
     axis.text.y = element_text(face = "italic", size = 8),
     legend.text = element_text(size = 7),
     legend.title = element_text(size = 7),
     legend.title.align = 0.5
     )
 
 return(plot)
}


###################################
#     GENERATE MATRIX TO PLOT     #
###################################

# plot SRG occurrence with also outgroups
species_order_withOUT <- as.factor(c("Agam", "Dbus", "Dalb", "Dgri", "Dari", "Dhyd", "Dwil", "Dmir", "Dpse", "Dbip", "Dana", "Dkik", "Dser", "Dele", "Dsuz", "Dere", "Dmel", "Dsec"))

species_fullNames <- as.factor(c("A. gambiae", "D. busckii", "D. albomicans", "D. grimshawi", "D. arizonae", "D. hydei", "D. willistoni", "D. miranda", "D. pseudoobscura", "D. bipectinata", "D. ananassae", "D. kikkawai", "D. serrata", "D. elegans", "D. suzukii", "D. erecta", "D. melanogaster", "D. sechellia"))

gene_order <- as.factor(c("Dsx", "Dmrt-99B", "Dmrt-93B", "Dmrt-11E", "dmrt5",
                          "Sox-15", "Sox-100B", "Sox-14", "Dichaete", "Sox-21B", "Sox-N", "Sox-21A",  
                          "foxZ", "Fox-3F", "Ches-1", "Jumeau", "Hcm-1", "Fox-P", "Fox-O", "Fd-102C", "Biniou", "Fd-2/Fox-L1", "Fd-19B", "Slp-1", "Slp-2", "Fd-3/Fd-59A", "Crocodile", "Fkh", "Fd-5/Fd-96Cb", "Fd-4/Fd-96Ca"))

occurrence_matrix_withOUT <- compute_occurrenceMatrix(occurrence_data, species_order_withOUT)

occurrence_matrix_withOUT <- occurrence_matrix_withOUT[occurrence_matrix_withOUT$species != "Neur",]

names(occurrence_matrix_withOUT)

# add plot spacer genes
occurrence_matrix_withOUT$dmrt5 <- NA
occurrence_matrix_withOUT$foxZ <- NA


############################
#     PLOT OCCURRENCES     #
############################

plot_withOUT <- plotOccurrence(occurrence_matrix_withOUT, species_order_withOUT, species_fullNames, gene_order)

plot_withOUT

# save legend to a separate object
legend <- ggpubr::get_legend(plot_withOUT)

# remove legend from plot
plot_withOUT <- plot_withOUT + guides(colour = "none", fill = "none")


######################
#     PLOT TREES     #
######################

# load species tree
speciesTree_plot <- ggtree(speciesTree_data) +
 # geom_tiplab(angle = 45, offset = 1, size = 2) +
 ggpubr::theme_transparent() 
speciesTree_plot


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


######################
#     PLOT PANEL     #
######################

# plot panel
empty_plot <- ggplot() + theme_void()

panel_geneTrees <- ggpubr::ggarrange(dmrtTree_plot, soxTree_plot, foxTree_plot, empty_plot,
                   ncol = 4, nrow = 1,
                   align = "hv",
                   widths = c(0.4, 0.7, 1.5, 0.7))

panel_geneTrees


panel_speciesTree <- ggpubr::ggarrange(empty_plot, speciesTree_plot,
                       ncol = 1, nrow = 2,
                       align = "hv",
                       heights = c(0.5,2.4))
panel_speciesTree

panel_full <- ggpubr::ggarrange(panel_speciesTree, plot_withOUT, legend, panel_geneTrees,
                ncol = 2, nrow = 2,
                align = "hv",
                heights = c(2,0.3),
                widths = c(0.2,2))

panel_full


########################
#     SAVE OUTPUTS     #
########################

write.table(occurrence_matrix_withOUT[,-c(30,31)], file = occurrence_matrix_withOutgroups,
            quote = FALSE, sep = "\t", row.names = FALSE)

ggsave(SRGs_occurrence_panel,
       plot = panel_full, device = "png",
       dpi = 300, height = 5, width = 9, units = ("in"), bg = 'transparent')
