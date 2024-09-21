#!/usr/bin/env Rscript

library(dplyr)
library(ggplot2)
library(ggtree)


######################
#     READ INPUT     #
######################

# occurrence matrix
occurrence_data <- "06_possvm_orthology/01_plot_occurrence/possvm_orthology_all_withOUT.tsv"

# trees
speciesTree_data <- read.tree(file = "00_input/species_tree_ALL.nwk")
dmrtTree_data <- read.tree(text = "(Dmrt-B1,((Dmrt-2,Dmrt-C2),(Dmrt-3,(Dmrt-1,(Dmrt-A1,Dmrt-A2)))));")
soxTree_data <- read.tree(text = "(Sox-18,(Sox-17,(Sox-7,((Sox-8,(Sox-9,Sox-10)),((Sox-30,(Sox-5,(Sox-13,Sox-6))),((Sox-4,(Sox-12,Sox-11)),(Sox-15,((Sox-3,Sry),(Sox-1,(Sox-2,(Sox-14,Sox-21)))))))))));")
foxTree_data <- read.tree(text = "(Fox-M1,(((Fox-O6,(Fox-O3/O3B,(Fox-O1,Fox-O4))),(((Fox-J3,Fox-J2),(Fox-J1,(Fox-K2,Fox-K1))),((Fox-P3,(Fox-P1,(Fox-P4,Fox-P2))),((Fox-N3,Fox-N2),(Fox-N1,(Fox-OG59/NA,Fox-R1)))))),((Fox-H1,(Fox-Q1,(Fox-F2,Fox-F1))),((((Fox-OG21/NA,Fox-L2),(Fox-I3,(Fox-I2,Fox-I1))),(Fox-L3,(Fox-L1,(Fox-S1,(Fox-C1,Fox-C2))))),(Fox-G1,(((Fox-B2,Fox-B1),(Hnf-3b/Fox-A2,(Hnf-3a/Fox-A1,Hnf-3g/Fox-A3))),((Fox-E3,Fox-E1),(Fox-D4,(Fox-D3,(Fox-D1,Fox-D2))))))))));") %>%
	ape::drop.tip(c("Fox-OG59/NA", "Fox-OG21/NA", "Fox-K2"))


############################
#     OUTPUT FILENAMES     #
############################

occurrence_matrix_withOutgroups <- "06_possvm_orthology/01_plot_occurrence/occurrence_matrix_withOutgroups.tsv"
SRGs_occurrence_panel <- "06_possvm_orthology/01_plot_occurrence/SRGs_occurrence_panel.png"

#####################
#     FUNCTIONS     #
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
 
 # rename headers
 names(gene_occurrence_filter) <- sub("dmrt", "Dmrt-", colnames(gene_occurrence_filter))
 names(gene_occurrence_filter) <- sub("fox", "Fox-", colnames(gene_occurrence_filter))
 names(gene_occurrence_filter) <- sub("sox", "Sox", colnames(gene_occurrence_filter))
 names(gene_occurrence_filter) <- sub("Sox-3/sry", "Sry", colnames(gene_occurrence_filter))
 names(gene_occurrence_filter) <- sub("/fox", "/", colnames(gene_occurrence_filter))
 names(gene_occurrence_filter) <- sub("hnf3a", "Hnf-3a/Fox-A1", colnames(gene_occurrence_filter))
 names(gene_occurrence_filter) <- sub("hnf3b", "Hnf-3b/Fox-A2", colnames(gene_occurrence_filter))
 names(gene_occurrence_filter) <- sub("hnf3g", "Hnf-3g/Fox-A3", colnames(gene_occurrence_filter))
 
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

  # annotate mammals
  annotate(xmin = -Inf, xmax = Inf,
       ymin = -Inf, ymax = 1.5,
       geom = "rect", fill = "#324542", alpha = 0.1) +
  
  
  geom_point(aes(size = ifelse(count==0, NA, count)), shape = 21, stroke = 0.9) +
  
  guides(colour = guide_colourbar(barwidth = 5, direction = "horizontal"),
         fill = guide_colourbar(barwidth = 5, direction = "horizontal")) +

  coord_cartesian(clip = "off") +
  
  geom_text(gene_occurrence_filter_longer,
            mapping = aes(x = gene_order_toplot, y = species_order_toplot, label = ifelse(count > 1, count, NA)),
            size = 2,
            colour = "white",
            fontface = 2,
            show.legend = FALSE) +
  
  scale_y_discrete(position = "right", labels = speciesNames_toPlot) +
  scale_x_discrete(position = "top",
                   breaks = gene_order[-c(8,29)]) +
  
  scale_size_continuous(range = c(4,4),
             # breaks = c(1,seq(4,16,4)),
             guide = "none") +
  
  scale_fill_gradient2(low = "#0a9396",
             mid = "#ee9b00",
             midpoint = 3.3,
             high = "#bb3e03",
             name = stringr::str_wrap("Number of genes", width = 8),
             breaks = c(1,4,6),
             limits = c(1,6)) +
  scale_colour_gradient2(low = "#0a9396",
              mid = "#ee9b00",
              midpoint = 3.3,
              high = "#bb3e03",
              name = stringr::str_wrap("Number of genes", width = 8),
              breaks = c(1,4,6),
              limits = c(1,6)) +
  theme_minimal() +
  theme(panel.border = element_blank(),
     axis.text = element_text(color = "black"),
     axis.title.x = element_blank(),
     axis.title.y = element_blank(),
     axis.text.x = element_text(face = "italic", size = 8, angle = 45, hjust = 0),
     axis.text.y = element_text(face = "italic", size = 8),
     legend.text = element_text(size = 7),
     legend.title = element_text(size = 10),
     legend.title.align = 0.5
     )
 
 return(plot)
}


###################################
#     GENERATE MATRIX TO PLOT     #
###################################


# plot SRG occurrence with also outgroups
species_order_withOUT <- as.factor(c("Ggal", "Oana", "Mdom", "Shar", "Cdid", "Dnov", "Oafe", "Casi", "Emax", "Tman", "Lcat", "Hsap", "Cimi", "Opri", "Scar", "Cpor", "Mmus", "Drot", "Pgig", "Rfer", "Cdro", "Pafr", "Bbub", "Hamp", "Bmus", "Ttru", "Csim", "Equa", "Mjav", "Ptig", "Clup", "Amel", "Mang"))

species_fullNames <- as.factor(c("G. gallus", "O. anatinus", "M. domestica", "S. harrisii", "C. didactylus", "D. novemcinctus", "O. afer", "C. asiatica", "E. maximus", "T. manatus", "L. catta", "H. sapiens", "C. imitator", "O. princeps", "S. carolinesis", "C. porcellus", "M. musculus", "D. rotundus", "P. giganteus", "R. ferrumequinum", "C. dromedarius", "P. africanus", "B. bubalis", "H. amphibius", "B. musculus", "T. truncatus", "C. simum", "E. quagga", "M. javanica", "P. tigris", "C. lupus", "A. melanoleuca", "M. angustirostris"))

gene_order <- as.factor(c(rev(c("dmrt5", "Dmrt-A2", "Dmrt-A1", "Dmrt-1", "Dmrt-3", "Dmrt-C2", "Dmrt-2", "Dmrt-B1")),
                          rev(c("Sox-21", "Sox-14", "Sox-2", "Sox-1", "Sry", "Sox-3", "Sox-15", "Sox-11", "Sox-12", "Sox-4", "Sox-5", "Sox-6", "Sox-13", "Sox-30", "Sox-8", "Sox-10", "Sox-9", "Sox-7", "Sox-17", "Sox-18")),  
                          rev(c("Fox-D2", "Fox-D1", "Fox-D3", "Fox-D4", "Fox-E1", "Fox-E3", "Hnf-3g/Fox-A3", "Hnf-3a/Fox-A1", "Hnf-3b/Fox-A2", "Fox-B1", "Fox-B2", "Fox-G1", "Fox-C2", "Fox-C1", "Fox-S1", "Fox-L1", "Fox-L3", "Fox-I1", "Fox-I2", "Fox-I3", "Fox-OG21/NA", "Fox-L2", "Fox-F1", "Fox-F2", "Fox-Q1", "Fox-H1", "Fox-OG59/NA", "Fox-R1", "Fox-N1", "Fox-N2", "Fox-N3", "Fox-P2", "Fox-P4", "Fox-P1", "Fox-P3", "Fox-K1", "Fox-K2", "Fox-J1", "Fox-J2", "Fox-J3", "Fox-O4", "Fox-O1", "Fox-O3/O3B", "Fox-O6", "Fox-M1", "foxZ"))))

occurrence_matrix_withOUT <- compute_occurrenceMatrix(occurrence_data, species_order_withOUT)

# switch sox3 and sry from platypus and chicken
occurrence_matrix_withOUT[c(1,2),63] <- 0
occurrence_matrix_withOUT[c(1,2),62] <- 1

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
                   widths = c(0.4, 0.9, 2, 0.3))

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

write.table(occurrence_matrix_withOUT[,-c(71,72)], file = occurrence_matrix_withOutgroups,
            quote = FALSE, sep = "\t", row.names = FALSE)

ggsave(SRGs_occurrence_panel,
       plot = panel_full, device = "png",
       dpi = 300, height = 7, width = 15, units = ("in"), bg = 'transparent')
