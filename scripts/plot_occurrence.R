library(tidyr)
library(ggplot2)
library(ggtree)

#####################
#     FUNCTIONS     #
#####################

compute_occurrenceMatrix <- function(possvm_orthology_file, species_order) {
  
  # read in the possvm orthology file
  possvm_orthology <- read.table(possvm_orthology_file, header = FALSE, sep = "\t")
  
  # rename columns
  names(possvm_orthology) <- c("species", "orthogroup", "reference_gene")
  
  # get a count matrix
  gene_occurrence <- as.data.frame.matrix(table(possvm_orthology[c("species","orthogroup")]))
  
  # filter out orthogroups with less than 50% of species
  gene_occurrence_filter <- gene_occurrence[, colSums(gene_occurrence == 0)/nrow(gene_occurrence) <= 0.6, drop = FALSE] %>%
    tibble::rownames_to_column(var = "species") %>%
    dplyr::arrange(factor(species, levels = species_order))
  
  # rename headers
  names(gene_occurrence_filter) <- sub(":NA", "ogNA", colnames(gene_occurrence_filter))
  names(gene_occurrence_filter) <- sub("OG.+:", "", colnames(gene_occurrence_filter))
  names(gene_occurrence_filter) <- sub("_fkh.+/", "_", colnames(gene_occurrence_filter))
  names(gene_occurrence_filter) <- sub("_dmrt-", "", colnames(gene_occurrence_filter))
  names(gene_occurrence_filter) <- sub("_sox-", "", colnames(gene_occurrence_filter))
  names(gene_occurrence_filter) <- sub("_fox-", "", colnames(gene_occurrence_filter))
  names(gene_occurrence_filter) <- sub("dmrt_dmd-8/mab-23", "dmrt1L", colnames(gene_occurrence_filter))
  names(gene_occurrence_filter) <- make.unique(names(gene_occurrence_filter))
  names(gene_occurrence_filter) <- sub("\\.", "", colnames(gene_occurrence_filter))
  names(gene_occurrence_filter) <- sub("fox_", "fox", colnames(gene_occurrence_filter))
  names(gene_occurrence_filter) <- sub("sox_", "sox", colnames(gene_occurrence_filter))
  names(gene_occurrence_filter) <- sub("4-5", "4/5", colnames(gene_occurrence_filter))
  names(gene_occurrence_filter) <- sub("dmrt", "Dmrt-", colnames(gene_occurrence_filter))
  names(gene_occurrence_filter) <- sub("fox", "Fox-", colnames(gene_occurrence_filter))
  names(gene_occurrence_filter) <- sub("sox", "Sox-", colnames(gene_occurrence_filter))
  names(gene_occurrence_filter) <- sub("N3", "N2/3", colnames(gene_occurrence_filter))
  names(gene_occurrence_filter) <- sub("N4", "N1/4", colnames(gene_occurrence_filter))
  names(gene_occurrence_filter) <- sub("2-3", "2/3", colnames(gene_occurrence_filter))
  names(gene_occurrence_filter) <- sub("ogNA", "/NA", colnames(gene_occurrence_filter))
  
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
    # annotate DSFGs
    # annotate(xmin = rep(-Inf, 2), xmax = rep(Inf, 2), 
    #          ymin = c(-Inf, 30.5), ymax = c(7.5, Inf), 
    #          geom = "rect", alpha = 0.2) +
     
    # annotate bivalves
    annotate(xmin = -Inf, xmax = Inf,
             ymin = 10.5, ymax = Inf,
             geom = "rect", alpha = 0.1) +
    
    
    geom_point(aes(size = ifelse(count==0, NA, count)), shape = 21, stroke = 0.9) +
    
    guides(colour = guide_colourbar(barheight = 5, title.vjust = 5),
           fill = guide_colourbar(barheight = 5, title.vjust = 5)) +

    coord_cartesian(clip = "off") +
    
    geom_text(gene_occurrence_filter_longer,
              mapping = aes(x = gene_order_toplot, y = species_order_toplot, label = ifelse(count > 1, count, NA)),
              size = 2,
              colour = "white",
              fontface = 2,
              show.legend = FALSE) +
    
    scale_y_discrete(position = "right", labels = speciesNames_toPlot) +
    scale_x_discrete(position = "top",
                     breaks = gene_order[-c(5,12)]) +
    
    scale_size_continuous(range = c(4,4),
                          # breaks = c(1,seq(4,16,4)),
                          guide = "none") +
    
    scale_fill_gradient2(low = "#0a9396",
                         mid = "#ee9b00",
                         midpoint = 8,
                         high = "#bb3e03",
                         name = stringr::str_wrap("Number of genes", width = 8),
                         breaks = c(1,6,11,16),
                         limits = c(1,16)) +
    scale_colour_gradient2(low = "#0a9396",
                           mid = "#ee9b00",
                           midpoint = 8,
                           high = "#bb3e03",
                           name = stringr::str_wrap("Number of genes", width = 8),
                           breaks = c(1,6,11,16),
                           limits = c(1,16)) +
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


#######################
#     ACTUAL CODE     #
#######################

# plot SRG occurrence with also outgroups
species_order_withOUT <- as.factor(c("Hsap", "Dmel", "Cele",
                                     "Osin","Obim",
                                     "Hruf","Gaeg","Pcan","Bgla","Acal",
                                     "Mmar","Pstr","Mner","Cpli","Hbia","Lorb","Tsqu","Pgen","Scon","Sgra","Pcor","Dpol","Mare","Mchi","Cflu","Poku","Amar","Csin","Mmer","Rphi","Rdec","Tgra","Sbro","Pyes","Pmax","Apur","Airc","Ppur","Mphi","Mmod","Pvir","Mcal","Mcor","Medu","Mgal","Pmar","Apec","Sglo","Oedu","Cvir","Cari","Cgig","Cang"))

species_fullNames <- as.factor(c("H. sapiens", "D. melanogaster", "C. elegans",
                                 "O. sinensis", "O. bimaculoides",
                                 "H. rufescens", "G. aegis", "P. canaliculata", "B. glabrata", "A. californica",
                                 "M. margaritifera", "P. streckersoni", "M. nervosa", "C. plicata*", "H. bialata", "L. orbiculatus*", "T. squamosa*", "P. generosa*", "S. constricta", "S. grandis*", "P. coreanum*", "D. polymorpha", "M. arenaria", "M. chinensis*", "C. fluminea*", "P. okutanii*", "A. marissinica", "C. sinensis", "M. mercenaria", "R. philippinarum", "R. decussatus*", "T. granosa", "S. broughtonii", "P. yessoensis", "P. maximus", "A. purpuratus", "A. irradians", "P. purpuratus*", "M. philippinarum", "M. modiolus*", "P. viridis", "M. californianus", "M. coroscus", "M. edulis", "M. galloprovincialis", "P. margaritifera*", "A. pectinata*", "S. glomerata", "O. edulis", "C. virginica", "C. ariakensis", "C. gigas", "C. angulata"))

gene_order <- as.factor(c("Dmrt-1L", "Dmrt-3", "Dmrt-2", "Dmrt-4/5", "dmrt5",
                          "Sox-H", "Sox-D", "Sox-B1/2", "Sox-C", "Sox-F", "Sox-E",   
                          "foxZ", "Fox-OG2/NA", "Fox-O", "Fox-P", "Fox-J2/3", "Fox-OG13/NA", "Fox-N2/3", "Fox-OG16/NA", "Fox-N1/4", "Fox-J1", "Fox-OG15/NA", "Fox-Q2", "Fox-OG28/NA", "Fox-G", "Fox-L2", "Fox-L1", "Fox-C", "Fox-F", "Fox-H", "Fox-E", "Fox-D", "Fox-OG39/NA", "Fox-B", "Fox-A"))

occurrence_matrix_withOUT <- compute_occurrenceMatrix("06_possvm_orthology/01_plot_occurrence/possvm_orthology_all_withOUT.tsv",
                                                      species_order_withOUT)
# remove dmrt-1l from Cele
occurrence_matrix_withOUT["3","Dmrt-1L"] <- 0

# add dmrt-1l to mytilida
occurrence_matrix_withOUT["38","Dmrt-1L"] <- 1 # Ppur
occurrence_matrix_withOUT["41","Dmrt-1L"] <- 1 # Pvir
occurrence_matrix_withOUT["42","Dmrt-1L"] <- 1 # Mcal
occurrence_matrix_withOUT["43","Dmrt-1L"] <- 1 # Mcor
occurrence_matrix_withOUT["44","Dmrt-1L"] <- 1 # Medu
occurrence_matrix_withOUT["45","Dmrt-1L"] <- 1 # Mgal

# add plot spacer genes
occurrence_matrix_withOUT$dmrt5 <- NA
occurrence_matrix_withOUT$foxZ <- NA

# merge sox-B11 and sox-B1 in same column and rename
occurrence_matrix_withOUT <- occurrence_matrix_withOUT %>%
  dplyr::mutate(sum = rowSums(dplyr::across(c("Sox-B11", "Sox-B1"))), .keep = "unused")
names(occurrence_matrix_withOUT)[36] <- "Sox-B1/2"

plot_withOUT <- plotOccurrence(occurrence_matrix_withOUT, species_order_withOUT, species_fullNames, gene_order)

plot_withOUT

# save legend to a separate object
legend <- ggpubr::get_legend(plot_withOUT)

# remove legend from plot
plot_withOUT <- plot_withOUT + guides(colour = "none", fill = "none")


ggsave("06_possvm_orthology/01_plot_occurrence/SRGs_occurrence.png",
       plot = plot_withOUT, device = "png",
       dpi = 300, height = 10, width = 8, units = ("in"), bg = 'transparent')

ggsave("06_possvm_orthology/01_plot_occurrence/SRGs_occurrence.pdf",
       plot = plot_withOUT, device = "pdf",
       dpi = 300, height = 10, width = 8, units = ("in"), bg = 'transparent')


# load species tree
speciesTree_data <- read.tree("00_input/species_tree_ALL.nwk")
speciesTree_plot <- ggtree(speciesTree_data) +
  # geom_tiplab(angle = 45, offset = 1, size = 2) +
  ggpubr::theme_transparent() 
speciesTree_plot


# load gene trees
dmrtTree_data <- read.tree("06_possvm_orthology/01_plot_occurrence/dmrt_tree_essential.txt")
dmrtTree_plot <- ggtree(dmrtTree_data) +
  # geom_tiplab() +
  coord_flip() +
  ggpubr::theme_transparent() +
  theme(plot.margin = unit(c(0, 0, 0, 0), "points"))

dmrtTree_plot

soxTree_data <- read.tree("06_possvm_orthology/01_plot_occurrence/sox_tree_essential.txt") %>%
  ape::drop.tip(c("soxOG1", "soxOG0"))
soxTree_plot <- ggtree(soxTree_data) +
  # geom_tiplab() +
  coord_flip() +
  ggpubr::theme_transparent() +
  theme(plot.margin = unit(c(0, 0, 0, 0), "points"))

soxTree_plot

foxTree_data <- read.tree("06_possvm_orthology/01_plot_occurrence/fox_tree_essential.txt") %>%
  ape::drop.tip("foxK")
foxTree_plot <- ggtree(foxTree_data) +
  # geom_tiplab() +
  coord_flip() +
  ggpubr::theme_transparent() +
  theme(plot.margin = unit(c(0, 0, 0, 0), "points"))

foxTree_plot

# plot panel
empty_plot <- ggplot() + theme_void()

panel_geneTrees <- ggpubr::ggarrange(dmrtTree_plot, soxTree_plot, foxTree_plot, empty_plot,
                                     ncol = 4, nrow = 1,
                                     align = "hv",
                                     widths = c(0.4, 0.6, 2, 0.7))

panel_geneTrees


panel_speciesTree <- ggpubr::ggarrange(empty_plot, speciesTree_plot,
                                              ncol = 1, nrow = 2,
                                              align = "hv",
                                              heights = c(0.15,2.4))
panel_speciesTree

panel_full <- ggpubr::ggarrange(panel_speciesTree, plot_withOUT, legend, panel_geneTrees,
                                ncol = 2, nrow = 2,
                                align = "hv",
                                heights = c(2,0.3),
                                widths = c(0.4,2))

panel_full

ggsave("06_possvm_orthology/01_plot_occurrence/SRGs_occurrence_panel.png",
       plot = panel_full, device = "png",
       dpi = 300, height = 12, width = 9, units = ("in"), bg = 'transparent')

ggsave("06_possvm_orthology/01_plot_occurrence/SRGs_occurrence_panel.pdf",
       plot = panel_full, device = "pdf",
       dpi = 300, height = 12, width = 9, units = ("in"), bg = 'transparent')
