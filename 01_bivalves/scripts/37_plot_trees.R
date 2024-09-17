#!/usr/bin/env Rscript

library(ggtree)
library(ggplot2)

####################
#     FUNCTION     #
####################

# function to load, plot, and save tree
plot.tree <- function(input_filename, output_filename, output_height, output_width) {
  
  tree <- read.tree(file = input_filename)

  tree$tip.label <- ifelse(grepl("Hsap|Cele|Dmel", tree$tip.label),
                           paste0("*", tree$tip.label),
                           tree$tip.label)
  
  tree.plot <- ggtree(tree, size = 0.2) +
    
    geom_tiplab(size = 1, offset = 0.1) +
    geom_point(aes(color = as.numeric(label)), size = 0.8) +
    
    scale_color_gradient2(name = "Bootstrap values",
                          low = "#fe4a49", mid = "#fed766", high = "#009fb7", midpoint = 50,
                          limits = c(0,100), breaks = seq(0,100,25), na.value = "transparent") +
    
    theme(legend.position = "top",
          legend.title.position = "top",
          legend.title = element_text(size = 7, face = "bold", hjust = 0.5),
          legend.text = element_text(size = 7),
          legend.key.size = unit(0.5, "cm"))
  
  ggsave(output_filename,
         plot = tree.plot, device = "pdf",
         dpi = 300, height = output_height, width = output_width, units = ("in"), bg = 'white',
         limitsize = FALSE)
  
  
  return(tree.plot)
  
}


#####################
#     PLOT TREES    #
#####################

dmrt.tree <- plot.tree("06_possvm_orthology/dmrt_ALL_reduced_aligned_trim04.faa.treefile.ortholog_groups.newick",
                       "06_possvm_orthology/02_plot_trees/supp_fig_S1.pdf",
                       12, 9)

sox.tree <- plot.tree("06_possvm_orthology/sox_ALL_reduced_aligned_trim04.faa_rooted.treefile.ortholog_groups.newick",
                      "06_possvm_orthology/02_plot_trees/supp_fig_S2.pdf",
                      24, 12)

fox.tree <- plot.tree("06_possvm_orthology/fox_ALL_reduced_aligned_trim04.faa_rooted.treefile.ortholog_groups.newick",
                      "06_possvm_orthology/02_plot_trees/supp_fig_S3.pdf",
                      60, 15)

dmrt.onlyBivalves.tree <- plot.tree("05b_family_phylogeny_dmrt_onlyBivalves/dmrt_ALL_allBivalves_reduced_aligned_trim04.faa.rotted.treefile",
                                    "05b_family_phylogeny_dmrt_onlyBivalves/supp_fig_S11.pdf",
                                    9, 9)

soxB1_soxB2.tree <- plot.tree("05c_family_phylogeny_soxB12_onlyBivalves/soxB12_aligned_trim04.faa.rooted.treefile",
			      "05c_family_phylogeny_soxB12_onlyBivalves/supp_fig_S12.pdf",
			      7, 7)

foxM.tree <- plot.tree("05d_family_phylogey_foxMOP/foxMOP_aligned_trim04.faa.rooted.treefile",
                       "05d_family_phylogey_foxMOP/supp_fig_S13.pdf",
                       10, 7)

foxY.tree <- plot.tree("05e_family_phylogeny_foxYZ/foxYZ_aligned_trim04.faa.rooted.treefile",
                       "05e_family_phylogeny_foxYZ/supp_fig_S14.pdf",
                       35, 7)
