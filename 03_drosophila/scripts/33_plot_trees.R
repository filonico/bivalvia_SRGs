#!/usr/bin/env Rscript

library(ggtree)
library(ggplot2)

####################
#     FUNCTION     #
####################

# function to load, plot, and save tree
plot.tree <- function(input_filename, output_filename, output_height, output_width) {
  
  tree <- read.tree(file = input_filename)

  tree$tip.label <- ifelse(grepl("Dhyd|Dpse|Dsuz|Dmel", tree$tip.label),
                           paste0("*", tree$tip.label),
                           tree$tip.label)
  
  tree.plot <- ggtree(tree, size = 0.2) +
    
    geom_tiplab(size = 1, offset = 0.05) +
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

dmrt.tree <- plot.tree("06_possvm_orthology/dmrt_ALL_aligned_trim04.faa.treefile.ortholog_groups.newick",
                       "06_possvm_orthology/02_plot_trees/supp_fig_S8.pdf",
                       5, 12)

sox.tree <- plot.tree("06_possvm_orthology/sox_ALL_aligned_trim04.faa_rooted.treefile.ortholog_groups.newick",
                      "06_possvm_orthology/02_plot_trees/supp_fig_S9.pdf",
                      7, 12)

fox.tree <- plot.tree("06_possvm_orthology/fox_ALL_aligned_trim04.faa_rooted.treefile.ortholog_groups.newick",
                      "06_possvm_orthology/02_plot_trees/supp_fig_S10.pdf",
                      15, 15)
