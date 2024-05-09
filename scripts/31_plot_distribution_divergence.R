suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))

# define the function to compute genetic distances, given an alignment file and a substitution model
computeGeneticDistances <- function(alignment, model) {
  # read in the file
  alignment <- phangorn::read.phyDat(alignment, format = "fasta", type = "AA")
  
  # compute the genetic distance matrix
  dist <- phangorn::dist.ml(alignment, model = model)
  
  # transform the matrix into a list
  dist.list <- dist[1:length(dist)]
  
  # return the list
  return(dist.list)
}


#######################
#     Actual code     #
#######################

write("", stdout())
write("Reading in the file with substitution models...", stdout())

# read in the tsv file with selected models per each alignment
modelData <- read.table("13_distribution_divergence/models_perOrthogroup_Rformatted.tsv", header = TRUE, sep = "\t")

# create a named list from dataframe
modelList <- with(modelData, split(model, alignment))

distList <- list()

write("Computing genetic distances per each orthogroup, this may take a while...", stdout())

# compute genetic distances for each alignment and store values in a named list
for (alignment in names(modelList)) {
  distList[[alignment]] <- computeGeneticDistances(alignment, modelList[[alignment]])
}

write("	Done!", stdout())


write("Computing median values of distances per each orthogroup...", stdout())
# create a dataframe with median values of each alignment
distMedian_values <- stack(distList) %>%
  group_by(ind) %>%
  summarise(median = median(values)) %>%
  rename(gene = ind)

# compute density function of median values
density.median <- density(distMedian_values$median)
density.median.df <- data.frame(x = density.median$x, y = density.median$y)

write("Computing quantile values...", stdout())
# compute the 95% quantile
quantile = c(0.90, 0.95, 0.99)
density.median.quant <- quantile(distMedian_values$median, prob = quantile)
density.median.df$quant <- factor(findInterval(density.median.df$x, density.median.quant))

write("Writing the median values per each orthogroup to file...", stdout())
# write median table to file
distMedian_values$quant <- factor(findInterval(distMedian_values$median, density.median.quant))
write.table(arrange(distMedian_values, median),
            file = "13_distribution_divergence/distance_median_values.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)

write("Plotting results...", stdout())
# plot the distribution
density.plot <- ggplot(density.median.df, aes(x, y)) +
  geom_ribbon(aes(ymin = 0, ymax = y, fill = quant)) +
  geom_line(linewidth = 1.5) +
  scale_fill_brewer(palette = "RdPu",
                    labels = c("", paste0("tail ", quantile*100, "%"))) +
  xlab("Genetic distance (aa)") +
  ylab("Density") +
  ggtitle(paste0("Distribution of amino acid divergence of ", nrow(distMedian_values), " genes")) +
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        panel.background = element_rect(fill = 'transparent'),
        plot.background = element_rect(fill = 'transparent', color = NA),
        axis.line = element_line(color = "grey50", linewidth = 0.8),
        axis.ticks = element_line(color = "grey50", linewidth = 0.8),
        axis.text = element_text(color = "grey50"),
        # legend.position = "none",
        legend.title = element_blank(),
        legend.background = element_rect(fill = 'transparent', color = NA),
        legend.text = element_text(color = "grey50"))

density.plot

ggsave("13_distribution_divergence/density_plot.png",
       plot = density.plot, device = "png",
       dpi = 300, height = 8, width = 12, units = ("in"), bg = 'transparent')

ggsave("13_distribution_divergence/density_plot.pdf",
       plot = density.plot, device = "pdf",
       dpi = 300, height = 8, width = 12, units = ("in"), bg = 'transparent')

write("DONE!", stdout())
