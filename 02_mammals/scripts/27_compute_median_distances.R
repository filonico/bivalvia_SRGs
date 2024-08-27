#!/usr/bin/env Rscript

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

# write median table to file
write("Writing the median values per each orthogroup to file...", stdout())

write.table(arrange(distMedian_values, median),
            file = "13_distribution_divergence/distance_median_values.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)
