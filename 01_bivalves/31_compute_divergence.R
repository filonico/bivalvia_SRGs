#!/usr/bin/env Rscript

suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))


#################
#     INPUT     #
#################

models_perOrthogroup_Rformatted <- "13a_distribution_divergence_onlyCrassostrea/models_perOrthogroup_Rformatted.tsv"


#################
#     OUTPUT    #
#################

distance_median_values <- "13a_distribution_divergence_onlyCrassostrea/distance_median_values.tsv"


#####################
#     FUNCTIONS     #
#####################

# define the function to compute genetic distances, given an alignment file and a substitution model
computeGeneticDistances <- function(alignment, model) {
  # read in the file
  alignment <- phangorn::read.phyDat(alignment, format = "fasta", type = "AA")

  # get the alignment length
  length <- nrow(as.data.frame(alignment))

  # get the number of species
  species <- length(alignment)

  # compute the genetic distance matrix
  dist <- phangorn::dist.ml(alignment, model = model)

  # create a list
  dist.list <- list("length" = length,
                    "species" = species,
                    "dist" = dist[1:length(dist)])

  # return the list
  return(dist.list)

}


######################################
#     COMPUTE PAIRWISE DISTANCES     #
######################################

write("", stdout())
write("Reading in the file with substitution models...", stdout())

# read in the tsv file with selected models per each alignment
modelData <- read.table(models_perOrthogroup_Rformatted, header = TRUE, sep = "\t")

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
dist_df <- do.call(rbind, lapply(names(distList), function(x) {
  # extract length, species, and dist from each list element
  length <- distList[[x]]$length
  species <- distList[[x]]$species
  dist <- distList[[x]]$dist
  
  # calculate the median of dist
  median <- median(dist)
  
  # create a data frame for each alignment
  data.frame(gene = x, length = length, species = species, median = median)
}))

# write median table to file
write("Writing the median values per each orthogroup to file...", stdout())

write.table(arrange(dist_df, median),
            file = distance_median_values,
            sep = "\t", quote = FALSE, row.names = FALSE)
