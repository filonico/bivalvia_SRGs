#!/usr/bin/env Rscript

library(dplyr)
library(topGO)

args = commandArgs(trailingOnly = TRUE)

# args1 gene universe
# args2 genes of interest
# args3 output prefix


######################
#     READ INPUT     #
######################

universe_GOannotation <- read.table(file = args[1], skip = 4, header = FALSE)
interest_GOannotation <- read.table(file = args[2], header = FALSE, sep = "\t")


#####################
#     FUNCTIONS     #
#####################


# function to perform GO enrichment and generate summary tables
performGOEnrichment <- function(GOannotation, ontology, algorithm) {
  
  # create topGOdata object
  GOdata <- new("topGOdata", ontology = ontology, allGenes = geneList, annot = annFUN.gene2GO, gene2GO = GOannotation)
  
  # perform GO enrichment
  result <- runTest(GOdata, algorithm = algorithm, statistic = "fisher")
  
  # get summary table
  results <- GenTable(GOdata, classicFisher = result, ranksOf = "classicFisher", topNodes = 500)
  
  return(results)
}


############################
#     DATA PREPARATION     #
############################

# prepare the gene universe dataset as a named list
geneID2GO <- universe_GOannotation %>% dplyr::select(c(V2, V4)) %>%
  group_by(V2) %>%
  summarize(V4 = paste(V4, collapse = ","), .groups = 'drop') %>%
  mutate(V4 = strsplit(V4, ",")) %>%
  tibble::deframe()
  
# get list of genes
geneNames <- names(geneID2GO)

# get the vector of the genes of interest
gene_int_list <- as.vector(interest_GOannotation$V2)

# get a factor list of interesting/not interesting genes 
geneList <- factor(as.integer(geneNames %in% gene_int_list))
names(geneList) <- geneNames


#########################
#     GO ENRICHMENT     #
#########################

# define GO ontologies
ontologies <- c("BP", "MF", "CC")

# define algorithms
algorithms <- c("classic", "elim")

# perform GO enrichment and generate summary tables for each ontology and each method
for (ontology in ontologies) {
  
  for (algorithm in algorithms) {
    
    result <- performGOEnrichment(geneID2GO, ontology, algorithm)
    file_name <- paste0(args[3], "topGO_", ontology, "_", algorithm, ".txt")
    write.table(result, file = file_name, quote = FALSE, row.names = FALSE, sep = "\t")
    
  }
}
