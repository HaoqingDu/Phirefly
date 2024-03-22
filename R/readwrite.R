library(ape)
library(phytools)
library(tidyverse)

extractContinuousChar <- function(phy, metadata, col_taxon, cols_contChars){
  data <- read.csv(metadata)
  tree <- read.tree(phy)

  taxa_present <- intersect(tree$tip.label, data[[col_taxon]])



  match_tree <- drop.tip(tree, tree$tip.label[-match(taxa_present, tree$tip.label)])
  match_contChars <- data %>%
    select(col_taxon, all_of(cols_contChars)) %>%
    filter()

  rownames(contChars) <- data[[col_taxon]]

  return(c(match_tree, match_contChars))
}


