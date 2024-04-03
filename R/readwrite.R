extractContinuousChar <- function(phy, metadata, col_taxon, cols_contChars){
  data <- read.csv(metadata)
  tree <- ape::read.tree(phy)

  taxa_present <- dplyr::intersect(tree$tip.label, data[[col_taxon]])

  match_tree <- ape::drop.tip(tree, tree$tip.label[-match(taxa_present, tree$tip.label)])
  match_contChars <- data %>%
    dplyr::select(col_taxon, all_of(cols_contChars)) %>%
    dplyr::filter() ## @priscilla: something missing here?

  rownames(contChars) <- data[[col_taxon]]

  return(c(match_tree, match_contChars))
}


