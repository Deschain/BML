bml <- function(file, ntree, pthres) {
  data <- read.table(file, header = TRUE, row.names=1, skip = 1)
  stopifnot(ntree > 0)
  stopifnot(pthres > 0.0, pthres<1.0)
  # Add sanity checks
  bml <- BML(as.matrix(data), ntree, pthres)
  attr(bml,'class') <- 'bml'
  return(bml)
}

print.bml <- function(bml){
  cat("Total Num of Edges ", bml$num_edges, "\n")
  cat("Num of unpruned edges ", bml$num_unpruned_edges, "\n")
  cat("Total number of genes ", bml$num_genes, "\n")
  cat("Number of genes with no parent after global pruning ", bml$num_noparent_after_global_pruning, "\n")
}

adm <- function(x){
  UseMethod("adm", x)
}

adm.bml <- function(bml) {
  tmp <- cbind(match(bml$DAG$edges_1, aux$DAG$nodes), match(bml$DAG$edges_2, aux$DAG$nodes))
  adm <- matrix(0, length(bml$DAG$nodes), length(bml$DAG$nodes))
  adm[tmp] <- 1
  colnames(adm) <- bml$DAG$labels
  rownames(adm) <- bml$DAG$labels
  return(adm)
}
