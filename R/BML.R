
bml <- function(data, ntree, threshold, rep = 0) {
  stopifnot(ntree > 0)
  stopifnot(threshold > 0.0, threshold<1.0)
  # Add sanity checks
  bml <- BML(as.matrix(data), ntree, threshold, rep)
  # Maybe store input data aswell
  class(bml) <- append(class(bml),"bml")
  return(bml)
}


print.bml <- function(bml) {
  cat("Total Num of Edges ", bml$num_edges, "\n")
  cat("Num of unpruned edges ", bml$num_unpruned_edges, "\n")
  cat("Total number of genes ", bml$num_genes, "\n")
  cat("Number of genes with no parent after global pruning ", bml$num_noparent_after_global_pruning, "\n")
}

summary.bml <- function(bml) {
  print(bml)
}

adjacency_matrix <- function(data) {
  UseMethod("adjacency_matrix", data)
}


adjacency_matrix.default <- function(data) {
  print("Adjancency matrix not implemented for class")
}


adjacency_matrix.bml <- function(bml) {
  tmp <- cbind(match(bml$DAG$edges_1, bml$DAG$nodes), match(bml$DAG$edges_2, bml$DAG$nodes))
  adm <- matrix(0, length(bml$DAG$nodes), length(bml$DAG$nodes))
  adm[tmp] <- 1
  colnames(adm) <- bml$DAG$labels
  rownames(adm) <- bml$DAG$labels
  return(adm)
}


writeDotFile <- function(data, file) {
  UseMethod("writeDotFile", data, file)
}


writeDotFile.bml <- function(data, file)
{
  writeDotFile(data, file)
}
