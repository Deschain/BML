
bml <- function(dataset, ntree, threshold, rep = 0) {
  
  #Asserts
  if (ntree <= 0) {
    stop("Ntree can't be less or equal 0")
  }
  
  if (threshold < 0.0 | threshold > 1.0){
    stop("Threshold must be between 0.0 and 1.0")
  }
  
  if (!(is.matrix(dataset))){
    stop("Data must be a matrix class")
  }
  
  if (any(is.na(dataset))){
    stop("Data can't have NAs values")
  }
  
  if (!(all(dataset == 0 | dataset == 1))){
    stop("Only 0 or 1 values are accepted (mutated or not)")
  }
  
  if(is.null(colnames(dataset))) {
    stop("Data must have gene names as column names")
  }
  
  bml <- BML(dataset, ntree, threshold, rep)
  # Maybe store input data aswell
  
  class(bml) <- append(class(bml),"bml")
  return(bml)
}


print.bml <- function(bml, ...) {
  cat("Total Num of Edges ", bml$num_edges, "\n")
  cat("Num of unpruned edges ", bml$num_unpruned_edges, "\n")
  cat("Total number of genes ", bml$num_genes, "\n")
  cat("Number of genes with no parent after global pruning ", bml$num_noparent_after_global_pruning, "\n")
}

summary.bml <- function(bml, ...) {
  print(bml)
}

adjacency_matrix <- function(bml) {
  UseMethod("adjacency_matrix", bml)
}


adjacency_matrix.default <- function(data) {
  cat(paste("Adjancency matrix is not implemented for class", class(data)))
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


writeDotFile.bml <- function(data, file) {
  writeDotFile(data, file)
}
