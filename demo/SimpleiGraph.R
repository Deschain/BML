library(BML)
library(igraph)

gmb_sm4 <- bml(BML::gbm_sm4, 20, 0.1)
gmb_sm4_dag <- graph_from_adjacency_matrix(adjacency_matrix(gmb_sm4), mode = "directed")
plot(gmb_sm4_dag)
plot(gmb_sm4_dag, layout= layout_as_tree(gmb_sm4_dag, flip.y = TRUE))
