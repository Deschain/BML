library(BML)
library(graph)

gmb_sm4 <- bml(BML::gbm_sm4, 20, 0.1)
gmb_sm4_graph <- graphAM(adjacency_matrix(gmb_sm4), edgemode = "directed")
plot(gmb_sm4_graph)
