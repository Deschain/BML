library(BML)
library(graph)

gmb_sm4 <- bml(BML::gbm_sm4, 20, 0.1)
gmb_sm4_dag <- new("graphNEL", nodes=gmb_sm4$DAG$labels, edgemode="directed")

for (edge_idx in 1:length(gmb_sm4$DAG$edges_1)) {

  gmb_sm4_dag <- addEdge(gmb_sm4$DAG$labels[match(gmb_sm4$DAG$edges_1[edge_idx], gmb_sm4$DAG$nodes)],
                        gmb_sm4$DAG$labels[match(gmb_sm4$DAG$edges_2[edge_idx], gmb_sm4$DAG$nodes)],
                        gmb_sm4_dag, 1)
}

plot(gmb_sm4_dag)
