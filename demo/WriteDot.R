library(BML)
col_resic <- bml(col_resic_data, 5, 0.5, 5)
writeDotFile(col_resic, 'COL_RESIC.dot')