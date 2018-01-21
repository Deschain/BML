library(BML)
col_resic <- bml(BML::col_resic, 5, 0.5, 5)
writeDotFile(col_resic, 'COL_RESIC.dot')