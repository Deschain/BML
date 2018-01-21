library(BML)

lung_sm4 <- bml(BML::lung_sm4,   10, 0.1, 100)

data <- t(lung_sm4$bootstrap$OBS_Probabilities[10:1,])
boxplot(data, horizontal = TRUE, yaxt="n", col="red")
axis(2, at = 1:10, las =2, labels = colnames(data))
title("P(g)")
