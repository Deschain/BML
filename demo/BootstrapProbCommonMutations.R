library(BML)
Lung_SM4 <- bml(lung_sm4_data,   10, 0.1, 100)

data <- t(Lung_SM4$bootstrap$OBS_Probabilities[10:1,])
boxplot(data, horizontal = TRUE, yaxt="n", col="red")
axis(2, at = 1:10, las =2, labels = colnames(data))
title("P(g)")
