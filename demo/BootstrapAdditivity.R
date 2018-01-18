library(BML)
Lung_SM4 <- bml(lung_sm4_data,   10, 0.1, 100)

data <- as.data.frame(t(log(Lung_SM4$bootstrap$EdgeProbabilities[(9+3):9,])))
boxplot(data, horizontal = TRUE, col = rev(c('red', 'lightblue', 'red','grey')), yaxt='n')
axis(2, at =1:4, labels = names(data), las = 2)
title("Log(P(g))")
