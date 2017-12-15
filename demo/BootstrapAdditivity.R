library(BML)
Lung_SM4 <- bml('demo/Data/Lung_SM4',   10, 0.1, 5)

data <- as.data.frame(t(log(Lung_SM4$bootstrap$EdgeProbabilities[9:(9+3),])))
boxplot(data, horizontal = TRUE, col = c('red', 'lightblue', 'red','grey'), yaxt='n')
axis(2, at =1:4, labels = names(data), las = 2)
