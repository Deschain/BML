library(BML)
Lung_SM4 <- bml('demo/Data/Lung_SM4',   10, 0.1, 5)
tPositives <- as.data.frame(Lung_SM4$bootstrap$ConfidenceTruePositives)

plot(tPositives$`Tree+DAG`, tPositives$Confidence, type = 'n', xlab = '', ylab = '')
title("Lung_SM4 True Positives")
lines(tPositives$`Tree+DAG`, tPositives$Confidence, type = 'o', col = 'blue', lty = 1, pch = 1)
lines(tPositives$DAG, tPositives$Confidence, type = 'o', col = 'red', lty = 2, pch = 0)