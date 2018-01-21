library("BML")

lung_sm4 <- bml(BML::lung_sm4,   10, 0.1, 5)
tPositives <- as.data.frame(lung_sm4$bootstrap$ConfidenceTruePositives)

plot(tPositives$`Tree+DAG`, tPositives$Confidence, type = 'n', xlab = '', ylab = '', yaxt = 'n')
axis(2, at = tPositives$Confidence, las = 2)
# axis(3, at = tPositives$`Tree+DAG`)
title("Lung_SM4 True Positives")
lines(tPositives$`Tree+DAG`, tPositives$Confidence, type = 'o', col = 'blue', lty = 1, pch = 1)
lines(tPositives$DAG, tPositives$Confidence, type = 'o', col = 'red', lty = 2, pch = 0)
