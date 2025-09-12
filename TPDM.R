source("Functions.R")

set.seed(0)
X = readRDS("SUMMER.rds")
marginal <- transform_Frechet(X, q = 0.9995)
X_tilde.s <- marginal$Y
matrix_data.s <- derive_TPDM( X_tilde.s, q=0.9995 )
zlim <- range(matrix_data.s, finite=TRUE)

## Derive TPDM
Sigma <- matrix_data.s

par(mar = c(5, 5, 2, 5)) 
layout(matrix(c(1,2), nrow=1), widths=c(4.8,0.8))  

image(1:nrow(Sigma), 1:ncol(Sigma), Sigma, 
      col = colorRampPalette(c("white", "red", "darkred"))(100), xaxt = "n", yaxt = "n", 
      xlab="Cell", ylab="Cell")

## 5x5 grids (generate 25 labels)
x_labels <- as.vector(outer(1:5, 1:5, FUN = function(i, j) paste0("(", i, ",", j, ")")))
y_labels <- x_labels  

axis(1, at = seq(1, nrow(matrix_data.s), length.out = 25), labels = x_labels, las = 2, cex.axis = 0.7)
axis(2, at = seq(1, ncol(matrix_data.s), length.out = 25), labels = y_labels, las = 2, cex.axis = 0.7)

## Bar
library(fields)
par(mar = c(6,0,3,0)) 
image.plot(zlim = zlim, col = colorRampPalette(c("white", "red", "darkred"))(100),
           legend.only = TRUE,
           legend.width = 0.8, legend.shrink = 1.0,
           axis.args   = list(cex.axis = 0.85, las = 1, tck = -0.5),
           legend.args = list(text = "", cex = 0.7))

