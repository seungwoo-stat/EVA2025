rm(list = ls())

source("Functions.R")

X1 = readRDS("SUMMER.rds")
X2 = readRDS("WINTER.rds")
X <- rbind(X1,X2)

X.Frechet.list <- transform_Frechet(X, q = 0.9995)
X.tilde        <- X.Frechet.list[[1]]

### Estimate the TPDM and derive eigendecomposition
Sigma <- derive_TPDM(X.tilde, q=0.9995)
estim <- eigen(Sigma)
U     <- estim$vectors
D     <- estim$values

### Derive extremal principal components 
V   <- t(t(U) %*% t(tf_inv(X.tilde)))

### Calculate values for W | ||V||_2 > r0
R   <- sqrt(rowSums(V^2))
r0  <- quantile(R, 0.9995)

X.min <- apply(X,1,min)
X.sixth <- apply(X,1,\(v)sort(v)[20])
X.third <- apply(X,1,\(v)sort(v)[23])

par(mfrow=c(1,3), mar=c(5,2,4,2)+0.1, oma=c(0,4,0,0))  
plot(R, X.min, xlab = "||V||", main = "(a) Minimum of all 25 cells", pch = 20, cex.lab=1.5, cex.main = 1.5)
abline(v = quantile(R, 0.9995))
abline(h = quantile(X.min, 0.99), col = "red")

plot(R, X.sixth, xlab = "||V||", main = "(b) Sixth largest among 25 cells", pch = 20, cex.lab=1.5, cex.main = 1.5)
abline(v = quantile(R, 0.9995))
abline(h = quantile(X.sixth, 0.99), col = "red")

plot(R, X.third, xlab = "||V||", main = "(c) Third largest among 25 cells", pch = 20, cex.lab=1.5, cex.main = 1.5)
abline(v = quantile(R, 0.9995))
abline(h = quantile(X.third, 0.99), col = "red")
mtext("Daily statistics (Leadbetters)", side=2, outer=TRUE, line=1)
