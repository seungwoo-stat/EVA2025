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

par(mfrow= c(1,3))
plot(R, X.min, xlab = "||V||", ylab = "Minimum of the 25 cells", pch = 20)
abline(v = quantile(R, 0.9995))
abline(h = quantile(X.min, 0.99), col = "red")

plot(R, X.sixth, xlab = "||V||", ylab = "Sixth largest of the 25 cells", pch = 20)
abline(v = quantile(R, 0.9995))
abline(h = quantile(X.sixth, 0.99), col = "red")

plot(R, X.third, xlab = "||V||", ylab = "Third largest of the 25 cells", pch = 20)
abline(v = quantile(R, 0.9995))
abline(h = quantile(X.third, 0.99), col = "red")