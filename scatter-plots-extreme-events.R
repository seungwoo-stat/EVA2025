source("Functions.R")

## load dataset
X = readRDS("SUMMER.rds")
SEASON = "SUMMER"
m       = 22     # num. of PCs    
ban     = 0.15  # bandwidth
noise_remove = TRUE

## extract dimensions
K <- length(X[1,])        ## num. of variables
N <- length(X[,1])/4      ## length of one run

## Transform observations to Frechet dist.
set.seed(0)
q=0.9995
marginal <- transform_Frechet(X, q = q)
X_tilde <- marginal$Y

## Derive TPDM
Sigma <- derive_TPDM( X_tilde, q = q )

## Extremal principal components
U     <- eigen( Sigma )$vectors
V     <- t( t( U ) %*% t( tf_inv( X_tilde ) ) )

### Calculate values for W | ||V||_2 > r_0
R   <- sqrt( apply( V^2, 1, sum ) )
r0  <- quantile( R, 0.99 )
ind <- which( R > r0 & apply(X,1,min) > quantile(apply(X,1,min),ifelse(noise_remove,0.99,0)))
W   <- V[ind,] / R[ind]

num_run = 1
X_tilde_samp <- sample_X_tilde( n=N*num_run, m=m, W, U, seed=500, q, R, ban, SEASON )

par(mfrow=c(1,2))

#Adjacent pair: (1,1) and (1,2)
k1=1
k2=6

## Plot synthetic points
ind <- which( sqrt( apply(X_tilde_samp^2, 1, sum) ) > r0 )

plot( X_tilde_samp[ind,k1], X_tilde_samp[ind,k2], 
      xlab=paste("Cell ","(1,1)"), ylab=paste("Cell ","(1,2)"), xlim=c(0,850), ylim=c(0,850), cex.lab=1,
      cex.axis=1, cex=1, col=8, pch=19 , main="(a) Adjacent pair")

## Add original points
ind <- which( sqrt( apply(X_tilde^2, 1, sum) ) > r0 )
points( X_tilde[ind,k1], X_tilde[ind,k2], col = alpha(1, 1), pch=19, cex=0.8  )
abline(0,1)
legend("topleft",legend=c("Original","Generated"),col=c("black","grey"), pch=16, cex=1.0)


#Distant pair: (1,1) and (3,3)
k1=1
k2=13

## Plot synthetic points
ind <- which( sqrt( apply(X_tilde_samp^2, 1, sum) ) > r0 )

plot( X_tilde_samp[ind,k1], X_tilde_samp[ind,k2], 
      xlab=paste("Cell ","(1,1)"), ylab=paste("Cell ","(3,3)"), xlim=c(0,850), ylim=c(0,850), cex.lab=1,
      cex.axis=1, cex=1, col=8, pch=19, main = "(b) Distant pair" )

## Add Original points
ind <- which( sqrt( apply(X_tilde^2, 1, sum) ) > r0 )
points( X_tilde[ind,k1], X_tilde[ind,k2], col = alpha(1, 1), pch=19, cex=0.8  )
abline(0,1)
legend("topleft",legend=c("Original","Generated"),col=c("black","grey"), pch=16, cex=1.0)

