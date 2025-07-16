## choose m (# of PCs) and bandwidth (in spherical KDE)

rm(list = ls())

## load functions
source("Functions.R")

## load datasets
Final.SUMMER  <- readRDS("SUMMER.rds")
Final.WINTER  <- readRDS("WINTER.rds")

################################################################################
## iterate the code below for the combination of {SUMMER, WINTER} and {Q1,Q2,Q3}
X <- Final.SUMMER
SEASON <- "SUMMER"
problem <- c("Q1","Q2","Q3")[1]

noise_remove <- TRUE
################################################################################

## transform each column of X to follow Frechet(loc = 0, scale = 1, shape = 2)
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
if(problem == "Q1") X.problem <- apply(X,1,min)
if(problem == "Q2") X.problem <- apply(X,1,\(v)sort(v)[20])
if(problem == "Q3") X.problem <- apply(X,1,\(v)sort(v)[23])
ind <- which(R > r0 & X.problem > quantile(X.problem, ifelse(noise_remove, 0.99, 0)))
W   <- V[ind,] / R[ind]

## m and bandwidth selection
## set a dense grid of possible values of m and bandwidth
bandwidth.seq <- seq(0.05,0.7,0.05)##seq(0.05,0.3,0.01)
m.seq <- c(3:15)##3:24

## Number of cores used for parallel computing
cl <- parallel::makeCluster(6, outfile = "")
registerDoParallel(cl)
pb <- txtProgressBar(min = 1, max = length(bandwidth.seq)*length(m.seq), style = 3)
bw.m.select <- foreach(i=seq_along(bandwidth.seq), .combine = rbind, .export = ls(globalenv()), .packages=c("evd", "Rfast", "ismev")) %dopar% {
  bw.m.select.vec <- vector(length = length(m.seq))
  for(j in seq_along(m.seq)){
    setTxtProgressBar(pb, length(m.seq)*(i-1) + j) 
    X.tilde.samp <- sample_X_tilde(n=2000000, m=m.seq[j], W, U, seed=300,
                                   q=0.9995, R=R, ban=bandwidth.seq[i], SEASON=SEASON)
    bw.m.select.vec[j] <- sum(abs(colMeans(X.tilde.samp >= qfrechet(1-nrow(W)/nrow(X),shape=2)) - (nrow(W)/nrow(X))))
  }
  bw.m.select.vec
}
stopCluster(cl)

rownames(bw.m.select) <- bandwidth.seq
colnames(bw.m.select) <- m.seq

bw.m.select
bw.m.select == min(bw.m.select)
par(mfrow = c(1,1))
image(bw.m.select, xaxt ="n", yaxt ="n", xlab ="bandwidth", ylab ="m",
      main=paste(SEASON,problem))
axis(1, seq(0,1,length.out=nrow(bw.m.select)), rownames(bw.m.select))
axis(2, seq(0,1,length.out=ncol(bw.m.select)), colnames(bw.m.select))

## select m and bandwidth with smallest value in the `bw.m.select` matrix


