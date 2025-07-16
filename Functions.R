library(Rfast)
library(progress)
library(evd)
library(ismev)
library(lattice)
library(RColorBrewer)
library(foreach)
library(doParallel)

# tau function
tf <- function( x )
{
  y        <- x 
  y[y<600] <- log( 1 + exp( y[y<600] ) )
  return( y )
}

# tau inverse function
tf_inv <- function( y )
{
  x        <- y
  x[x<600] <- log( exp(x[x<600]) - 1 )
  return( x )
}

# Transform each cell's marginal distribution to a unit Frechet distribution
transform_Frechet <- function(X, q){
  Y <- colRanks(X)/(nrow(X)+1)
  gpd.mle.mat <- matrix(ncol = 3, nrow = ncol(X))
  gpd.mle.mat[,3] <- apply(X, 2, quantile, probs = q)
  for(i in 1:ncol(X)){
    x <- X[,i]
    u <- gpd.mle.mat[i,3]
    x.trans <- Y[,i]
    gpd.mle.mat[i,1:2] <- gpd.mle <- gpd.fit(x, u, show = FALSE)$mle
    ind <- (x >= u)
    x.trans[ind] <- q + (1-q) * 
      evd::pgpd(x[ind] - u, scale = gpd.mle[1], shape = gpd.mle[2])
    Y[,i] <- (-log(x.trans))^(-1/2) ## to unit Frechet margin
  }
  colnames(gpd.mle.mat) <- c("scale","shape","threshold")
  return(list(Y = Y, MLE = gpd.mle.mat))
}

# compute TPDM 
derive_TPDM <- function( X, q, u=0 )
{
  K <- length( X[1,] )  
  rad <- sqrt( apply( X^2, 1, sum ) )
  if( u == 0 ) 
    u <- quantile( rad, q )
  ind <- which( rad > u )
  ang <- X[ind,]
  for( k in 1:K )
    ang[,k] <- ang[,k] / rad[ind]

  S <- matrix( 0, K, K )
  for(k1 in 1:K){
    for(k2 in 1:K){
      S[k1,k2] <- 1 / length(ind) * sum( ang[,k1] * ang[,k2] ) 
    }
  }
  return( S )
} 

# sample from a mixture of von Mises-Fisher distributions
# Fixed some errors in Directional::rmixvmf
rmixvmf_new <- function (n, probs, mu, k) 
{
  p2 <- c(0, cumsum(probs))
  p <- ncol(mu)
  u <- runif(n)
  g <- length(k)
  ina <- as.numeric(cut(u, breaks = p2))
  ina <- sort(ina)
  nu <- tabulate(ina)
  # y <- array(dim = c(n, p, g))
  x <- matrix(0, nrow = n, ncol = p)
  start.index <- 0
  for (j in 1:length(nu)){
    if(nu[j]>=1){
      x[start.index + 1:nu[j],] <- Rfast::rvmf(nu[j], mu[j,], k[j])
      start.index <- start.index + nu[j]
    }
  }
  list(id = ina, x = x)
}

sample_W <- function( n, m, W, ban, seed=500 )
{
  ## Extract fixed values
  K <- length( W[1,] ) 
  J <- length( W[,1] )
  
  ## Derive values for Z
  Z <- cbind( W[,1:m], -(W[,m+1]<0) + (W[,m+1]>0) )
  for( i in 1:length(W[,1]) )
    Z[i,m+1] <- sqrt( 1 - sum( Z[i,1:m] ^ 2 ) ) * Z[i,m+1]
  
  ## Sample values for lower-dimensional variable Z
  set.seed( seed )
  z_samp <- rmixvmf_new( n = n, probs = rep( 1 / J, J ), 
                     mu = Z, k = rep( 1 / ban^2, J ) )$x
  
  ## Derive nearest neighbour 
  q <- apply( z_samp %*% t(Z), 1, which.max )
  
  ## Map Z to high-dimensional space of W
  w_samp       <- matrix( 0, n, K )
  w_samp[,1:m] <- z_samp[,1:m]
  for( i in 1:n )
    w_samp[i,(m+1):K] <- 
      abs( z_samp[i,(m+1)] / Z[q[i],(m+1)] ) * W[q[i],(m+1):K]
  
  return( w_samp )
}

sample_X_tilde <- function( n, m, W, U, seed=500, q, R, ban, SEASON )
{
  ## Derive fixed values
  K <- length( W[1,] )     
  
  ## Sample radial component ||V||_2 
  ## This is different from the original code provided by Rohrbeck and Cooley (2023)
  set.seed( seed )
  r_samp <- vector(length = n)
  if(SEASON == "SUMMER"){
    high.ind <- (runif(n) >= 1-nrow(W)/121440)#q
  }else{ ## WINTER
    high.ind <- (runif(n) >= 1-nrow(W)/119460)#q
  }
  # above threshold : regular variation
  r_samp[high.ind] <- quantile(R, q)/sqrt(1-runif(sum(high.ind)))
  # below threshold : empirical cdf
  r_samp[!high.ind] <- sample(R[R <= quantile(R, q)], sum(!high.ind), replace = TRUE)
  
  ## Sample values for W
  w_samp <- sample_W( n, m, W, ban=ban, seed )
  
  ## Derive generated samples for V and X_tilde
  v_samp <- w_samp * r_samp
  x_samp <- tf( v_samp %*% t(U) )
  
  return( x_samp )
}

### block bootstrap for the target quantity #3
block_bstrp <- function(X, L){
  # L : block size
  K <- length( X[1,] )        ## Number of variables / gauges
  N <- length( X[,1] )        ## Number of observations (184 x 165 x 4 or 181 x 165 x 4)
  
  if(N %% L != 0) stop("Block size should divide the total length!")
  else{
    blknum <- N / L # number of blocks
    ind <- sample( x=1:blknum, size=blknum, replace = T )
    row_indices <- unlist(lapply(ind, function(i) seq((i - 1) * L + 1, i * L)))
    return(X[row_indices,])
  }
}


### Get point estimates

get_estimates_Q123 = function(X, m, ban, num_run, N, q = 0.9995, SEASON, problem, noise_remove){ # num_run: 생성할 세트 개수 / N: run 1개 길이 
  n = num_run * N
  
  ## Transform observations to Frechet margins
  marginal <- transform_Frechet(X, q = q)
  X_tilde <- marginal$Y
  
  ## Derive TPDM
  Sigma <- derive_TPDM( X_tilde, q = q )
  
  ## Perform eigendecomposition and derive extremal components
  U     <- eigen( Sigma )$vectors
  V     <- t( t( U ) %*% t( tf_inv( X_tilde ) ) )
  
  ### Calculate values for W | ||V||_2 > r0
  R   <- sqrt( apply( V^2, 1, sum ) )
  r0  <- quantile( R, q )
  if(problem %in% c("Q1")){
    ind <- which( R > r0 & apply(X,1,min) > quantile(apply(X,1,min),ifelse(noise_remove,0.99,0)))
  }else if(problem == "Q2"){
    ind <- which( R > r0 & apply(X,1,\(v)sort(v)[20]) > quantile(apply(X,1,\(v)sort(v)[20]),ifelse(noise_remove,0.99,0)))
  }else{
    ind <- which( R > r0 & apply(X,1,\(v)sort(v)[23]) > quantile(apply(X,1,\(v)sort(v)[23]),ifelse(noise_remove,0.99,0)))
  }
  W   <- V[ind,] / R[ind]
  
  ### Generate samples 
  ## set number of cores for parallel computing
  cl <- parallel::makeCluster(min(detectCores()-2, 15), outfile = "")
  registerDoParallel(cl)
  pb <- txtProgressBar(min = 1, max = num_run, style = 3)
  
  X_samp <- foreach(i=1:num_run, .combine = rbind, .export = ls(globalenv()), .packages=c("evd", "Rfast", "ismev")) %dopar% {
    setTxtProgressBar(pb, i)   
    set.seed(i)
    x_samp <- sample_X_tilde( n=N, m=m, W, U, seed=i, q=q, R=R, ban=ban, SEASON=SEASON )
    
    x_samp_f <- evd::pfrechet( x_samp, shape=2 )
    for( k in 1:length( X[1,] )  ){
      id <- which( x_samp_f[,k] <= q )
      
      ## For samples not in the tail
      x_samp_f[id,k] <- quantile( X[,k], x_samp_f[id,k], type=6 )
      
      ## For samples in the tail
      x_samp_f[-id,k] <- marginal$MLE[k,3] + 
        evd::qgpd( (x_samp_f[-id,k]-q) / (1-q), 
                   scale = marginal$MLE[k,1], 
                   shape = marginal$MLE[k,2] )    
    }
    
    if(problem == "Q1"){
      apply(x_samp_f, 1, min)
    }else if(problem == "Q2"){
      apply(x_samp_f, 1, \(v)sort(v)[20])
    }else{
      apply(x_samp_f, 1, \(v)sort(v)[23])
    }
  }
  stopCluster(cl)
  return(X_samp)
}