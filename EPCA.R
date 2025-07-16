rm(list = ls())
## Modify lines 5 to 8. Use `m-bandwidth-select.R` to get optimal values of `m` and `ban`.
## SUMMER: Mar.--Oct. / WINTER: Nov.--Apr.

SEASON  = c("SUMMER","WINTER")[1]
problem = c("Q1","Q2","Q3")[1]
m       = 7    # num. of PCs    
ban     = 0.3  # bandwidth

noise_remove = TRUE
num_run = 1000    # number of runs
B       = 199     # number of bootstraps

## load dataset
X = readRDS(paste("Final",SEASON,"rds",sep="."))

## Extract dimensions
K <- length( X[1,] )        ## Number of variables / gauges
N <- length(X[,1])/4        ## Length of each run

## load functions
source("Functions.R")

################################################################################
## 1. Point estimate
# (1) Save synthetic dataset
X_samp = get_estimates_Q123(X, m, ban, num_run, N, q=0.9995,
                            SEASON=SEASON, problem=problem, noise_remove=noise_remove)
saveRDS(X_samp, paste("point",problem,SEASON,"rds",sep="."))

# (2) Point estimate for each season. To get a correct CI, use `Confidence-interval.R` file.
X_samp = readRDS(paste("point",problem,SEASON,"rds",sep="."))
if(problem == "Q1"){
  paste0(problem," estimate: ",sum(X_samp >= 1.7)/num_run) |> cat()
}else if(problem == "Q2"){
  paste0(problem," estimate: ",sum(X_samp >= 5.7)/num_run) |> cat()
}else{
  Q3.estimate.daily <- sum(X_samp >= 5)/num_run
  paste0(problem," estimate (daily exceedance): ",Q3.estimate.daily, "\n") |> cat()

  X.23 <- apply(X,1,\(v) sort(v)[23])
  chi.hat <- extRemes::taildep(X.23[-ifelse(SEASON=="SUMMER",121440,119460)], X.23[-1], type = "chi", u = ecdf(X.23)(5))
  if(chi.hat == 0){
    paste0(problem," estimate (consecutive exceedance): ",(Q3.estimate.daily)^2/ifelse(SEASON=="SUMMER",121440,119460)*4) |> cat()
  }else{
    paste0(problem," estimate (consecutive exceedance): ",Q3.estimate.daily*chi.hat) |> cat()
  }
}


################################################################################
## 2. Confidence Interval
## (1) Save point estimates from bootstrapped datasets
ANSWERS_B = c()
index = 1
if(problem %in% c("Q1","Q2")){
  for(b in 1:(10*B)){ 
    set.seed(b)
    X_b = X[sample(1:nrow(X),nrow(X), replace = T),]
    skip_to_next = FALSE
    X_samp_b = tryCatch(get_estimates_Q123(X_b, m, ban, num_run, N, q=0.9995,
                                           SEASON=SEASON, problem=problem, noise_remove=noise_remove), error = \(e){skip_to_next <<- TRUE})
    if(skip_to_next){ next }
    if(problem == "Q1"){
      ANSWERS_B[index] = sum(X_samp_b >= 1.7)/num_run
    } else if(problem == "Q2"){
      ANSWERS_B[index] = sum(X_samp_b >= 5.7)/num_run
    }
    print(c(index,b))
    if(index==B){break}
    index = index + 1
  }
}else{
  L = ifelse(SEASON=="SUMMER",184,181)
  for(b in 1:(10*B)){ 
    set.seed(b)
    X_b = block_bstrp(X, L)
    skip_to_next = FALSE
    X_samp_b = tryCatch(get_estimates_Q123(X_b, m, ban, num_run, N, q=0.9995,
                                           SEASON=SEASON, problem=problem, noise_remove=noise_remove), error = \(e){skip_to_next <<- TRUE})
    if(skip_to_next){ next }
    Q3.estimate.daily <- sum(X_samp_b >= 5)/num_run
    X_b.23 = apply(X_b, 1, \(v) sort(v)[23])
    chi.hat <- extRemes::taildep(X_b.23[-ifelse(SEASON=="SUMMER",121440,119460)], X_b.23[-1], type = "chi", u = ecdf(X_b.23)(5))
    if(is.nan(chi.hat) | chi.hat == 0){
      ANSWERS_B[index] = (Q3.estimate.daily)^2/ifelse(SEASON=="SUMMER",121440,119460)*4
    }else{
      ANSWERS_B[index] = Q3.estimate.daily*chi.hat
    }
    print(c(index,b))
    if(index==B){break}
    index = index + 1
  }
}
saveRDS(ANSWERS_B, paste("boots",problem,SEASON,"rds",sep="."))


## (2) CI for each season. To get a correct CI, use `Confidence-interval.R` file.
ANSWERS_B = readRDS(paste("boots",problem,SEASON,"rds",sep="."))
X_samp = readRDS(paste("point",problem,SEASON,"rds",sep="."))
if(problem == "Q1"){
  est = sum(X_samp >= 1.7)/num_run
}else if(problem == "Q2"){
  est  = sum(X_samp >= 5.7)/num_run
}else{
  Q3.estimate.daily <- sum(X_samp >= 5)/num_run
  X.23 <- apply(X,1,\(v) sort(v)[23])
  chi.hat <- extRemes::taildep(X.23[-ifelse(SEASON=="SUMMER",121440,119460)], X.23[-1], type = "chi", u = ecdf(X.23)(5))
  if(chi.hat == 0){
    est = (Q3.estimate.daily)^2/ifelse(SEASON=="SUMMER",121440,119460)*4
  }else{
    est = Q3.estimate.daily * chi.hat
  }
}
ci = quantile(ANSWERS_B, c(0.025, 0.975))
paste0(problem," ",SEASON," confidence interval: (", 2*est - ci[2], ",",2*est-ci[1],")") |> cat()



