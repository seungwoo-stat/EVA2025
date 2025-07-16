problem = c("Q1","Q2","Q3")[1]

num_run = 1000
est = 0
for(SEASON in c("SUMMER","WINTER")){
  X = readRDS(paste0("Final.",SEASON,".rds"))
  X_samp = readRDS(paste0("results/point.",problem,".",SEASON,".rds")) #For point estimate
  if(problem == "Q1"){
    est = est + sum(X_samp >= 1.7)/num_run
  }else if(problem == "Q2"){
    est = est + sum(X_samp >= 5.7)/num_run
  }else{
    Q3.estimate.daily <- sum(X_samp >= 5)/num_run
    X.23 <- apply(X,1,\(v) sort(v)[23])
    chi.hat <- extRemes::taildep(X.23[-ifelse(SEASON=="SUMMER",121440,119460)], X.23[-1], type = "chi", u = ecdf(X.23)(5))
    if(chi.hat == 0){
      est = est + (Q3.estimate.daily)^2/ifelse(SEASON=="SUMMER",121440,119460)*4
    }else{
      est = est + Q3.estimate.daily*chi.hat
    }
  }
}

ANSWERS_B_SUMMER = readRDS(paste0("results/boots.",problem, ".SUMMER.rds"))
ANSWERS_B_WINTER = readRDS(paste0("results/boots.",problem, ".WINTER.rds"))
ci.pair = sort(ANSWERS_B_SUMMER + ANSWERS_B_WINTER)[c(5,195)]

paste0(problem," Total estimate: ", est,"\n") |> cat()
paste0(problem," Total confidence interval from 199 pairs: (", 2*est - ci.pair[2], ",",2*est-ci.pair[1],")\n") |> cat()