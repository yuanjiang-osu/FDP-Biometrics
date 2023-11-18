############# Simulation Settings ######################

# Specify model M1-M7 here
sigma <- c("independence", "equal", "fan", "Cauchy", "3f", "2f", "nonlinear")[1];

# Specify the number of replicates here
REP.N <- 1000;

########################################################

### Library for parallel computing
if(!require(snow)) { install.packages("snow"); library(snow); }


### Main function
evaluation <- function(input.list, param){
  
  ### libraries needed by RcPP
  if(!require(gtools)) { install.packages("gtools"); library(gtools); }
  if(!require(Rcpp)) { install.packages("Rcpp"); library(Rcpp); }
  if(!require(devtools)) { install.packages("devtools");  libray(devtools); }
  
  library(MASS)
  library(pbivnorm)
  library(qvalue)
  library(limma)
  
  sourceCpp("CPP_Eval_3.cpp");
  source("SLIMcodes.R")
  
  sigma = param$sigma;
  method = param$method;
  alpha = param$alpha;
  nzero.num.true = param$nzero.num.true;
  n.sim = param$n.sim;
  
  set.seed(input.list$seed.input);
  
  n <- 400
  p <- 2000
  
  mu <- read.csv(file = paste("2000 400", nzero.num.true, "mu.csv"))$x
  mu[mu != 0] <- 2 * qnorm(1 - alpha/2)
  
  if(sigma == "independence") {
    Sigma <- as.matrix(diag(p))
  } else {
    Sigma <- as.matrix(read.csv(file = paste("2000 sigma", sigma, "weak.csv")))
  }
  
  A <- Sigma
  A.standard <- cov2cor(A)
  A.eigen <- eigen(A.standard)
  
  nzero.num.all = NULL;
  FDP.var.all = NULL;
  FDP.all = NULL;
  
  for (i in 1 : n.sim) {
    temp <- proc.time()
    Z.A.raw <- matrix(rnorm(p), 1)
    Z.A <- drop(mu) + A.eigen$vectors %*% diag(sqrt(pmax(A.eigen$values, 0)), p) %*% t(Z.A.raw)
    p.val.raw <- 2 * (1 - pnorm(abs(Z.A)))
    print(proc.time() - temp)
    if (method == "empirical") {
      FDP.all <- c(FDP.all, sum((abs(Z.A) >= qnorm(1-alpha/2)) & (mu == 0))/sum(abs(Z.A) >= qnorm(1-alpha/2)))
    } 
    
    else {
      # Estimating alternative number
      if (method == "true") pi0 <- 1 - nzero.num.true/p
      if (method == "smoother") pi0 <- pi0est(p.val.raw, lambda = seq(0.1, 0.99, 0.01), pi0.method = "smoother")$pi0
      if (method == "bootstrap") pi0 <- pi0est(p.val.raw, lambda = seq(0.1, 0.99, 0.01), pi0.method = "bootstrap")$pi0
      if (method == "Langaas") pi0 <- convest(p.val.raw)
      if (method == "SLIM") pi0 <- SLIMfunc(p.val.raw, STA = 0.1, Divi = 90)$pi0_Est
      
      nzero.num.result <- p - round(p * pi0)
      
      # mu.result <- sum(mu != 0)
      mu.1 <- Z.A
      if (nzero.num.result <= 0) mu.1 <- rep(0, p) else mu.1[abs(mu.1) < (sort(abs(mu.1), decreasing = TRUE)[nzero.num.result])] <- 0
      
      T.cov.1 <- matrix(0, nrow = p, ncol = p)
      
      for (m in 1:p) {
        T.cov.1[m,m] <- (pnorm(qnorm(alpha/2) + mu.1[m]) + pnorm(qnorm(alpha/2) - mu.1[m])) * (1 - pnorm(qnorm(alpha/2) + mu.1[m]) - pnorm(qnorm(alpha/2) - mu.1[m]))
      }
      R0.1 <- sum(pnorm(qnorm(alpha/2)+mu.1) + pnorm(qnorm(alpha/2)-mu.1))
      V0.1 <- sum((pnorm(qnorm(alpha/2)+mu.1) + pnorm((qnorm(alpha/2)-mu.1)))[mu.1 == 0])
      R.var.1 <- sum(diag(T.cov.1))
      V.var.1 <- sum(diag(T.cov.1)[mu.1 == 0])
      
      ####### RCPP fitting
      
      temp_1 = CPP_Eval(p, T.cov.1, A, mu.1, alpha, R.var.1, V.var.1);
      
      R.var.1 = temp_1[[1]];
      V.var.1 = temp_1[[2]];
      RV.cov.1 = temp_1[[3]];
      
      ####################
      
      # print(proc.time() - temp)
      FDP.var.result <- V.var.1/R0.1^2 + R.var.1 * V0.1^2/R0.1^4 - 2 * RV.cov.1 * V0.1/R0.1^3
      
      nzero.num.all = c(nzero.num.all, nzero.num.result);
      FDP.var.all = c(FDP.var.all, FDP.var.result);
    
    # print(FDP.var.result);
    }
  }
  
  if (method == "empirical") {
    return(list(FDP.var = var(FDP.all)));
  } else {
    return(list(nzero.num = nzero.num.all, FDP.var = FDP.var.all));
  }
  
}


### Record the starting time
time_1 = as.numeric(Sys.time());
print(paste("Start at ", Sys.time(), sep = ""))

### Create the data frame to save the results
data.results <- list(sigma=NULL, method=NULL, alpha=NULL, nzero.num.true=NULL, nzero.num=NULL, FDP.var=NULL)

# Specify number of CPU for parallel computing here
# REP.N needs to be a multiple of CPU.N
CPU.N = 1; 

### Set random seeds
input.list = NULL;

for (i in 1 : CPU.N) {
  seed.i = sample(c(1:9999999), 1);
  input.i = list(seed.input = seed.i, CPU.i = i);
  input.list = c(input.list, list(input.i));
}

### Loop through all methods
for (method in c("smoother", "bootstrap", "Langaas", "SLIM")) {
  n.sim = REP.N/CPU.N;
  print(method)

  ### Loop over t and p1
  for (alpha in c(0.005, 0.02, 0.05)) {
      for (nzero.num.true in c(50, 100, 200)) {
      print(c(alpha, nzero.num.true))
      param = list(sigma = sigma, method = method, alpha = alpha, nzero.num.true = nzero.num.true, n.sim = n.sim);

      if(CPU.N == 1) {

        results = evaluation(input.i, param);
        nzero.num = results$nzero.num;
        FDP.var = results$FDP.var;
        data.results$sigma <- c(data.results$sigma, rep(sigma, REP.N))
        data.results$method <- c(data.results$method, rep(method, REP.N))
        data.results$alpha <- c(data.results$alpha, rep(alpha, REP.N))
        data.results$nzero.num.true <- c(data.results$nzero.num.true, rep(nzero.num.true, REP.N))
        data.results$nzero.num <- c(data.results$nzero.num, nzero.num)
        data.results$FDP.var <- c(data.results$FDP.var, FDP.var)

      } else {

        clust = makeCluster(CPU.N, type = "SOCK");
        param = list(sigma = sigma, method = method, alpha = alpha, nzero.num.true = nzero.num.true, n.sim = n.sim);
        temp = clusterApply(clust, input.list, evaluation, param);
        stopCluster(clust);
        nzero.num = NULL;
        FDP.var = NULL;

        for(i in 1 : CPU.N){
          nzero.num = c(nzero.num, temp[[i]]$nzero.num);
          FDP.var = c(FDP.var, temp[[i]]$FDP.var);
        }
        print(nzero.num)
        print(FDP.var)
        data.results$sigma <- c(data.results$sigma, rep(sigma, REP.N))
        data.results$method <- c(data.results$method, rep(method, REP.N))
        data.results$alpha <- c(data.results$alpha, rep(alpha, REP.N))
        data.results$nzero.num.true <- c(data.results$nzero.num.true, rep(nzero.num.true, REP.N))
        data.results$nzero.num <- c(data.results$nzero.num, nzero.num)
        data.results$FDP.var <- c(data.results$FDP.var, FDP.var)
      }

    }
  }
}

CPU.N = 1;

### Set random seeds
input.list = NULL;

for (i in 1 : CPU.N) {
  seed.i = sample(c(1:9999999), 1);
  input.i = list(seed.input = seed.i, CPU.i = i);
  input.list = c(input.list, list(input.i));
}

REP.N <- 1000;

for (method in c("true", "empirical"))
{
  print(method)
  if(method == "true") {
    n.sim = 1;
  } else if(method == "empirical") {
    n.sim = REP.N
  }
  
  for (alpha in c(0.005, 0.02, 0.05)) {
    for (nzero.num.true in c(50, 100, 200)) {
      print(c(alpha, nzero.num.true))
      param = list(sigma = sigma, method = method, alpha = alpha, nzero.num.true = nzero.num.true, n.sim = n.sim);
      FDP.var = evaluation(input.i, param)$FDP.var;
      data.results$sigma <- c(data.results$sigma, sigma)
      data.results$method <- c(data.results$method, method)
      data.results$alpha <- c(data.results$alpha, alpha)
      data.results$nzero.num.true <- c(data.results$nzero.num.true, nzero.num.true)
      data.results$nzero.num <- c(data.results$nzero.num, NA)
      data.results$FDP.var <- c(data.results$FDP.var, FDP.var)
    }
  }
}


write.csv(as.data.frame(data.results), paste("results_", sigma, ".csv", sep = ""))

time_2 = as.numeric(Sys.time());

print(paste("Ends at ", Sys.time(), sep = ""))

print(paste("The computation spends ", round(time_2-time_1), " seconds.", sep = ""))

