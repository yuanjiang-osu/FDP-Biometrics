library(MASS)
library(pbivnorm)
library(qvalue)
library(limma)
library(pbivnorm)
library("quantreg")
source("SLIMcodes.R")

# Specify population group
group <- c("JPTCHB", "CEU", "YRI")[1]

# Specify a method to estimate pi_0
method <- c("Langaas", "SLIM")[1]

# Specify t based on Table 5 for each group
alpha <- 0.005

# read the SNP data
A <- read.table(paste("CCT8_", group, "_012.raw", sep = ""), header = T)
A[is.na(A)] <- 0     # missing SNP measurement being imputed as zero
A <- A[,7:ncol(A)]
A <- as.matrix(A)

n <- nrow(A)
p <- 2 * ncol(A)
index0 <- c(1:p)
X <- matrix(rep(0,n*p),nrow=n)

# derive the recoded design matrix of the SNPs
for (i in 1:n){
  for (j in 1:(p/2)){
    if (A[i,j]==1){
      X[i,(2*j-1)] <- 1
    }
    if (A[i,j]==2){
      X[i,(2*j)] <- 1
    }
  }
}

# exclude redundant SNP data
redun <- rep(0,p)
for (i in 1:p){
  redun[i] <- sum((X[,i]- mean(X[,i]))^2)
}
del <- index0[redun==0]
index1 <- index0[redun!=0]
X <- X[,redun!=0]
p <- length(index1)

# response
y <- log(t(read.table(paste("CCT8_", group, ".txt", sep = ""))))

# fit marginal linear regression
betahat <- rep(0,p)
for (i in 1:p){
  betahat[i] <- sum((X[,i]- mean(X[,i]))*(y-mean(y)))/sum((X[,i]- mean(X[,i]))^2)
}
alphahat <- rep(0,p)
for (i in 1:p){
  alphahat[i] <- mean(y) - betahat[i]*mean(X[,i])
}
resid <- rep(0,p)
for (i in 1:p){
  resid[i] <- sum((y - alphahat[i] - betahat[i]*X[,i])^2)
}
R <- matrix(rep(0,n*p),nrow=n)
for (i in 1:p){
  R[,i] <- (X[,i]- mean(X[,i]))/sqrt(sum((X[,i]- mean(X[,i]))^2))
}

# Principal Factor Approximation procedure
if (group == "JPTCHB") {
  sigma <- 0.0135
} else if (group == "CEU") {
  sigma <- 0.015
} else if (group == "YRI") {
  sigma <- 0.013
}

pca <- svd(t(R))
sqrt_eigval <- pca$d
eigvec <- pca$u
K <- 10
sqrt_lambda <- sqrt_eigval[1:K]
L <- eigvec[,1:K]
# number of factors
for (i in 1:K){
  L[,i] <- L[,i]*sqrt_lambda[i]	   # factor loadings
}

Zhat <- rep(0,p)
for (i in 1:p){
  Zhat[i] <- betahat[i]*sqrt(sum((X[,i]- mean(X[,i]))^2))/sqrt(sigma^2)
  # Z-values
}
z2<-Zhat

P <- 2*(1-pnorm(abs(Zhat)))	 # P-values
o <- order(abs(Zhat))
Zperm <- Zhat[o]
Lperm <- L[o,]
Z.reduce <- Zperm[1:(p*0.95)]
L.reduce <- Lperm[1:(p*0.95),]
W.hat <- rq(Z.reduce ~ L.reduce, 0.5)$coef	# L_1 regression
W.hat <- W.hat[2:(K+1)]	   # least absolute deviation estimates of the factors
LW.est <- L%*%(W.hat)
rs <- rowSums(L^2)

Z_adj <- (Zhat - LW.est)/sqrt(1 - rs)

Sigma_eigen <- eigen(cov2cor(cor(R) - L %*% t(L)))

for (i in 1:20){
  sqrt_lambda.1 <- sqrt(Sigma_eigen$values[1:(n-K-1-i)])
  L.1 <- Sigma_eigen$vectors[, 1:(n-K-1-i)]
  K.1 <- n-K-1-i
  for (j in 1:K.1){
    L.1[,j] <- L.1[,j]*sqrt_lambda.1[j]	   # factor loadings
  }
  P <- 2*(1-pnorm(abs(Z_adj)))	 # P-values
  o <- order(abs(Z_adj))
  Zperm <- Z_adj[o]
  Lperm <- L.1[o,]
  Z.reduce <- Zperm[1:(p*0.95)]
  L.reduce <- Lperm[1:(p*0.95),]
  W.hat <- rq(Z.reduce ~ L.reduce, 0.5)$coef	# L_1 regression
  W.hat <- W.hat[2:(K.1+1)]	   # least absolute deviation estimates of the factors
  LW.est <- L.1%*%(W.hat)
  rs <- rowSums(L.1^2)
  t <- alpha
  Rt <- sum(P <= t)
  a1 <- pnorm((qnorm(t/2) + (LW.est))/sqrt(1 - rs))
  a2 <- pnorm((qnorm(t/2) - (LW.est))/sqrt(1 - rs))
  d <- sum(a1+a2)
  if (Rt > 0){
    FDP <- min(d,Rt)/Rt	
  }
  Vt.hat <- Rt*FDP
  print(i)
  print(FDP)
}

# Estimate pi_0

if (method == "Langaas")
{
  p0.result <- round(p - p * convest(2 * (1 - pnorm(abs(Z_adj)))))
} else if (method == "SLIM") {
  p0.result <- round(p - p * SLIMfunc((2 * (1 - pnorm(abs(Z_adj)))), STA = 0.1, Divi = 90)$pi0_Est)
}

# Calculate the FDP variance under independent assumption
mu <- Z_adj
if (p0.result <= 0) mu <- rep(0, p) else mu[abs(mu) < (sort(abs(mu), decreasing = TRUE)[p0.result])] <- 0
T.cov <- matrix(0, nrow = p, ncol = p)
for (i in 1:p)
  T.cov[i,i] <- (pnorm(qnorm(alpha/2)+mu[i]) + pnorm(qnorm(alpha/2)-mu[i]))*
  (1 - pnorm(qnorm(alpha/2)+mu[i]) - pnorm(qnorm(alpha/2)-mu[i]))
R0 <- sum(pnorm(qnorm(alpha/2)+mu) + pnorm(qnorm(alpha/2)-mu))
V0 <- sum((pnorm(qnorm(alpha/2)+mu) + pnorm((qnorm(alpha/2)-mu)))[mu == 0])
R.var <- sum(diag(T.cov))
V.var <- sum(diag(T.cov)[mu == 0])
RV.cov <- V.var
var.theoretical.ind <- V.var/R0^2 + R.var * V0^2/R0^4 - 2 * RV.cov * V0/R0^3

# Calculate the FDP variance under dependent assumption
Sigma <- cov2cor(cor(R) - L %*% t(L))

Sigma[Sigma >= 1] <- 1
Sigma[Sigma <= -1] <- -1
for (i in 2:p)
  for (j in 1:(i-1)){
    T.cov[i,j] <- 1 - pnorm(qnorm(1-alpha/2) - mu[i]) - pnorm(qnorm(1-alpha/2) - mu[j]) + 
      pbivnorm((qnorm(1-alpha/2) - mu[i]), (qnorm(1-alpha/2) - mu[j]), Sigma[i, j]) + 
      pnorm(qnorm(alpha/2) - mu[i]) - pbivnorm((qnorm(alpha/2) - mu[i]), (qnorm(1-alpha/2) - mu[j]), Sigma[i, j]) + 
      pnorm(qnorm(alpha/2) - mu[j]) - pbivnorm((qnorm(alpha/2) - mu[j]), (qnorm(1-alpha/2) - mu[i]), Sigma[i, j]) +
      pbivnorm((qnorm(alpha/2) - mu[j]), (qnorm(alpha/2) - mu[i]), Sigma[i, j]) - 
      (pnorm(qnorm(alpha/2)+mu[i]) + pnorm(qnorm(alpha/2)-mu[i])) * 
      (pnorm(qnorm(alpha/2)+mu[j]) + pnorm(qnorm(alpha/2)-mu[j]))
    R.var <- R.var + 2 * T.cov[i,j]
    if ((mu[i] == 0) & (mu[j] == 0)){
      V.var <- V.var + 2 * T.cov[i, j]
      RV.cov <- RV.cov + 2 * T.cov[i, j]
    }
    if ((mu[i] * mu[j] == 0) & (mu[i]+mu[j] != 0)) {
      RV.cov <- RV.cov + T.cov[i, j]
    }
  }

var.theoretical <- V.var/R0^2 + R.var * V0^2/R0^4 - 2 * RV.cov * V0/R0^3

# Print the FDP variances under independent and dependent assumptions
print(sqrt(var.theoretical.ind))
print(sqrt(var.theoretical))