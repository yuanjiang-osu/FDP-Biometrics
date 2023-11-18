### Specify model M1-M7 here
sigma <- c("independence", "equal", "fan", "Cauchy", "3f", "2f", "nonlinear")[1]

### Read the file with the intermediate results
data.results <- read.csv(paste("results_", sigma, ".csv", sep = ""))
data.results$FDP.var[data.results$FDP.var < 0] <- 0

summary.results <- list(method=NULL, alpha=NULL, nzero.num.true=NULL,
                        mean.nzero.num=NULL, sd.nzero.num=NULL,
                        mean.FDP.sd=NULL, sd.FDP.sd=NULL)

### Summarize the simulation results
for (method.1 in c("true", "empirical", "smoother", "bootstrap", "Langaas", "SLIM")) {
  for (alpha.1 in c(0.005, 0.02, 0.05)) {
    for (nzero.num.true.1 in c(50, 100, 200)) {
      results <- subset(data.results, (method == method.1) & (alpha == alpha.1) & (nzero.num.true == nzero.num.true.1))
      # print(mean(results$nzero.num))
      summary.results$method <- c(summary.results$method, method.1)
      summary.results$alpha <- c(summary.results$alpha, alpha.1)
      summary.results$nzero.num.true <- c(summary.results$nzero.num.true, nzero.num.true.1)
      summary.results$mean.nzero.num <- c(summary.results$mean.nzero.num, mean(results$nzero.num))
      summary.results$sd.nzero.num <- c(summary.results$sd.nzero.num, sd(results$nzero.num))
      summary.results$mean.FDP.sd <- c(summary.results$mean.FDP.sd, mean(sqrt(results$FDP.var))*100)
      summary.results$sd.FDP.sd <- c(summary.results$sd.FDP.sd, sd(sqrt(results$FDP.var))*100)
    }
  }
}

print(as.data.frame(summary.results))
