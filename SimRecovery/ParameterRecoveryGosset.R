## Recovering nu in Gosset Regression using glmx
## Philipp Doebler - doebler@uni-muenster.de

## Design of Simulation Experiment

## 1) Generate Binary Data according to rbinom(1,1,p)
##    where p = pt(q, nu) and q = sum(beta*x) +  eps, where eps~N(0,sd_eps)
## 2a)nu varies, the number of covariates varies, and their variance as well
##    and of course sample size varies, 
##    beta: intercept is 0, all other betas are 1
##    mean of all covariates is 0
## 2b) sd_deps could vary, but fixed at 0 for now
## 3) 1000 replications
## 4) mean, bias, rmse are calculated

source("glmxhack.R") # for hessian = FALSE
source("plinkshack.R") # would run with plinks.R

nu <- exp(c(0,1,2))
n <- c(500, 1000, 5000)
sd_X <- c(1,2)
type = c("unif") # if type == "unif", then uniform on [-2*sd_X,2*sd_X]
sd_eps <- c(0)
repl = 1000

generate_data <- function(n, nu, sd_X, sd_eps = 1, type = "norm"){
if(type == "norm"){X <-matrix(rnorm(n*length(sd_X),0, sd_X), 
                              ncol = length(sd_X), byrow = TRUE)}
if(type == "unif"){X <-matrix(runif(n*length(sd_X),-2*sd_X, 2*sd_X), 
                              ncol = length(sd_X), byrow = TRUE)}
colnames(X) <- rep("beta", length(sd_X))
 return(list(data = rbinom(n,1,pt(rowSums(X)+rnorm(n, sd = sd_eps),nu)), X = X))
}

gosset_family <- function(nu) binomial(link = gosset(nu))

fit_gosset <- function(data_list){
  glmx.fit(model.matrix(data_list$data~data_list$X), data_list$data, family = gosset_family,
           xlink = "log", xstart = 0, hessian = FALSE)
}

set.seed(1984)

do.sims <- function(nu, n, sd_X, type, sd_eps, repl = 1){
  # only save coefs and  rate of convergence
  output <- expand.grid(nu = nu, n = n, sd_X = sd_X, type = type, sd_eps = sd_eps)
  output$converged <- numeric(nrow(output))
  output$Intercept_mean <- numeric(nrow(output))
  output$beta_mean <- numeric(nrow(output))
  output$lognu_mean <- numeric(nrow(output))
  output$Intercept_bias <- numeric(nrow(output))
  output$beta_bias <- numeric(nrow(output))
  output$lognu_bias <- numeric(nrow(output))
  output$Intercept_rmse <- numeric(nrow(output))
  output$beta_rmse <- numeric(nrow(output))
  output$lognu_rmse <- numeric(nrow(output))

  try_error_counter <- rep(0, nrow(output))
  for(i in 1:nrow(output)){
    cat("Simulating", i, "of", nrow(output), "\n")
    coefs <- numeric(3)
    rmse <- numeric(3)
    n <- output[i,"n"] 
    nu <- output[i, "nu"]
    lognu <- log(nu)
    sd_X <- output[i, "sd_X"]
    sd_eps <- output[i,"sd_eps"]
    type <- output[i, "type"]
    true_coefs <- c(0,1,lognu)
    
    converged <- numeric(repl)
    Z <- matrix(NA, nrow = repl, ncol = 3)
    for(j in 1:repl){
      fff <- numeric()
      class(fff) <- "try-error"
      while(class(fff) == "try-error"){
      fff <- try(fit_gosset(generate_data(n,nu,sd_X, sd_eps,type)))
      if(class(fff) == "try-error"){try_error_counter[i] <- try_error_counter[i] + 1}
      }
      Z[j, ] <- coef(fff)
      converged[j] <- as.numeric(fff$converged)
      if(j %% 10 == 0){cat("j is",j,"\n")}

    }
    output[i,"converged"] <- mean(converged)
    
    if(sum(converged) < 2){
      means <- rmse <- bias <- rep(NA, 3)
    }else{
    means <- colMeans(matrix(Z[converged == 1, ], ncol = 3))
    bias <- colMeans(matrix(Z[converged == 1,], ncol = 3) - matrix(true_coefs, byrow = TRUE, ncol = 3, nrow = sum(converged)))
    rmse <- sqrt(colMeans((matrix(Z[converged == 1,], ncol = 3) - matrix(true_coefs, byrow = TRUE, ncol = 3, nrow = sum(converged)))^2))
   }
    output[i, c("Intercept_mean", "beta_mean", "lognu_mean")] <- means
    output[i, c("Intercept_bias", "beta_bias", "lognu_bias")] <- bias
    output[i, c("Intercept_rmse", "beta_rmse", "lognu_rmse")] <- rmse
    save(output, file = "gosset_sim_incomplete.Robj")
  }
  output <- list(output = output, repl = repl, te = try_error_counter)
  save(output, file = "gosset_sim.Robj")
}

do.sims(nu, n, sd_X, type, sd_eps, repl)
