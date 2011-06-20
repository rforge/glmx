## some form of generalized logistic link
## constructed as right-censored negative
## binomial distribution with parameter theta

nblogit <- function(theta) {
  structure(list(
    linkfun = function(mu) log(1 - (1 - mu)^(1/theta)) - log(1 - mu)/theta + log(theta),
    linkinv = function(eta) 1 - (1 / (1 + exp(eta)/theta))^theta, 
    mu.eta = function(eta) exp( -(theta + 1) * log(1 + exp(eta)/theta) + eta ),
    valideta = function(eta) TRUE,
    name = "nblogit"), 
    class = "link-glm")
}
