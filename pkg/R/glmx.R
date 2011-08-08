glmx.fit <- function(x, y, weights = NULL, offset = NULL, start = NULL,
  family = negative.binomial, xlink = "log", xstart = 0,
  profile = TRUE, nuisance = FALSE,           ## FIXME: optionally: extra pars as nuisance pars
  control = glm.control(), xcontrol = list())
{
  ## process link for extra parameters
  if(!inherits(xlink, "link-glm")) xlink <- make.link(xlink)  

  ## starting family
  family_start <- family(xlink$linkinv(xstart))

  ## response
  nobs <- n <- NROW(x)
  if(is.null(weights)) weights <- rep.int(1L, n)
  if(is.null(offset)) offset <- rep.int(0L, n)
  start <- etastart <- mustart <- NULL
  eval(family_start$initialize)

  ## update starting values for coefficients
  start <- glm.fit(x, y, weights = weights, offset = offset,
      start = NULL, control = control, family = family_start)$coefficients

  ## regressors and parameters
  nobs <- sum(weights > 0L)
  k <- NCOL(x)
  q <- length(xstart)
  kstar <- if(family_start$family %in% c("gaussian", "Gamma", "inverse.gaussian")) k + q + 1L else k + q ## FIXME?
  ## FIXME: indicate whether dispersion needs to be estimated or not

  ## objective function
  profile_loglik <- function(par) {
    aic <- glm.fit(x, y, weights = weights, offset = offset,
      start = start, control = control, family = family(xlink$linkinv(par)))$aic
    (aic - kstar) / 2
  }

  full_loglik <- function(par) {
    beta <- par[1:k]
    gamma <- par[-(1:k)]
    f <- family(xlink$linkinv(gamma))
    mu <- f$linkinv(drop(x %*% beta + offset))    
    dev <- sum(f$dev.resids(y, mu, weights))
    (f$aic(y, n, mu, weights, dev) - kstar) / 2
  }

  if(profile) {
    ## optimize profile likelihood first
    opt1 <- optim(par = xstart, fn = profile_loglik, method = "BFGS", hessian = FALSE, control = xcontrol)
    if(opt1$convergence > 0L) warning("optimization of profile likelihood failed to converge")

    ## extract optimal parameters
    xpar <- opt1$par

    ## optimal coefficients at chosen extra parameters
    cf <- glm.fit(x, y, weights = weights, offset = offset,
      start = start, control = control, family = family(xlink$linkinv(xpar)))$coefficients
  } else {
    ## use starting values of extra parameters
    opt1 <- NULL
    xpar <- xstart
    cf <- glm.fit(x, y, weights = weights, offset = offset,
      start = start, control = control, family = family(xlink$linkinv(xpar)))$coefficients
  }
  
  ## optimize full likelihood
  opt2 <- optim(par = c(cf, xpar), fn = full_loglik, method = "BFGS", hessian = TRUE, control = xcontrol)
  if(opt2$convergence > 0L) warning("optimization of full likelihood failed to converge")
    
  ## extract fitted values/parameters
  vc <- solve(as.matrix(opt2$hessian))
  beta <- as.vector(opt2$par[1:k])
  gamma <- as.vector(opt2$par[-(1:k)])
  f <- family(xlink$linkinv(gamma))
  mu <- f$linkinv(drop(x %*% beta + offset))    
  dev <- f$dev.resids(y, mu, weights)

  ## names
  if(!is.null(colnames(x))) {
    names(beta) <- colnames(x)
    gnam <- names(formals(family))[1]
    names(gamma) <- if(length(gamma) > 1L) paste(gnam, 1:length(gamma), sep = "") else gnam
    if(xlink$name != "identity") names(gamma) <- paste(xlink$name, "(", names(gamma), ")", sep = "")
    rownames(vc) <- colnames(vc) <- c(names(beta), names(gamma))
  }

  ## set up return value
  rval <- list(  
    coefficients = c(beta, gamma), ## list(model = beta, extra = gamma),
    residuals = dev,
    fitted.values = mu,
    optim = list(profile = opt1, full = opt2),
    n = n,
    nobs = nobs,
    df = kstar,
    loglik = -opt2$value,
    vcov = vc,
    family = family,
    xlink = xlink,
    converged = opt2$convergence < 1    
  )
  class(rval) <- "glmx"

  return(rval)
}

coef.glmx <- function(object, ...) object$coefficients
vcov.glmx <- function(object, ...) object$vcov
nobs.glmx <- function(object, ...) object$nobs
logLik.glmx <- function(object, ...) structure(object$loglik, df = object$df, class = "logLik")

