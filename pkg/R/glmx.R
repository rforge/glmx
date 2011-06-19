glmx.fit <- function(x, y, weights = NULL, offset = NULL,
  family = negative.binomial, xlink = "log", xstart = 0,
  control = glm.control(), xcontrol = list()) ## FIXME: maybe merge control arguments?
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

  ## regressors and parameters
  n <- NROW(x)
  k <- NCOL(x)
  q <- length(xstart)
  kstar <- if(family_start$family %in% c("gaussian", "Gamma", "inverse.gaussian")) k + q + 1 else k + q ## FIXME?
  ## FIXME: indicate whether dispersion needs to be estimated or not

  ## objective function
  loglik <- function(par) {
    aic <- glm.fit(x, y, weights = weights, offset = offset, control = control, family = family(xlink$linkinv(par)))$aic
    (aic - kstar) / 2
  }

  ## optimize likelihood  
  opt <- optim(par = xstart, fn = loglik, method = "BFGS", hessian = FALSE, control = xcontrol)
  if(opt$convergence > 0) warning("optimization failed to converge")

  return(opt)
}
