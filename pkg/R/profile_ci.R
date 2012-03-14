profile_ci <- function(object, ...)UseMethod("profile_ci")

# only two-sided intervals as of now (i.e. argument type is ignored)
profile_ci.glmx <- function(object, x, y,  linkmin, linkmax, level = 0.95,
                            type = "twosided", weights = NULL, offset = NULL,                 
                            control = glm.control(), xcontrol = list(), start = c(0,0))
{
  if(!object$converged){stop("glmx did not converge")}  
 
  ## regressors and parameters
  n <- nobs <- object$nobs
  if(is.null(weights)) weights <- rep.int(1L, n)
  if(is.null(offset)) offset <- rep.int(0L, n)
  k <- NCOL(x)
  q <- length(object$optim$opt1$par)
  dpar <- object$df - k - q
  beta <- coef(object)[1:k]
  gamma <- coef(object)[-(1:k)]
  if(length(gamma) > 1){stop("only one dimensional parameters")}
  xlink <- object$xlink
  family <- object$family
  
  ## objective function
  profile_loglik <- function(par) {
    aic <- glm.fit(x, y, weights = weights, offset = offset,
      start = start, control = control, family = family(xlink$linkinv(par)))$aic
    -(aic/2 - dpar)
  }

  loglik <- profile_loglik(gamma)
  
  chi2quant <- qchisq(level, 1)
  a <- (1 - level)/2
  a <- c(a, 1 - a)
  pct <- stats:::format.perc(a, 3)
  
  ## define function for which roots need to be found:
  fn <- function(par){2*(profile_loglik(par) - loglik) + chi2quant}

  ## find lower and upper bounds, but first
  ## check if interval includes linkmin or linkmax:
  if(fn(linkmin) > 0){warning("lower bound of interval set equal to linkmin");
                      lower <- linkmin}else{
                      lower <- uniroot(fn, c(linkmin, gamma))$root  
                      }
  if(fn(linkmax) > 0){warning("upper bound of interval set equal to linkmax");
                      upper <- linkmax}else{
                      upper <- uniroot(fn, c(gamma, linkmax))$root  
                      }
  ret <- matrix(c(lower, upper), nrow = 1)
  colnames(ret) <- pct
  rownames(ret) <- names(coef(object))[-(1:k)]
  ret
}



## Test case 1: talpha Regression with
#set.seed(1981)
#alpha <- 1
#x <- runif(10000, -1.2,1.2)
#y <- rbinom(10000,1,talpha(alpha)$linkinv(x))
#x <- model.matrix(y~x)
#
#talpha2_family <- function(alpha)binomial(link = talpha(alpha*2))
#
#fit <- glmx.fit(x,y,xlink = "logit", xstart = 0,  family = talpha2_family)
#pci <- profile_ci.glmx(fit,x,y,linkmin = -5, linkmax = 5)
#pci #wider than the interval based on normal approximation
#
## Test case 2: Gosset Regression with
#nu <- 2
#x <- runif(10000, -1.5,1.5)
#y <- rbinom(10000,1,pt(x,nu))
#x <- model.matrix(y~x)
#
#gosset_family <- function(nu)binomial(link = gosset(nu))
#
#fit <- glmx.fit(x,y,xlink = "log", xstart = 0,  family = gosset_family)
#pci <- profile_ci.glmx(fit,x,y,linkmin = -10, linkmax = 10)
#pci # again wider than the interval based on normal approximation

