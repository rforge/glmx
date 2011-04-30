MoS <- function(object, maxit = 10) {

  ## extract data
  Y <- object$y ## quick&dirty  
  X <- model.matrix(object, model = "mean")
  Z <- model.matrix(object, model = "scale")
  offset <- if(is.null(object$offset)) 0 else object$offset
  dev.resids <- object$family$dev.resids
  wts <- weights(object)
  if(is.null(wts)) wts <- 1

  ## pretend to be in starting values
  par <- object$start
  
  for(i in 1:maxit) {
    ## split up parameters
    beta <- par[1:NCOL(X)]
    gamma <- par[-(1:NCOL(X))]  

    ## process linear predictors, links, etc.
    scale_eta <- object$link$scale$linkfun(1) + drop(Z %*% gamma)
    scale <- object$link$scale$linkinv(scale_eta)
    eta <- drop(X %*% beta + offset) / scale
    fog <- object$link$mean$mu.eta(eta) / scale
    mu <- object$link$mean$linkinv(eta)
    varmu <- object$family$variance(mu)   
 
    ## compute weights
    gbeta <- sqrt(wts) * (1 / varmu) * fog
    ggamma <- - gbeta * eta * object$link$scale$mu.eta(scale_eta)
    Hbeta <- sqrt(wts) * (1 / sqrt(varmu)) * fog
    Hgamma <- - Hbeta * eta * object$link$scale$mu.eta(scale_eta)

    ## scaled cross products  
    A <- crossprod(cbind(Hbeta * X, Hgamma * Z))
    b <- crossprod(cbind(gbeta * X, ggamma * Z), sqrt(wts) * (Y - mu))
    step <- drop(solve(A, b))
    ## should probably have a step halving loop here
    par <- par + step
    dev <- sum(dev.resids(Y,mu,wts))
  }
  
  names(par) <- names(coef(object))  
  return(par)  
}
