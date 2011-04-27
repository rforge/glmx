hetglm <- function(formula, data, subset, na.action, weights, offset,
                   family = binomial,
                   link = c("probit", "logit", "cloglog", "cauchit", "log"),
                   link.scale = c("log", "sqrt", "identity"),
 		   control = hetglm.control(...),
		   model = TRUE, y = TRUE, x = FALSE, ...)
{  
  ## call
  cl <- match.call()
  if(missing(data)) data <- environment(formula)
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "na.action", "weights", "offset"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  keep_y <- y
  keep_x <- x
  
  ## links
  if(is.character(link)) link <- match.arg(link)
  if(is.character(link.scale)) link.scale <- match.arg(link.scale)
  link.mean <- link

  ## family
  if(is.character(family)) family <- get(family, mode = "function", envir = parent.frame())
  if(is.function(family)) family <- family(link.mean)
  if(!identical(family$family, "binomial")) stop("currently only 'binomial' family supported")  

  ## formula
  oformula <- as.formula(formula)
  formula <- as.Formula(formula)
  if(length(formula)[2] < 2L) {
    formula <- as.Formula(formula(formula), formula(formula, lhs = 0L))
  } else {
    if(length(formula)[2] > 2L) {
      formula <- Formula(formula(formula, rhs = 1:2))
      warning("formula must not have more than two RHS parts")
    }
  }
  mf$formula <- formula

  ## evaluate model.frame
  mf[[1L]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())

  ## extract terms, model matrix, response
  mt <- terms(formula, data = data)
  mtX <- terms(formula, data = data, rhs = 1L)
  mtZ <- delete.response(terms(formula, data = data, rhs = 2L))
  attr(mtZ, "intercept") <- 1L
  Y <- model.response(mf, "any")
  X <- model.matrix(mtX, mf)
  Z <- model.matrix(mtZ, mf)[, -1, drop = FALSE]

  ## process response
  if(length(dim(Y)) == 1L) {
    nm <- rownames(Y)
    dim(Y) <- NULL
    if(!is.null(nm)) names(Y) <- nm
  }
  n <- NROW(Y)
  if(n < 1) stop("empty model")

  ## weights
  weights <- model.weights(mf)
  if(is.null(weights)) weights <- 1L
  if(length(weights) == 1L) weights <- rep(weights, n)
  weights <- as.vector(weights)
  names(weights) <- rownames(mf)
  
  ## offset  
  offset <- model.offset(mf)
  if(is.null(offset)) offset <- 0L
  if(length(offset) == 1L) offset <- rep.int(offset, n)
  offset <- as.vector(offset)

  ## initialize family (essentially: process response)
  nobs <- n
  y <- Y
  eval(family$initialize)
  Y <- y  

  ## call the actual workhorse: hetglm.fit()
  rval <- hetglm.fit(X, Y, Z, weights, offset, family, link.scale, control)

  ## further model information
  rval$call <- cl
  rval$formula <- oformula
  rval$terms <- list(mean = mtX, scale = mtZ, full = mt)
  rval$levels <- list(mean = .getXlevels(mtX, mf), scale = .getXlevels(mtZ, mf), full = .getXlevels(mt, mf))
  rval$contrasts <- list(mean = attr(X, "contrasts"), scale = attr(Z, "contrasts"))
  if(model) rval$model <- mf
  if(keep_y) rval$y <- Y
  if(keep_x) rval$x <- list(mean = X, scale = Z)

  class(rval) <- "hetglm"  
  return(rval)
}

hetglm.control <- function(method = "BFGS", maxit = 5000, hessian = TRUE, trace = FALSE, start = NULL, ...)
{
  rval <- list(method = method, maxit = maxit, hessian = hessian, trace = trace, start = start)
  rval <- c(rval, list(...))
  if(!is.null(rval$fnscale)) warning("fnscale must not be modified")
  if(!isTRUE(rval$hessian)) warning("hessian must not be modified")
  rval$fnscale <- -1
  if(is.null(rval$reltol)) rval$reltol <- .Machine$double.eps^(1/1.2)
  rval
}

hetglm.fit <- function(x, y, z = NULL, weights = NULL, offset = NULL,
  family = binomial(), link.scale = "log", control = hetglm.control())
{
  ## response and regressor matrix
  nobs <- n <- NROW(x)
  if(is.null(weights)) weights <- rep.int(1L, n)
  if(is.null(offset)) offset <- rep.int(0L, n)
  if(is.null(z)) z <- x
  eval(family$initialize)
  n <- NROW(x)
  k <- NCOL(x)
  m <- NCOL(z)  

  ## link processing
  linkinv <- family$linkinv
  linkobj <- family[c("linkfun", "linkinv", "mu.eta", "valideta", "link")]
  names(linkobj)[5] <- "name"
  class(linkobj) <- "link-glm"
  scale_linkstr <- link.scale
  variance <- family$variance
  mu.eta <- family$mu.eta
  scale_linkobj <- make.link(scale_linkstr)
  scale_linkfun <- scale_linkobj$linkfun
  scale_linkinv <- scale_linkobj$linkinv
  scale_mu.eta <- scale_linkobj$mu.eta

  ## control parameters
  method <- control$method
  start <- control$start
  hessian <- control$hessian
  control <- control[-which(names(control) %in% c("method", "start", "hessian"))]

  ## null model and starting values
  nullreg <- glm.fit(x = x, y = y, weights = weights, offset = offset, family = family)
  if(is.null(start)) {
    start <- list(
      mean = nullreg$coefficients,
      scale = rep.int(0, m)
    )
  }
  if(is.list(start)) start <- do.call("c", start)

  ## objective function and gradient
  loglikfun <- function(par) {
    beta <- par[1:k]
    gamma <- par[-(1:k)]
    eta <- drop(x %*% beta + offset)
    scale_eta <- scale_linkfun(1) + drop(z %*% gamma)
    scale <- scale_linkinv(scale_eta)
    prob <- linkinv(eta / scale)
    ll <- if(NCOL(y) > 1L) {
      dbinom(y[, 1L], size = rowSums(y), prob = prob, log = TRUE)
    } else {
      dbinom(y, size = 1L, prob = prob, log = TRUE)
    }
    sum(weights * ll)
  }
  gradfun <- function(par) { ## hard-coded for binary response
    beta <- par[1:k]
    gamma <- par[-(1:k)]
    scale_eta <- scale_linkfun(1) + drop(z %*% gamma)
    scale <- scale_linkinv(scale_eta)
    eta <- drop(x %*% beta + offset) / scale
    mu <- linkinv(eta)
    resid <- (y - mu) / variance(mu)
    gbeta <- resid * mu.eta(eta) / scale
    ggamma <- - gbeta * eta * scale_mu.eta(scale_eta)
    colSums(cbind(weights * gbeta * x, weights * ggamma * z))
  }

  ## optimize likelihood  
  opt <- optim(par = start, fn = loglikfun, gr = gradfun,
    method = method, hessian = hessian, control = control)
  if(opt$convergence > 0) warning("optimization failed to converge")

  ## extract fitted values/parameters
  vc <- solve(-as.matrix(opt$hessian))
  beta <- as.vector(opt$par[1:k])
  gamma <- as.vector(opt$par[-(1:k)])
  eta <- drop(x %*% beta + offset)
  scale_eta <- scale_linkfun(1) + drop(z %*% gamma)
  scale <- scale_linkinv(scale_eta)
  prob <- linkinv(eta / scale)
  nobs <- sum(weights > 0L)

  ## names
  names(beta) <- colnames(x)
  names(gamma) <- colnames(z)
  rownames(vc) <- colnames(vc) <- c(colnames(x),
    if(m > 0L) paste("(scale)", colnames(z), sep = "_") else NULL)

  ## set up return value
  rval <- list(  
    coefficients = list(mean = beta, scale = gamma),
    residuals = y - prob,
    fitted.values = structure(prob, .Names = names(y)),
    optim = opt,
    method = method,
    control = control,
    start = start,
    weights = if(identical(as.vector(weights), rep(1, n))) NULL else weights,
    offset = if(identical(offset, rep(0, n))) NULL else offset,
    n = n,
    nobs = nobs,
    df.null = nobs - k,
    df.residual = nobs - k - m,
    loglik = opt$value,
    loglik.null = -nullreg$deviance/2,
    vcov = vc,
    family = family,
    link = list(mean = linkobj, scale = scale_linkobj),
    converged = opt$convergence < 1    
  )

  return(rval)
}

print.hetglm <- function(x, digits = max(3, getOption("digits") - 3), ...) 
{
  cat("\nCall:", deparse(x$call, width.cutoff = floor(getOption("width") * 0.85)), "", sep = "\n")

  if(!x$converged) {
    cat("model did not converge\n")
  } else {
    if(length(coef(x))) {
      cat(paste("Coefficients (", x$family$family, " model with ", x$link$mean$name, " link):\n", sep = ""))
      print.default(format(x$coefficients$mean, digits = digits), print.gap = 2, quote = FALSE)
      cat("\n")
      cat(paste("Latent scale model coefficients (with ", x$link$scale$name, " link):\n", sep = ""))
      if(length(x$coefficients$scale)) {
        print.default(format(x$coefficients$scale, digits = digits), print.gap = 2, quote = FALSE)
      } else {
        cat("None (constant scale = 1).")
      }
      cat("\n")
    }
    else cat("No coefficients\n\n")
  }
  
  invisible(x)
}

summary.hetglm <- function(object, vcov. = NULL, type = "deviance", ...)
{
  ## residuals
  type <- match.arg(type, c("deviance", "pearson", "response"))
  object$residuals <- residuals(object, type = type)
  object$residuals.type <- type
  
  ## extend coefficient table
  k <- length(object$coefficients$mean)
  cf <- as.vector(do.call("c", object$coefficients))
  se <- if(is.null(vcov.)) sqrt(diag(object$vcov)) else sqrt(diag(vcov.(object)))
  cf <- cbind(cf, se, cf/se, 2 * pnorm(-abs(cf/se)))
  colnames(cf) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
  cf <- list(mean = cf[1:k, , drop = FALSE], scale = cf[-(1:k), , drop = FALSE])
  rownames(cf$mean) <- names(object$coefficients$mean)
  rownames(cf$scale) <- names(object$coefficients$scale)
  object$coefficients <- cf
  
  ## LR test against null model
  object$lrtest <- if(length(object$coefficients$scale)) {
    c("LR" = -2 * (object$loglik.null - object$loglik),
      "Df" = object$df.null - object$df.resid,
      "p-value" = pchisq(-2 * (object$loglik.null - object$loglik),
        object$df.null - object$df.resid, lower.tail = FALSE))
  } else {
    rep(NA, 3)
  }
  
  ## number of iterations
  object$iterations <- as.vector(tail(na.omit(object$optim$count), 1))  
  
  ## delete some slots
  object$fitted.values <- object$terms <- object$model <- object$y <-
    object$x <- object$levels <- object$contrasts <- object$start <- NULL

  ## return
  class(object) <- "summary.hetglm"
  object
}

print.summary.hetglm <- function(x, digits = max(3, getOption("digits") - 3), ...)
{
  cat("\nCall:", deparse(x$call, width.cutoff = floor(getOption("width") * 0.85)), "", sep = "\n")
  
  if(!x$converged) {
    cat("model did not converge\n")
  } else {
    types <- c("deviance", "pearson", "response")
    Types <- c("Deviance residuals", "Pearson residuals", "Raw response residuals")
    cat(sprintf("%s:\n", Types[types == match.arg(x$residuals.type, types)]))
    print(structure(round(as.vector(quantile(x$residuals)), digits = digits),
      .Names = c("Min", "1Q", "Median", "3Q", "Max")))
  
    cat(paste("\nCoefficients (", x$family$family, " model with ", x$link$mean$name, " link):\n", sep = ""))
    printCoefmat(x$coefficients$mean, digits = digits, signif.legend = FALSE)
  
    cat(paste("\nLatent scale model coefficients (with ", x$link$scale$name, " link):\n", sep = ""))
    if(length(x$coefficients$scale)) {
      printCoefmat(x$coefficients$scale, digits = digits, signif.legend = FALSE)
    } else {
      cat("None (constant scale = 1).\n")
    }
    
    if(getOption("show.signif.stars") & any(do.call("rbind", x$coefficients)[,4] < 0.1))
      cat("---\nSignif. codes: ", "0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1", "\n")
  
    cat("\nLog-likelihood:", formatC(x$loglik, digits = digits),
      "on", sum(sapply(x$coefficients, NROW)), "Df")
    if(!is.na(x$lrtest[1])) cat("\nLR test for homoskedasticity:",
      formatC(x$lrtest[1], digits = digits), "on", x$lrtest[2], "Df, p-value:", format.pval(x$lrtest[3], digits = digits))
    cat(paste("\nNumber of iterations in", x$method, "optimization:", x$iterations, "\n"))
  }
  
  invisible(x)
}

predict.hetglm <- function(object, newdata = NULL,
  type = c("response", "link", "scale"), na.action = na.pass, ...) 
{
  type <- match.arg(type)
  
  if(missing(newdata)) {

    rval <- switch(type,
      "response" = {
        object$fitted.values
      },
      "link" = {
        object$link$mean$linkfun(object$fitted.values)
      },      
      "scale" = {
        gamma <- object$coefficients$scale
        z <- if(is.null(object$x)) model.matrix(object, model = "scale") else object$x$scale
	object$link$scale$linkinv(object$link$scale$linkfun(1) + drop(z %*% gamma))
      }
    )
    names(rval) <- names(object$fitted.values)
    return(rval)

  } else {

    mf <- model.frame(delete.response(object$terms[["full"]]), newdata, na.action = na.action, xlev = object$levels[["full"]])
    X <- model.matrix(delete.response(object$terms$mean), mf, contrasts = object$contrasts$mean)
    Z <- model.matrix(object$terms$scale, mf, contrasts = object$contrasts$scale)[, -1, drop = FALSE]
    offset <- if(!is.null(off.num <- attr(object$terms$full, "offset"))) {
      eval(attr(object$terms$full, "variables")[[off.num + 1]], newdata)
    } else {
      if(!is.null(object$offset)) eval(object$call$offset, newdata)
    }
    if(is.null(offset)) offset <- rep(0, NROW(X))

    eta_mean <- drop(X %*% object$coefficients$mean + offset)
    pred_scale <- object$link$scale$linkinv(object$link$scale$linkfun(1) + drop(Z %*% object$coefficients$scale))

    rval <- switch(type,    
      "response" = {
        object$link$mean$linkinv(eta_mean / pred_scale)
      },      
      "link" = {
        eta_mean / pred_scale
      },      
      "scale" = {
        pred_scale
      }
    )
    names(rval) <- rownames(mf)
    return(rval)

  }
}

coef.hetglm <- function(object, model = c("full", "mean", "scale"), ...) {
  cf <- object$coefficients

  model <-  match.arg(model)  
  switch(model,
    "mean" = {
      cf$mean
    },
    "scale" = {
      cf$scale
    },
    "full" = {
      cf <- c(cf$mean, cf$scale)
      names(cf) <- colnames(object$vcov)
      cf
    }
  )
}

vcov.hetglm <- function(object, model = c("full", "mean", "scale"), ...) {
  vc <- object$vcov
  k <- length(object$coefficients$mean)

  model <- match.arg(model)
  switch(model,
    "mean" = {
      vc[1:k, 1:k, drop = FALSE]
    },
    "scale" = {
      vc <- vc[-(1:k), -(1:k), drop = FALSE]
      colnames(vc) <- rownames(vc) <- names(object$coefficients$scale)
      vc
    },
    "full" = {
      vc
    }
  )
}

bread.hetglm <- function(x, ...) {
  vcov(x) * x$nobs
}

estfun.hetglm <- function(x, ...)
{
  ## extract response y and regressors X and Z
  y <- if(is.null(x$y)) model.response(model.frame(x)) else x$y
  xmat <- if(is.null(x$x)) model.matrix(x, model = "mean") else x$x$mean
  zmat <- if(is.null(x$x)) model.matrix(x, model = "scale") else x$x$scale
  offset <- if(is.null(x$offset)) 0 else x$offset
  wts <- weights(x)
  if(is.null(wts)) wts <- 1

  beta <- x$coefficients$mean
  gamma <- x$coefficients$scale
  scale_eta <- x$link$scale$linkfun(1) + drop(zmat %*% gamma)
  scale <- x$link$scale$linkinv(scale_eta)
  eta <- drop(xmat %*% beta + offset) / scale
  mu <- x$link$mean$linkinv(eta)
  resid <- (y - mu) / x$family$variance(mu)
  gbeta <- resid * x$link$mean$mu.eta(eta) / scale
  ggamma <- - gbeta * eta * x$link$scale$mu.eta(scale_eta)

  rval <- cbind(wts * gbeta * xmat, wts * ggamma * zmat)
  rownames(rval) <- names(y)
  colnames(rval) <- colnames(x$vcov)

  attr(rval, "assign") <- NULL
  return(rval)
}

coeftest.hetglm <- function(x, vcov. = NULL, df = Inf, ...)
  coeftest.default(x, vcov. = vcov., df = df, ...)  

logLik.hetglm <- function(object, ...) {
  structure(object$loglik, df = sum(sapply(object$coefficients, length)), class = "logLik")
}

terms.hetglm <- function(x, model = c("mean", "scale"), ...) {
  x$terms[[match.arg(model)]]
}

model.frame.hetglm <- function(formula, ...) {
  if(!is.null(formula$model)) return(formula$model)
  if(is.Formula(formula$formula)) formula$call$formula <- formula$formula <-
    formula(formula$formula, collapse = TRUE)
  formula$terms <- formula$terms$full
  NextMethod()
}

model.matrix.hetglm <- function(object, model = c("mean", "scale"), ...) {
  model <- match.arg(model)
  rval <- if(!is.null(object$x[[model]])) {
    object$x[[model]]
  } else {
    mm <- model.matrix(object$terms[[model]], model.frame(object), contrasts = object$contrasts[[model]])
    if(model == "scale") mm[, -1, drop = FALSE] else mm
  }
  return(rval)
}

residuals.hetglm <- function(object,
  type = c("deviance", "pearson", "response"), ...)
{
  ## desired type
  type <- match.arg(type)

  ## extract fitted information
  res <- object$residuals
  y <- if(is.null(object$y)) model.response(model.frame(object)) else object$y
  mu <- fitted(object)
  wts <- weights(object)
  if(is.null(wts)) wts <- rep.int(1, length(mu))
  
  res <- switch(type,  
    "deviance" = {
      sign(res) * sqrt(object$family$dev.resids(y, mu, wts))
    },
    "pearson" = {
      sqrt(wts) * res / sqrt(object$family$variance(mu))
    },
    "response" = {
      sqrt(wts) * res
    }
  )
  names(res) <- names(mu)

  return(res)
}

update.hetglm <- function (object, formula., ..., evaluate = TRUE)
{
  call <- object$call
  if(is.null(call)) stop("need an object with call component")
  extras <- match.call(expand.dots = FALSE)$...
  if(!missing(formula.)) call$formula <- formula(update(Formula(formula(object)), formula.))
  if(length(extras)) {
    existing <- !is.na(match(names(extras), names(call)))
    for (a in names(extras)[existing]) call[[a]] <- extras[[a]]
    if(any(!existing)) {
      call <- c(as.list(call), extras[!existing])
      call <- as.call(call)
    }
  }
  if(evaluate) eval(call, parent.frame())
  else call
}
