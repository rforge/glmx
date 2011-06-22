pglm <- function(formula,link="Pregibon", data = NULL, method = "optim", 
		theta = NULL, start = NULL, lower = NULL, upper = NULL, plevel = 0){
  switch(link,
    Gosset = {
      if(is.null(lower)) lower <- .25
      if(is.null(upper)) upper <-  30
      # Caveat Emptor -- no checking for validity of limits!  Fixme!
      g <- function(nu, formula, data)
        glm(formula, family=binomial(link=Gosset(nu)),data=data)$deviance
      h <- function(nu, formula, data, gopt)
        g(nu, formula, data) - gopt -  qchisq(.95,1)
      f <- optimize(g, c(lower,upper), formula = formula, data = data)
      nuhat <- f$minimum
      gopt <- f$objective
      fopt <- glm(formula,family=binomial(link=Gosset(nuhat)),data=data)
      if(h(lower,formula,data,gopt) > 0 )
        nulo <- uniroot(h, c(lower,nuhat), formula = formula, 
		data = data, gopt = gopt)$root
      else
        nulo <-  lower
      if(h(upper,formula,data,gopt) > 0 )
        nuhi <- uniroot(h,c(nuhat,upper),formula = formula, data = data, gopt = gopt)$root
      else
        nuhi <- upper
      return(list(f=fopt, nuhat=nuhat, nulo=nulo, nuhi=nuhi)) }, 
    Pregibon = {
      g <- function(theta,  start = start, formula, data)
        glm(formula, family=binomial(link=Pregibon(a=theta[1],b=theta[2])), 
		start = start, data=data)$deviance
      if(is.null(theta))
	theta <- rep(0,2)
      f0  <- glm(formula, family=binomial(link=Pregibon(a=theta[1],b=theta[2])), 
                start = start, data=data)
      if(f0$converged == FALSE) stop("glm doesn't converge at initial values")
      if(is.null(start)) start <- coef(f0)
      if(is.null(lower)) lower <- theta + c(-2,-2)
      if(is.null(upper))  upper <- theta + c(2,2)
      #Should perhaps check for glm convergence at the corners?
      f <- optim(theta,g,method = "L-BFGS-B", lower = lower, upper = upper, 
			start = start, formula = formula, data = data)
      theta <- f$par
      fn  <- glm(formula, family=binomial(link=Pregibon(a=theta[1],b=theta[2])), 
                start = start, data=data)
      coef <- fn$coef
      mess <- f$message
      gopt <- f$value
      gnull <- glm(formula, family=binomial(link="logit"), data = data)$deviance
      LRtest <- gnull - gopt 
      return(list(theta = theta, f = fn, LRtest = LRtest)) 
    }, stop(paste(link, "link not recognized")))
}
