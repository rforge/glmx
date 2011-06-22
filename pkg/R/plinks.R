### angular transformation (also called arcsine transformation)
angular <- function(truncation.warning = FALSE){
	linkfun <- function(mu){
		return(asin(sqrt(mu)))
		}
	
	linkinv <- function(eta){
		etastar <- pmin(asin(1)-.Machine$double.eps, pmax(.Machine$double.eps, eta))
		if(truncation.warning)
			{if(!all(eta == etastar)){
					warning("truncation in inverse of angular link")}
			}		
		return((sin(etastar))^2)
		}
			
	mu.eta <- function(eta){2*cos(eta)*sin(eta)}
	
	valideta <- function(eta) {
		if(truncation.warning){if(!all(eta <= asin(1) & 0 <= eta)){warning("Some of the current etas are out of range")}}
		TRUE
		}
	name <- "angular"
    structure(list(linkfun=linkfun, linkinv=linkinv, 
       mu.eta=mu.eta, valideta=valideta, name=name), 
       class = "link-glm")}
       
arcsine <- angular



### transformations from Aranda-Ordaz (1981)
### ao1 is the symmetric one from Section 2.1

ao1 <- function(phi, truncation.warning = FALSE){
	if(phi < 0){warning("sign of phi ignored, since ao1(phi) = ao1(-phi)!"); phi <- -phi}
	linkfun <- function(mu){
		if(phi == 0){return(log(mu)-log(1-mu))}
		else{return((2/phi)*(mu^phi-(1-mu)^phi)/(mu^phi+(1-mu)^phi))}
		}
	
	linkinv <- function(eta){
		if(phi == 0){return(exp(eta)/(1+exp(eta)))}
		etastar <- pmin(2/phi-.Machine$double.eps, pmax(-2/phi+.Machine$double.eps, eta))
		if(truncation.warning){if(!all(eta == etastar)){warning("truncation in inverse of link")}}
		else{return(((1+0.5*phi*etastar)^(1/phi))/((1+0.5*phi*etastar)^(1/phi)+(1-0.5*phi*etastar)^(1/phi)))}
		}
			
	mu.eta <- function(eta){
		if(phi == 0){return(exp(eta)/((1+exp(eta))^2))}
		else{
			return((1 + 0.5 * phi * eta)^((1/phi) - 1) * (0.5)/((1 + 
    0.5 * phi * eta)^(1/phi) + (1 - 0.5 * phi * eta)^(1/phi)) - 
    ((1 + 0.5 * phi * eta)^(1/phi)) * ((1 + 0.5 * phi * eta)^((1/phi) - 
        1) * (0.5) - (1 - 0.5 * phi * eta)^((1/phi) - 
        1) * (0.5))/((1 + 0.5 * phi * eta)^(1/phi) + 
        (1 - 0.5 * phi * eta)^(1/phi))^2)
			}
		}
	
	valideta <- function(eta){
		if(truncation.warning){
			if(!all(abs(0.5*phi*eta)<1)){warning("Some of the current etas are out of range")}}
		TRUE
		}
	name <- "Aranda-Ordaz symmetric"
    structure(list(linkfun=linkfun, linkinv=linkinv, 
       mu.eta=mu.eta, valideta=valideta, name=name), 
       class = "link-glm")}
    
### transformations from Aranda-Ordaz (1981)
### ao2 is the assymmetric one from 2.2

ao2 <- function(phi, truncation.warning = FALSE){
	linkfun <- function(mu){
		if(phi == 1){return(log(mu)-log(1-mu))}
		if(phi == 0){return(log(-log(1-mu)))}
		else{return(log(((1-mu)^(-phi)-1)/phi))}
		}
		
	linkinv <- function(eta){
		if(phi == 1){return(exp(eta)/(1+exp(eta)))}
		if(phi == 0){return(1-exp(-exp(eta)))}
		else{
			if(phi < 0){	etastar <- pmin(eta, log(-1/phi)-.Machine$double.eps)
							if(truncation.warning){if(!all(etastar == eta)){warning("truncation in inverse of link")}}
							eta <- etastar
						}
			return(1-(1+phi*exp(eta))^(-1/phi))}		}
			
	mu.eta <- function(eta){
		if(phi == 0){return(exp(eta)*exp(-exp(eta)))}
		else{return(exp(eta)*(1+phi*exp(eta))^(-(1+phi)/phi))}}

	valideta <- function(eta) {if(phi == 1 | phi == 0){TRUE}
		else{	if(truncation.warning){if(all(phi*exp(eta) > -1)){warning("Some of the current etas are out of range")}}
				TRUE
			}}
	name <- "Aranda-Ordaz asymmetric"
    structure(list(linkfun=linkfun, linkinv=linkinv, 
       mu.eta=mu.eta, valideta=valideta, name=name), 
       class = "link-glm")}

### folded exponential transformation (Piepho 2003) as link function

foldexp <- function(phi){
	
	upper <- ifelse(phi == 0, 0.5 , (exp(phi)-1)/(2*phi))
	lower <- ifelse(phi == 0, -0.5, (1-exp(phi))/(2*phi))

	linkfun <- function(mu){
		if(phi == 0){return(mu-0.5)}
		return((exp(phi*mu)-exp(phi*(1-mu)))/(2*phi))
		}
	
			
	dfoldexp <- function(mu){
		if(phi == 0){return(rep(1,length(mu)))}
		else{return((exp(phi*mu)+exp(phi*(1-mu)))/2)}
		}
	
	linkinv <- function(eta){
		etastar <- pmin(upper - .Machine$double.eps, pmax(lower+.Machine$double.eps, eta))
		if(phi == 0){return(etastar+0.5)}
		return(0.5 + asinh( (phi*etastar) / exp(phi/2) )/phi)
		}
			
	mu.eta <- function(eta){1/dfoldexp(linkinv(eta))}
	
	valideta <- function(eta) {
		## old version: 	if(all(eta <= upper & eta >= lower)){TRUE}
		## 					FALSE
		TRUE
		}
	name <- "folded exponential"
    structure(list(linkfun=linkfun, linkinv=linkinv, 
       mu.eta=mu.eta, valideta=valideta, name=name), 
       class = "link-glm")}

### Guerrero-Johnson (1982) transformation 

gj <- function(phi){
	linkfun <- function(mu){
		if(phi == 0){log(mu)-log(1-mu)}
		else{return(((mu/(1-mu))^phi-1)/phi)}
		}
	
	linkinv <- function(eta){
		if(phi == 0){exp(eta)/(1+exp(eta))}		
		etastar <- pmax((ifelse(phi > 0, -1, 1)/phi)+.Machine$double.eps, eta)
		return(((phi*etastar+1)^(1/phi))/(1+((phi*etastar+1)^(1/phi))))
		}
			
	mu.eta <- function(eta){
		if(phi == 0){return(exp(eta)/((1+exp(eta))^2))}
		else{return((phi * eta + 1)^((1/phi) - 1)/(1 + ((phi * 
    eta + 1)^(1/phi))) - ((phi * eta + 1)^(1/phi)) * ((phi * 
    eta + 1)^((1/phi) - 1)/(1 + ((phi * eta + 
    1)^(1/phi)))^2))}
		}
	
	valideta <- function(eta){
		TRUE
		## Old version:
		##		if(phi == 0){return(TRUE)}
		##	else{
		##	if(all((ifelse(phi > 0, -1, 1)/phi)<eta)){TRUE}
		##	else{FALSE}
		##	}
			}
	name <- "Guerrero-Johnson"
    structure(list(linkfun=linkfun, linkinv=linkinv, 
       mu.eta=mu.eta, valideta=valideta, name=name), 
       class = "link-glm")}
       




### t_alpha as link-function

t_alpha <- function(alpha, eps = 2*.Machine$double.eps, maxiter = 100){
	if(alpha > 2 | alpha < 0){stop("alpha must be between 0 and 2")}
	linkfun <- function(mu){
			alpha*log(mu)-(2-alpha)*log(1-mu)}
	
			
	dtrafo <- function(alpha,x){
			alpha/x + (2-alpha)/(1-x)}
	
	linkinv <- function(eta){		
		# simple newton algorithm
		# maxiter is rarely reached, only when alpha close to 0 or 2
		# 8 times faster than uniroot
		newton.inv.trafo <- function(alpha,x, eps = 2*.Machine$double.eps, 
									maxiter = 100)
			{
			y <- exp(x)/(1+exp(x))
			n <- 0
			while(n < maxiter){
			y.prev <- y
			y <- pmax(pmin(y - (linkfun(y)-x)/(dtrafo(alpha,y)-1), 1-eps), eps)
			if(max(abs(y.prev - y)) < eps){return(y)}
			n <- n+1
			}
			warning("linkinv: maxiter reached")
			return(y)
			}
				
		if(alpha == 2){
			etastar <- pmin(eta, -.Machine$double.eps)
			return(exp(etastar/2))}
		if(alpha == 0){
			etastar <- pmax(eta, .Machine$double.eps)
			return(1-exp(-etastar/2))}
		return(newton.inv.trafo(alpha,eta, eps = eps, maxiter = maxiter))
	} # end of linkinv
	
	mu.eta <- function(eta){1/dtrafo(alpha,linkinv(eta))}
	valideta <- function(eta) {
## old version:
##		if(alpha == 0){	if(any(eta <= 0)){FALSE}
##						if(any(eta > 0)){TRUE}}
##		if(alpha == 2){	if(any(eta >= 0)){FALSE}
##						if(any(eta < 0)){TRUE}}
		TRUE}

	
	name <- "t_alpha"
	
    structure(list(linkfun=linkfun, linkinv=linkinv, 
       mu.eta=mu.eta, valideta=valideta, name=name), 
       class = "link-glm")}


## Rocke (1993) Beta Transformation as link function

rocke <- function(shape1, shape2 = shape1, truncation.warning = FALSE){
	
	if(shape1 == 0 & shape2 == 0){return(make.link("logit"))}

	if(shape1 < 0 | shape2 < 0){stop("Parameters of Rocke transformation need to both strictly positive or need to both equal 0")}
	

# what to do if shape1 == 0, and shape2 > 0?
# nothing: then the integral does not exist!
		
	 shift <- pbeta(0.5, shape1, shape2)
	 betafactor <- beta(shape1, shape2)
			
	linkfun <- function(mu){
		betafactor*(pbeta(mu,shape1,shape2)-shift)
		}
	
	linkinv <- function(eta){
		etastar <- pmin((1-shift)*betafactor, pmax(eta, -shift*betafactor))	
		if(truncation.warning){if(!all(eta == etastar)){warning("truncation in inverse of Rocke-link")}}	
		pmin(1-.Machine$double.eps, pmax(.Machine$double.eps, qbeta(etastar/betafactor + shift, shape1, shape2)))
		}
	

	mu.eta <- function(eta){
		1/(betafactor*(dbeta(linkinv(eta), shape1, shape2)))
		}

	valideta <- function(eta) {
		if(truncation.warning){
			if(any(	eta < -shift*betafactor | eta > (1-shift)*betafactor )){warning("Some of the current etas are out of range")}
			}
		TRUE
		}
	name <- "Rocke's Beta"
    structure(list(linkfun=linkfun, linkinv=linkinv, 
       mu.eta=mu.eta, valideta=valideta, name=name), 
       class = "link-glm")}

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


Gosset <- function(nu) {
	if(nu <= 0){stop("Gosset link: nu must be positive")}
	if(nu < 0.2){warning("Gosset link might be unreliable for nu < 0.2")}
   qqt <- function(p, nu) 
      sign(p-0.5)*sqrt(qf(1-2*pmin(p,1-p), 1, nu))
   linkfun <- function(mu) qqt(mu,nu)
   linkinv <- function(eta) {
      thresh <- -qqt(.Machine$double.eps,nu)
      eta <- pmin(thresh, pmax(eta, -thresh))
      pt(eta, nu)}
    mu.eta <- function(eta) 
         pmax(dt(eta, nu), .Machine$double.eps)
    valideta <- function(eta) TRUE
    name <- "Gosset"
    structure(list(linkfun=linkfun, linkinv=linkinv, 
       mu.eta=mu.eta, valideta=valideta, name=name), 
       class = "link-glm")}
       
Pregibon <- function(a, b) {
   linkfun <- function(mu) 
      - qPregibon(1 - mu,a = a, b = b)
   linkinv <- function(eta) {
      eps <- .Machine$double.eps^.5
      tlo <- qPregibon(eps, a = a, b = b)
      thi <- qPregibon(1 - eps, a = a, b = b)
      eta <- -pmin(thi, pmax(-eta, tlo))
      1 - pPregibon(-eta, a = a, b = b)}
   mu.eta <- function(eta)
      pmax(dPregibon(-eta, a = a, b = b), 
         .Machine$double.eps^.5)
   valideta <- function(eta) TRUE
   name <- "Pregibon"
   structure(list(linkfun = linkfun, linkinv = linkinv, 
      mu.eta = mu.eta, valideta = valideta, name = name), 
      class = "link-glm")}
qPregibon <- function(x,a = 0,b = 0){
   s <- (qgl(3/4,c(0,1,a-b,a+b)) - 
      qgl(1/4,c(0,1,a-b,a+b)))/2.197224
   qgl(x,c(0,1,a-b,a+b))/s}
pPregibon <- function(x,a = 0,b = 0,tol=1e-12){
    s <- (qgl(3/4,c(0,1,a-b,a+b)) - 
       qgl(1/4,c(0,1,a-b,a+b)))/2.197224
    pgl(x*s, c(0,1,a-b,a+b),inverse.eps=tol)}
dPregibon <- function(x,a = 0,b = 0,tol=1e-12){
    s <- (qgl(3/4,c(0,1,a-b,a+b)) - 
       qgl(1/4,c(0,1,a-b,a+b)))/2.197224
    dgl(x*s, c(0,1,a-b,a+b),inverse.eps=tol)*s}
rPregibon <- function(n,a = 0,b = 0){
    qPregibon(runif(n),a=a,b=b)}
