
Gosset <- function(nu) {
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
